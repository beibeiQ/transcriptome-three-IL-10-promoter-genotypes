# Project paths.
.il10_start <- Sys.getenv("IL10_PROJECT_ROOT", unset = getwd())
.il10_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!nzchar(Sys.getenv("IL10_PROJECT_ROOT")) && length(.il10_file))
  .il10_start <- dirname(sub("^--file=", "", .il10_file[1]))
.il10_root <- normalizePath(.il10_start, winslash = "/", mustWork = FALSE)
while (!file.exists(file.path(.il10_root, ".il10-project"))) {
  .il10_parent <- dirname(.il10_root)
  if (identical(.il10_parent, .il10_root)) stop("Run from the checkout or set IL10_PROJECT_ROOT.")
  .il10_root <- .il10_parent
}
source(file.path(.il10_root, "config", "paths.R"))

# 数据分析处理步骤
# 设置工作路径
setwd(project_root())

# 文件读入
library(data.table)
fpkm <- fread(project_path("fpkm"))

# 对数据的行名进行替换
library(dplyr)
library(stringr)
library(tidyverse)

# 原列名
old_names <- names(fpkm)[2:55]

# 自定义函数进行重命名
rename_sample <- function(name) {
  geno_code <- str_sub(name, 1, 1)
  tissue_code <- str_sub(name, 2, 2)
  day <- str_extract(name, "\\d(?=_)")
  rep <- str_extract(name, "(?<=_)\\d")
  
  genotype <- case_when(
    geno_code == "A" ~ "HomoIns",
    geno_code == "B" ~ "HomoDel",
    geno_code == "C" ~ "Hetero",
    TRUE ~ "Unknown"
  )
  
  tissue <- case_when(
    tissue_code == "C" ~ "I",
    tissue_code == "P" ~ "Sp",
    TRUE ~ "Unknown"
  )
  
  paste0(genotype, "_", tissue, "_", day, "d_Rep", rep)
}

# 应用
new_names <- sapply(old_names, rename_sample)

# 替换列名
fpkm_new <- fpkm
names(fpkm_new)[2:55] <- new_names


# 创建一个绘图的长数据
fpkm_new_long <- fpkm_new |> dplyr::select(-c(1,57:64)) |>
  pivot_longer(
  cols = -gene_name,
  names_to = c("Genotype", "Tissue", "Day", "Rep"),
  names_pattern = "(HomoIns|HomoDel|Hetero)_([A-Za-z]+)_([0-9]+)d_Rep([0-9]+)",
  values_to = "FPKM"
)

# 准备PCA分析数据
pca_df <- fpkm_new |> dplyr::select(1:55)
pca_df <- pca_df |> column_to_rownames('gene_id')

# 过滤掉0值的数据
pca_df_filter <- pca_df[rowSums(pca_df) > 0,]

library(showtext)
il10_font_family <- register_plot_fonts()
showtext_auto()


fpkm_tmp <- fpkm_new |> dplyr::select(-c(56:64))
fpkm_long <- fpkm_tmp %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "Sample",
    values_to = "FPKM"
  )
p1 <- ggplot(fpkm_long, aes(x = Sample, y = FPKM, fill = Sample)) +
  geom_boxplot(outlier.shape = NA, size = 0.2) +  # 隐藏异常值点
  scale_y_log10() +
  labs(title = "", x = "Sample", y = "FPKM (log10)") +
  theme_bw() +
  theme(text = element_text(family = il10_font_family),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = 'None',
        panel.border = element_rect(color = "black", size = 0.1),   # 细边框
        axis.line = element_line(color = "black", size = 0.1))
p1
ggsave(filename = '01.PCA/FPKM_Distribution_Across_Samples.pdf',
       plot = p1,
       device = 'pdf',
       width = 84.8,
       height = 50.6,
       units = 'mm')

library(pheatmap)
cor <- cor(pca_df_filter)
p2 <- pheatmap(cor)
p2

library(dendextend)

# 计算相关性矩阵
distance_matrix <- as.dist(1 - cor)  # 相关性越高，距离越小
hc <- hclust(distance_matrix, method = "average")  # 可选 complete / ward.D

# 2. 转为 ggplot 支持格式
dendro_data <- ggdendro::dendro_data(hc)
dendro_data$labels$label <- hc$labels  # 覆盖为原始样本名

# 3. 提取标签并识别组织类型
dendro_data$labels$tissue <- ifelse(grepl("_Sp_", dendro_data$labels$label),
                                    "Spleen", "Intestine")
# 4. 画 ggplot 树图（x 为样本名，y 为高度）
library(ggrepel)   # 防止标签重叠（可选）
p3 <- ggplot() +
  geom_segment(data = dendro_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = dendro_data$labels,
            aes(x = x,
                y = -0.05 * max(dendro_data$segments$y),  # 更远离树干
                label = label,
                color = tissue),
            hjust = 1, size = 2.5, family = il10_font_family) +
  scale_color_manual(values = c("Spleen" = "blue", "Intestine" = "red")) +
  theme_minimal(base_family = il10_font_family) +
  theme(
    text = element_text(family = il10_font_family),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.margin = margin(20, 10, 10, 100)  # 增加左边距
  ) +
  labs(title = "",
       x = NULL, y = NULL) +
  coord_flip(clip = "off")  # 允许标签溢出范围显示

ggsave(filename = '01.PCA/sample_cluster.pdf',
       plot = p3,
       device = 'pdf',
       width = 8,
       height = 6,
       dpi = 300,
       units = 'in')


library(ggtree)
library(ape)
# 假设 hc 是 hclust 对象
tree <- as.phylo(hc)
p <- ggtree(tree, layout = "rectangular")
meta <- str_split(p$data$label, "_", n=4,simplify = T) |> as.data.frame()
meta$label <- p$data$label
names(meta) <- c("Genotype", "tissue", "Days", "Rep","lable1")
df <- p$data
p$data <- p$data |>
  left_join(meta, by = c("label" = "lable1"))


# 绘制树状图
p + 
  geom_tiplab(aes(color = tissue), size = 2.5, align = TRUE) +
  scale_color_manual(values = c("Sp" = "blue", "I" = "red")) +
  theme_tree2(base_family = il10_font_family) +
  theme(legend.position = "right")
  

# 创建pcametadata信息
metadata <- names(fpkm_new)[2:55] |> str_split_fixed(pattern = '_', n=4) |> as.data.frame()
rownames(metadata) <- names(fpkm_new)[2:55]
names(metadata) <- c('Genotype', 'Tissue', 'Day', 'Rep')
metadata$Challenge <- ifelse(metadata$Day == '0d', "Mock", "Infected")

# 开始PCA分析
# 先看单因素
library(PCAtools)
library(ggplot2)
library(patchwork)
metadata$Tissue <- factor(metadata$Tissue,
                          levels = c("I", "Sp"),
                          labels = c("Intestine", "Spleen"))

pca <- pca(pca_df_filter, metadata = metadata)
p1 <- biplot(pca,colby = 'Genotype', legendPosition = 'right', encircle = T, title = 'Genotype')

p2 <- biplot(pca,colby = 'Tissue',
             encircle = T, 
             title = '',
             pointSize = 2,
             lab = '')+
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(family = il10_font_family, size = 7))
p2


p3 <- biplot(pca,colby = 'Day', legendPosition = 'right', encircle = T, title = 'Day')
p4 <- biplot(pca,colby = 'Challenge', legendPosition = 'right', encircle = T, title = 'Challenge')

ggsave(filename = '01.PCA/pca_plot_allSample.pdf',
       plot = p2,
       device = 'pdf',
       width = 70,
       height = 65,
       units = 'mm')

(p1 + p2) / (p3 + p4)

eigencorplot(pca, metavars = names(metadata)[-4],posColKey = 'bottom')

df <- metadata[,names(metadata)[-4]]
df$Day <- str_sub(df$Day, 1,1)
df$Day <- as.numeric(df$Day)
df$Genotype <- if_else(df$Genotype == 'HomoIns', 1,if_else(
  df$Genotype == 'HomoDel', 2,3
))
df$Tissue <- if_else(df$Tissue == 'Intestine', 1,2)
df$Challenge <- if_else(df$Challenge == 'Mock', 1,2)

cor_df <- data.frame(PC = character(),
                     Variable = character(),
                     Estimate = numeric(),
                     P_value = numeric(),
                     stringsAsFactors = FALSE)
for (i in 1:10) {
  pca_num <- pca$rotated[[paste0("PC", i)]]
  
  for (l in 1:4) {
    meta_num <- df[[l]]
    
    res <- cor.test(pca_num, meta_num, method = 'pearson')
    est <- res$estimate
    p <- res$p.value
    
    # 添加一行到 cor_df
    cor_df <- rbind(cor_df, 
                    data.frame(PC = paste0("PC", i),
                               Variable = names(df)[l],
                               Estimate = est,
                               P_value = p))
  }
}

library(colorspace)

# 先加上星号 label
cor_df <- cor_df %>%
  mutate(star = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01  ~ "**",
    P_value < 0.05  ~ "*",
    TRUE            ~ ""
  ),
  label = sprintf("%.2f%s", Estimate, star))  # 保留2位小数 + 星号

cor_df$PC <- factor(cor_df$PC, levels = unique(cor_df$PC))
# 画图
ggplot(cor_df, aes(x = PC, y = Variable, fill = Estimate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 1.7) +  # 数字 + 星号
  scale_fill_gradient2(low = "#67a9cf", mid = "white", high = "#ef8a62", 
                       midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  theme_minimal(base_family = il10_font_family) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 6),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.5, 0.05),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 3),
        legend.key.width = unit(0.2, "cm"),        # legend bar 宽
        legend.key.height = unit(0.2, "cm"))+
   guides(fill = guide_colorbar(direction = "horizontal")) -> p1

ggsave(filename = '01.PCA/eigencorplot.pdf',
       plot = p1,
       device = 'pdf',
       width = 70,
       height = 45,
       units = 'mm')

spl_pac_df <- pca_df_filter[,28:54]
spl_metadata <- metadata |> filter(Tissue == 'Sp')
pca_spl <- pca(spl_pac_df, metadata = spl_metadata)
p1 <- biplot(pca_spl, colby = 'Genotype', legendPosition = 'right', encircle = T, title = 'Genotype')
p2 <- biplot(pca_spl, colby = 'Day', legendPosition = 'right', encircle = T, title = 'Day')
p3 <- biplot(pca_spl, colby = 'Challenge', legendPosition = 'right', encircle = T, title = 'Challenge')
(p1 + p2) / (p3 + plot_spacer())


eigencorplot(pca_spl, metavars = names(spl_metadata)[c(-2,-4)])

plotloadings(pca_spl, labSize = 3)

library(ggplot2)

library(ggplot2)
library(dplyr)

plot_gene_expression_boxplot <- function(data, gene_id, gene_name = NULL, tissue = NULL) {
  # 默认用 gene_id 作为标题
  if (is.null(gene_name)) gene_name <- gene_id

  # 子集数据
  gene_plot <- data %>%
    filter(gene_name == gene_id)

  if (!is.null(tissue)) {
    gene_plot <- gene_plot %>%
      filter(Tissue == tissue)
  }

  # 判断是否为空数据
  if (nrow(gene_plot) == 0) {
    return(ggplot() + 
             annotate("text", x = 0.5, y = 0.5, label = paste("No data for", gene_name), size = 6) +
             theme_void())
  }

  # 设置 Genotype 因子顺序
  gene_plot <- gene_plot %>%
    mutate(Genotype = factor(Genotype, levels = c("HomoIns", "HomoDel", "Hetero")))

  # 判断是否可以 facet
  facet_formula <- NULL
  if (length(unique(gene_plot$Tissue)) > 1 | length(unique(gene_plot$Day)) > 1) {
    facet_formula <- ~ Tissue + Day
  }

  # 绘图
  p <- ggplot(gene_plot, aes(x = Genotype, y = FPKM, fill = Genotype)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    scale_y_log10() +
    labs(
      title = paste("Expression of", gene_name),
      x = "Genotype",
      y = "FPKM (log10)"
    ) +
    scale_fill_manual(values = c("HomoIns" = "#619CFF", "HomoDel" = "#00BA38", "Hetero" = "#F8766D")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))

  # 加 facet_wrap 如果需要
  if (!is.null(facet_formula)) {
    p <- p + facet_wrap(facet_formula)
  }

  return(p)
}



p1 <- plot_gene_expression_boxplot(fpkm_new_long, 'LOC127509220', tissue = 'Sp', gene_name = 'fthl27')
p2 <- plot_gene_expression_boxplot(fpkm_new_long, 'LOC127500753', tissue = 'Sp', gene_name = 'progranulin')

p1 / p2

int_pac_df <- pca_df_filter[,1:27]
int_metadata <- metadata |> filter(Tissue == 'I')
pca_int <- pca(int_pac_df, metadata = int_metadata)
p1 <- biplot(pca_int, colby = 'Genotype', legendPosition = 'right', encircle = T,title = 'Genotype')
p2 <- biplot(pca_int, colby = 'Day', legendPosition = 'right', encircle = T,title = 'Day')
p3 <- biplot(pca_int, colby = 'Challenge', legendPosition = 'right', encircle = T,title = 'Challenge')
(p1 + p2) / (p3 + plot_spacer())

eigencorplot(pca_int, metavars = names(int_metadata)[c(-2,-4)])

plotloadings(pca_int, labSize = 3)

genelist =  paste0('LOC',c('127499742','127518377','127494072','127499012','127497684'))
genelist = c(genelist, 'COX3')
names(genelist) = c('Arf4','tmsb4x','aldob','ctsl.1','prss59.1','COX3')
plot_list <- Map(function(gid, gname) {
  plot_gene_expression_boxplot(data = fpkm_new_long, gene_id = gid, tissue = 'I', gene_name = gname)
}, gid = genelist, gname = names(genelist))

# 组合成 2 行 × 3 列
library(cowplot)
combined_plot  <- plot_grid(plotlist = plot_list, ncol = 3, nrow = 2)
legend <- get_legend(plot_list[[1]])
plots_no_legend <- lapply(plot_list, function(p) p + theme(legend.position = "none"))
# 拼所有图（不带图例）
combined <- plot_grid(plotlist = plots_no_legend, ncol = 3)

# 拼主图 + 图例：设置图例比例较小（0.2宽度）
final_plot <- plot_grid(combined, legend, rel_widths = c(1, 0.15), ncol = 2)
final_plot

library(dplyr)
library(tidyr)
pca_sp_0d <- pca(pca_df_filter|> dplyr::select(contains('_Sp_0d_')),
                 metadata = spl_metadata[str_detect(rownames(spl_metadata), "_Sp_0d_"), ],
                 removeVar = 0.1)
pca_sp_1d <- pca(pca_df_filter|> dplyr::select(contains('_Sp_1d_')),
                 metadata = spl_metadata[str_detect(rownames(spl_metadata), "_Sp_1d_"), ],
                 removeVar = 0.1)
pca_sp_5d <- pca(pca_df_filter|> dplyr::select(contains('_Sp_5d_')),
                 metadata = spl_metadata[str_detect(rownames(spl_metadata), "_Sp_5d_"), ],
                 removeVar = 0.1)

p1 <- biplot(pca_sp_0d,
       colby = 'Genotype',
       encircle = TRUE,
       encircleFillKey = c("Hetero" = "#825ca6",
                           "HomoDel" = "#a8dad8",
                           "HomoIns" = "#e7872b"),
       lab = '', pointSize = 2) +
  scale_color_manual(values = c("Hetero" = "#825ca6",
                                "HomoDel" = "#a8dad8",
                                "HomoIns" = "#e7872b")) +
  theme_minimal(base_family = il10_font_family) +
  theme(legend.position = "none",
        text = element_text(family = il10_font_family,size = 7))
p2 <- add_uniform_grid(p1)
ggsave(filename = '01.PCA/new_pca_plot/pca_sp_0d.pdf',
       plot = p2,
       device = 'pdf',
       width = 54,
       height = 41,
       units = 'mm')




p1 <- biplot(pca_sp_1d,
       colby = 'Genotype',
       encircle = TRUE,
       encircleFillKey = c("Hetero" = "#825ca6",
                           "HomoDel" = "#a8dad8",
                           "HomoIns" = "#e7872b"),
       lab = '', pointSize = 2) +
  scale_color_manual(values = c("Hetero" = "#825ca6",
                                "HomoDel" = "#a8dad8",
                                "HomoIns" = "#e7872b")) +
  theme_minimal(base_family = il10_font_family) +
  theme(legend.position = "none",
        text = element_text(family = il10_font_family,size = 7))
p2 <- add_uniform_grid(p1)
ggsave(filename = '01.PCA/new_pca_plot/pca_sp_1d.pdf',
       plot = p2,
       device = 'pdf',
       width = 54,
       height = 41,
       units = 'mm')

p1 <- biplot(pca_sp_5d,
       colby = 'Genotype',
       encircle = TRUE,
       encircleFillKey = c("Hetero" = "#825ca6",
                           "HomoDel" = "#a8dad8",
                           "HomoIns" = "#e7872b"),
       lab = '', pointSize = 2) +
  scale_color_manual(values = c("Hetero" = "#825ca6",
                                "HomoDel" = "#a8dad8",
                                "HomoIns" = "#e7872b")) +
  theme_minimal(base_family = il10_font_family) +
  theme(legend.position = "none",
        text = element_text(family = il10_font_family,size = 7))
p2 <- add_uniform_grid(p1)
ggsave(filename = '01.PCA/new_pca_plot/pca_sp_5d.pdf',
       plot = p2,
       device = 'pdf',
       width = 54,
       height = 41,
       units = 'mm')



pca_I_0d <- pca(pca_df_filter|> dplyr::select(contains('_I_0d_')),
                 metadata = int_metadata[str_detect(rownames(int_metadata), "_I_0d_"), ],
                 removeVar = 0.1)
pca_I_1d <- pca(pca_df_filter|> dplyr::select(contains('_I_1d_')),
                 metadata = int_metadata[str_detect(rownames(int_metadata), "_I_1d_"), ],
                 removeVar = 0.1)
pca_I_5d <- pca(pca_df_filter|> dplyr::select(contains('_I_5d_')),
                 metadata = int_metadata[str_detect(rownames(int_metadata), "_I_5d_"), ],
                 removeVar = 0.1)

p1 <- biplot(pca_I_0d,
       colby = 'Genotype',
       encircle = TRUE,
       encircleFillKey = c("Hetero" = "#825ca6",
                           "HomoDel" = "#a8dad8",
                           "HomoIns" = "#e7872b"), 
       lab = '',pointSize = 2) +
  scale_color_manual(values = c("Hetero" = "#825ca6",
                                "HomoDel" = "#a8dad8",
                                "HomoIns" = "#e7872b")) +
  theme_minimal(base_family = il10_font_family) +
  theme(legend.position = "none",
        text = element_text(family = il10_font_family,size = 7))
p2 <- add_uniform_grid(p1)
ggsave(filename = '01.PCA/new_pca_plot/pca_I_0d.pdf',
       plot = p2,
       device = 'pdf',
       width = 55,
       height = 41,
       units = 'mm')

p1 <- biplot(pca_I_1d,
       colby = 'Genotype',
       encircle = TRUE,
       encircleFillKey = c("Hetero" = "#825ca6",
                           "HomoDel" = "#a8dad8",
                           "HomoIns" = "#e7872b"),
       lab = '',pointSize = 2) +
  scale_color_manual(values = c("Hetero" = "#825ca6",
                                "HomoDel" = "#a8dad8",
                                "HomoIns" = "#e7872b")) +
  theme_minimal(base_family = il10_font_family) +
  theme(legend.position = "none",
        text = element_text(family = il10_font_family,size = 7))
p2 <- add_uniform_grid(p1)
ggsave(filename = '01.PCA/new_pca_plot/pca_I_1d.pdf',
       plot = p2,
       device = 'pdf',
       width = 55,
       height = 41,
       units = 'mm')

p1 <- biplot(pca_I_5d,
       colby = 'Genotype',
       encircle = TRUE,
       encircleFillKey = c("Hetero" = "#825ca6",
                           "HomoDel" = "#a8dad8",
                           "HomoIns" = "#e7872b"),
       lab = '',pointSize = 2) +
  scale_color_manual(values = c("Hetero" = "#825ca6",
                                "HomoDel" = "#a8dad8",
                                "HomoIns" = "#e7872b")) +
  theme_minimal(base_family = il10_font_family) +
  theme(legend.position = "none",
        text = element_text(family = il10_font_family,size = 7))
p2 <- add_uniform_grid(p1)
ggsave(filename = '01.PCA/new_pca_plot/pca_I_5d.pdf',
       plot = p2,
       device = 'pdf',
       width = 55,
       height = 41,
       units = 'mm')

# 绘图网格设置函数
add_uniform_grid <- function(p, n = 6,
                             color = "gray90", size = 0.3,
                             x_range = NULL, y_range = NULL) {
  
  ## 1. 取出 x / y 轴范围（或手动指定）
  gb <- ggplot_build(p)
  xr <- if (is.null(x_range)) gb$layout$panel_params[[1]]$x.range else x_range
  yr <- if (is.null(y_range)) gb$layout$panel_params[[1]]$y.range else y_range
  
  ## 2. 计算网格线位置
  x_breaks <- seq(xr[1], xr[2], length.out = n)
  y_breaks <- seq(yr[1], yr[2], length.out = n)
  
  ## 3. 先做两层网格线（inherit.aes = FALSE 避免继承映射）
  grid_layers <- list(
    geom_vline(xintercept = x_breaks, colour = color, size = size,
               inherit.aes = FALSE),
    geom_hline(yintercept = y_breaks, colour = color, size = size,
               inherit.aes = FALSE)
  )
  
  ## 4. 把网格线插到最前面 → 处于“底层”
  p$layers <- append(grid_layers, p$layers)
  
  ## 5. 关闭默认主/次网格线，返回
  p + theme(panel.grid = element_blank())
}

source('load_data.R')
df_sp <- fpkm_plot |>
  filter(gene_id == '127522978' & Tissue == 'Sp')
df_sp$Genotype <- factor(df_sp$Genotype,
                         levels = c('HomoDel','Hetero','HomoIns'))
df_sp$days <- paste0('Days ',df_sp$Day)
p1 <- ggplot(df_sp, aes(x=Day,y = FPKM,fill = Genotype))+
  geom_boxplot(linewidth = 0.2) +
  scale_fill_manual(values = c("HomoDel" = "#a8dad8" ,
                                "Hetero" = "#825ca6", 
                                "HomoIns" = "#e7872b"))+
  theme_bw()+
  facet_wrap(~Genotype,ncol = 3)+
  scale_y_continuous(limits = c(5, 20), breaks = seq(0, 20, by = 5))+
  theme(text = element_text(family = il10_font_family),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = 'none',
        strip.text = element_text(size = 6),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.1),
        panel.spacing = unit(0.1, "lines"))

trend_df <- df_sp %>% 
  group_by(Genotype, Day) %>% 
  summarise(medianFPKM  = median(FPKM), .groups = "drop")
p2 <- p1 +
  geom_line(
    data   = trend_df,
    aes(x = Day, y = medianFPKM, group = Genotype),
    colour = "black",
    linewidth = 0.3
  ) +
  geom_point(
    data   = trend_df,
    aes(x = Day, y = medianFPKM),
    colour = "black",
    size   = 0.2
  )
p2
ggsave(filename = '01.PCA/il10 expression/sp_il10.pdf',
       plot = p2,
       device = 'pdf',
       width = 44.82,
       height = 33.32,
       units = 'mm')
  
df_int <- fpkm_plot |>
  filter(gene_id == '127522978' & Tissue == 'I')
df_int$Genotype <- factor(df_int$Genotype,
                         levels = c('HomoDel','Hetero','HomoIns'))
df_int$days <- paste0('Days ',df_int$Day)
p1 <- ggplot(df_int, aes(x=Day,y = FPKM,fill = Genotype))+
  geom_boxplot(linewidth = 0.2) +
  scale_fill_manual(values = c("HomoDel" = "#a8dad8" ,
                               "Hetero" = "#825ca6", 
                               "HomoIns" = "#e7872b"))+
  theme_bw()+
  facet_wrap(~Genotype,ncol = 3)+
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, by = 5))+
  theme(text = element_text(family = il10_font_family),
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.position = 'none',
        strip.text = element_text(size = 6),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.grid.minor = element_line(linewidth = 0.1),
        panel.spacing = unit(0.1, "lines"))
trend_df <- df_int %>% 
  group_by(Genotype, Day) %>% 
  summarise(medianFPKM  = median(FPKM), .groups = "drop")
p2 <- p1 +
  geom_line(
    data   = trend_df,
    aes(x = Day, y = medianFPKM, group = Genotype),
    colour = "black",
    linewidth = 0.3
  ) +
  geom_point(
    data   = trend_df,
    aes(x = Day, y = medianFPKM),
    colour = "black",
    size   = 0.2
  )
p2
ggsave(filename = '01.PCA/il10 expression/int_il10.pdf',
       plot = p2,
       device = 'pdf',
       width = 44.82,
       height = 33.32,
       units = 'mm')
