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

# 设置工作路径
setwd(project_root())

# 文件读入
library(data.table)
fpkm <- fread(project_path("fpkm"))

# 对数据的行名进行替换
library(dplyr)
library(stringr)

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
names(fpkm)
names(fpkm_new)
pca_df <- fpkm_new |> dplyr::select(1:55)
pca_df <- pca_df |> column_to_rownames('gene_id')

# 过滤掉0值的数据
pca_df_filter <- pca_df[rowSums(pca_df) > 0,]


# 创建pcametadata信息
metadata <- names(fpkm_new)[2:55] |> str_split_fixed(pattern = '_', n=4) |> as.data.frame()
rownames(metadata) <- names(fpkm_new)[2:55]
names(metadata) <- c('Genotype', 'Tissue', 'Day', 'Rep')
metadata$Challenge <- ifelse(metadata$Day == '0d', "Mock", "Infected")
names(metadata)

# 开始PCA分析
# 先看单因素
library(PCAtools)
pca <- pca(pca_df_filter, metadata = metadata)
biplot(pca,colby = 'Genotype', legendPosition = 'right', encircle = T)
biplot(pca,colby = 'Tissue', legendPosition = 'right', encircle = T)
biplot(pca,colby = 'Day', legendPosition = 'right', encircle = T)
biplot(pca,colby = 'Rep', legendPosition = 'right', encircle = T)
biplot(pca,colby = 'Challenge', legendPosition = 'right', encircle = T)
biplot(pca,x = 'PC1',y = 'PC10',colby = 'Genotype', legendPosition = 'right',encircle = T)
eigencorplot(pca, metavars = names(metadata))
pairsplot(pca, 
          components = getComponents(pca, c(1:10)), 
          triangle = TRUE, 
          trianglelabSize = 12, 
          hline = 0, vline = 0, 
          pointSize = 0.4, 
          gridlines.major = FALSE, 
          gridlines.minor = FALSE, 
          colby = 'Genotype', 
          title = 'Pairs plot', plotaxes = FALSE, 
          margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))

# 查看基因表达水平
library(ggplot2)
# 将数据框转为长格式
pca_df_filter_1 <- pca_df_filter
pca_df_filter_1$gene_id <- rownames(pca_df_filter)
fpkm_long <- melt(pca_df_filter_1, variable.name = "Sample", value.name = "FPKM")

# 画箱线图
ggplot(fpkm_long, aes(x = Sample, y = FPKM)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_y_log10() +  # Log10 FPKM scale.
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "FPKM Expression Distribution per Sample",
       y = "FPKM (log10)", x = "Sample")



# 热图
library(pheatmap)
cor <- cor(pca_df_filter)
p <- pheatmap(cor)


# 树状聚类
distance_matrix <- as.dist(1 - cor)  # 相关性越高，距离越小
hc <- hclust(distance_matrix, method = "average")  # 可选 complete / ward.D
plot(hc, main = "Hierarchical Clustering of Samples",
     xlab = "", sub = "", cex = 0.8, hang = -1)


######    脾脏-pca      ######
names(pca_df_filter)
names(metadata)
spl_pac_df <- pca_df_filter[,28:54]
spl_metadata <- metadata |> filter(Tissue == 'Sp')
pca_spl <- pca(spl_pac_df, metadata = spl_metadata)
biplot(pca_spl, colby = 'Genotype', legendPosition = 'right', encircle = T)
screeplot(pca_spl, axisLabSize = 18, titleLabSize = 22)
biplot(pca_spl, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
biplot(pca_spl, colby = 'Day', legendPosition = 'right', encircle = T)
biplot(pca_spl, colby = 'Challenge', legendPosition = 'right', encircle = T)


######    肠道-pca      ######
int_pac_df <- pca_df_filter[,1:27]
int_metadata <- metadata |> filter(Tissue == 'I')
pca_int <- pca(int_pac_df, metadata = int_metadata)
biplot(pca_int, colby = 'Genotype', legendPosition = 'right', encircle = T)
biplot(pca_int, colby = 'Day', legendPosition = 'right', encircle = T)
biplot(pca_int, colby = 'Challenge', legendPosition = 'right', encircle = T)


fpkm |> filter(gene_id == '127509220')
fpkm |> filter(gene_id == '127500753')
fpkm |> filter(gene_name %in% 'il')
fpkm_new |>
  filter(str_detect(gene_name, "il"))
fpkm_new[,'gene_name'] |> head()
