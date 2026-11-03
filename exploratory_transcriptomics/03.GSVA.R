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

source(file.path(project_root(), "load_data.R"))

# kegg数据集
library(data.table)
kegg_df <- fread(annotation_path("kegg_pathway_output.tsv"), sep = "\t")

# GO数据集
library(readr)
go_df <- read_table(annotation_path("go.txt"))



library(GSVA)
library(limma)
library(tidyverse)

# 创建kegg的基因集
kegg_df_1 <- kegg_df[-c(1:10),]
kegg_df_1$Gene_list <- strsplit(kegg_df_1$Gene_id, ",")
# 构建命名 list
kegg_gene_sets <- setNames(kegg_df_1$Gene_list, kegg_df_1$Name)

# 基因表达数据
expr_matrix <- fpkm_new |> 
  dplyr::select(1:55) |> 
  column_to_rownames('gene_id') |> 
  as.matrix()
log_fpkm <- log2(expr_matrix + 1)


# 开始富集分析
gsvaPar_kegg <- gsvaParam(exprData = log_fpkm, geneSets = kegg_gene_sets)
gsva.es.kegg <- gsva(gsvaPar_kegg, verbose=T)
write.csv(gsva.es.kegg,result_path("figure1", "gsva_es_kegg.csv"), row.names = T)

# === 第一步：准备环境和数据 ===
library(limma)

# GSVA matrix: pathways by samples.
# gsva.es.kegg

# 提取样本名
samples <- colnames(gsva.es.kegg)

# 从样本名中提取分组，例如 "Hetero_I_0d"
group <- gsub("_Rep[0-9]+$", "", samples)

# 构建 metadata
meta <- data.frame(sample = samples,
                   group = factor(group))
rownames(meta) <- meta$sample

# === 第二步：构建设计矩阵 ===
design <- model.matrix(~ 0 + group, data = meta)
colnames(design) <- levels(meta$group)

# === 第三步：定义比较（所有组 vs Hetero 组） ===

# 提取所有 group 名称
all_groups <- levels(meta$group)

# 找出 Hetero 相关组
hetero_groups <- grep("^Hetero", all_groups, value = TRUE)

# 找出非 Hetero 的组
non_hetero_groups <- setdiff(all_groups, hetero_groups)

# 构建比较列表：非hetero组 vs 对应 hetero组（时间点和处理相同）
contrast_list <- c()

for (h in hetero_groups) {
  # 构建匹配组名
  base_name <- gsub("^Hetero", "", h)  # 如 "_I_0d"
  other_match <- grep(base_name, non_hetero_groups, value = TRUE)
  
  for (other in other_match) {
    cname <- paste0(h, "_vs_", other)
    contrast_line <- paste0(h, " - ", other)
    contrast_list <- c(contrast_list, contrast_line)
    names(contrast_list)[length(contrast_list)] <- cname
  }
}

# === 第四步：创建对比矩阵 ===
contrast.matrix <- makeContrasts(contrasts = contrast_list, levels = design)

# === 第五步：线性建模与统计 ===
fit <- lmFit(gsva.es.kegg, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# === 第六步：提取所有差异分析结果 ===
gsva_diff_results <- list()

for (i in 1:ncol(contrast.matrix)) {
  comp_name <- colnames(contrast.matrix)[i]
  res <- topTable(fit2, coef = i, number = Inf, adjust = "fdr")
  res$comparison <- comp_name
  gsva_diff_results[[comp_name]] <- res
}

# ✅ 所有结果合并为一个 data.frame
gsva_diff_df <- do.call(rbind, gsva_diff_results)
head(gsva_diff_df)

# === 第七步：筛选显著差异通路（可选） ===
sig_gsva <- subset(gsva_diff_df, adj.P.Val < 0.05 & abs(logFC) > 0.3)

# 保存
write.csv(sig_gsva, file = "GSVA_hetero_vs_others_differential_pathways.csv", row.names = FALSE)



gsva_long <- gsva.es.kegg %>%
  as.data.frame() %>%
  rownames_to_column(var = "Pathway") %>%
  pivot_longer(-Pathway, names_to = "Sample", values_to = "GSVA_score")
gsva_long$group <- gsub("_Rep[0-9]+$", "", gsva_long$Sample)
gsva_long$Days <- str_split(gsva_long$Sample,pattern = '_',n=4,simplify = T)[,3]
gsva_long$Tissue <- str_split(gsva_long$Sample,pattern = '_',n=4,simplify = T)[,2]
gsva_long$Genetype <- str_split(gsva_long$Sample,pattern = '_',n=4,simplify = T)[,1]
gsva_long$Genetype <- factor(gsva_long$Genetype, levels = c('HomoIns','Hetero','HomoDel'))



for (i in seq_along(unique(gsva_long$Pathway))) {
  
  # 当前通路名
  Pathway_id <- unique(gsva_long$Pathway)[i]
  df_plot <- gsva_long |> filter(Pathway == Pathway_id)
  
  # 转换通路名为合法文件名（将 / 替换为 _，还可以替换其他特殊字符）
  safe_name <- Pathway_id |> 
    str_replace_all("[/\\\\:*?\"<>|]", "_") |>   # 替换所有非法字符为 "_"
    str_replace_all(" ", "_")                   # 可选：空格也替换成 "_"

  # 绘图
  p <- ggplot(df_plot, aes(x = Days, y = GSVA_score, fill = Genetype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    theme_bw(base_size = 13) +
    scale_fill_manual(values = c(Hetero = "#d95f02", HomoDel = "#7570b3", HomoIns = "#1b9e77"),
                      name = "Genotype") +
    labs(title = Pathway_id, x = "Days", y = "GSVA enrichment score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ Tissue, nrow = 1)

  # 保存图像
  ggsave(filename = paste0("GSVA_kegg/", safe_name, ".png"),
         plot = p,
         device = "png",
         width = 12,
         height = 8)

  # 打印进度
  message(sprintf("✅ 已完成：%d / %d", i, length(unique(gsva_long$Pathway))))
}



# 1. 提取通路名称（行为显著通路）
sig_pathways <- rownames(sig_gsva)
sig_pathways_1 <- str_split(sig_pathways, pattern = '\\.',n = 2,simplify = T)

# 2. 子集化 GSVA enrichment matrix
gsva_subset <- gsva.es.kegg[sig_pathways_1[,2], ]

# 3. 构建样本注释（基因型 + 时间）
samples <- colnames(gsva_subset)
genotype <- sub("_.*", "", samples)
timepoint <- sub(".*_I_|.*_Sp_", "", samples)
timepoint <- sub("_Rep[0-9]+$", "", timepoint)

annotation_col <- data.frame(
  Genotype = genotype,
  Time = timepoint
)
rownames(annotation_col) <- samples

# 4. 绘图
library(pheatmap)

pheatmap(gsva_subset,
         annotation_col = annotation_col,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row = 8,
         fontsize_col = 6,
         show_rownames = TRUE,
         show_colnames = FALSE)


library(stringr)
library(pheatmap)

# 1. 提取所有比较名称（如：Hetero_I_0d - HomoDel_I_0d）
comparisons <- unique(sig_gsva$comparison)

# 2. 遍历每个比较，画热图
for (comp in comparisons) {
  # 提取该比较的通路名
  comp_rows <- sig_gsva[sig_gsva$comparison == comp, ]
  path_names <- rownames(comp_rows)
  path_names <- str_split(path_names, pattern = '\\.', n = 2, simplify = TRUE)[,2]
  
  # 提取比较组名称（左、右组）
  comp_clean <- gsub(" ", "", comp)
  groups <- str_split(comp_clean, "-", simplify = TRUE)
  g1 <- trimws(groups[1])
  g2 <- trimws(groups[2])
  
  # 获取这两个组的样本列名
  group_samples <- grep(paste0("^(", g1, "|", g2, ")_"), colnames(gsva.es.kegg), value = TRUE)
  
  # 过滤缺失路径
  path_names <- path_names[path_names %in% rownames(gsva.es.kegg)]
  
  # 检查非空
  if (length(path_names) == 0 || length(group_samples) < 2) {
    cat("❗ 跳过空热图：", comp, "\n")
    next
  }
  
  # 提取 GSVA enrichment 子集
  gsva_subset <- gsva.es.kegg[path_names, group_samples]
  
  # 构建注释
  genotype <- sub("_.*", "", group_samples)
  timepoint <- sub(".*_I_|.*_Sp_", "", group_samples)
  timepoint <- sub("_Rep[0-9]+$", "", timepoint)
  
  annotation_col <- data.frame(
    Genotype = genotype,
    Time = timepoint
  )
  rownames(annotation_col) <- group_samples
  
  # 保存为 PDF 或直接画
  pdf(paste0("heatmap_", g1, "_vs_", g2, ".pdf"), width = 8, height = 10)
  pheatmap(gsva_subset,
           annotation_col = annotation_col,
           scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "euclidean",
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           fontsize_row = 8,
           fontsize_col = 6,
           main = comp,
           show_rownames = TRUE,
           show_colnames = TRUE)
  dev.off()
}


grass_carp_hallmark_list <- readRDS(annotation_path("grass_carp_hallmark_list.rds"))

gsvaPar_hallmark <- gsvaParam(exprData = log_fpkm, geneSets = grass_carp_hallmark_list)
gsva.es.hallmark <- gsva(gsvaPar_hallmark, verbose=T)
write.csv(gsva.es.hallmark,result_path("figure1", "gsva.es.hallmark.csv"), row.names = T)


gsva_long <- gsva.es.hallmark %>%
  as.data.frame() %>%
  rownames_to_column(var = "Pathway") %>%
  pivot_longer(-Pathway, names_to = "Sample", values_to = "GSVA_score")
gsva_long$group <- gsub("_Rep[0-9]+$", "", gsva_long$Sample)
gsva_long$Days <- str_split(gsva_long$Sample,pattern = '_',n=4,simplify = T)[,3]
gsva_long$Tissue <- str_split(gsva_long$Sample,pattern = '_',n=4,simplify = T)[,2]
gsva_long$Genetype <- str_split(gsva_long$Sample,pattern = '_',n=4,simplify = T)[,1]
gsva_long$Genetype <- factor(gsva_long$Genetype, levels = c('HomoIns','Hetero','HomoDel'))



for (i in seq_along(unique(gsva_long$Pathway))) {
  
  # 当前通路名
  Pathway_id <- unique(gsva_long$Pathway)[i]
  df_plot <- gsva_long |> filter(Pathway == Pathway_id)
  
  # 转换通路名为合法文件名（将 / 替换为 _，还可以替换其他特殊字符）
  safe_name <- Pathway_id |> 
    str_replace_all("[/\\\\:*?\"<>|]", "_") |>   # 替换所有非法字符为 "_"
    str_replace_all(" ", "_")                   # 可选：空格也替换成 "_"

  # 绘图
  p <- ggplot(df_plot, aes(x = Days, y = GSVA_score, fill = Genetype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    theme_bw(base_size = 13) +
    scale_fill_manual(values = c(Hetero = "#d95f02", HomoDel = "#7570b3", HomoIns = "#1b9e77"),
                      name = "Genotype") +
    labs(title = Pathway_id, x = "Days", y = "GSVA enrichment score") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ Tissue, nrow = 1)

  # 保存图像
  ggsave(filename = paste0("GSVA_hallmark/", safe_name, ".png"),
         plot = p,
         device = "png",
         width = 12,
         height = 8)

  # 打印进度
  message(sprintf("✅ 已完成：%d / %d", i, length(unique(gsva_long$Pathway))))
}


target_df <- kegg_df |> filter(Name == c('Ferroptosis', 'Efferocytosis', 'Autophagy - animal','Lysosome'))
target_gene_id <- str_split(target_df$Gene_id,pattern = ',') |> unlist() |>unique()
for (i in 1:length(target_gene_id)) {
  target_df <- fpkm_plot |> filter(gene_id == target_gene_id[i] &Tissue == 'Sp')
  df_plot <- target_df |> left_join(id_mapping_clean, by = c('gene_id' = 'GeneID'))
  df_plot <- df_plot |> mutate(Gene_name = if_else(is.na(df_plot$GeneName), gene_id,GeneName))
  gene_name <- unique(df_plot$Gene_name)
  df_plot$Genotype <- factor(df_plot$Genotype, levels = c('HomoIns', 'Hetero','HomoDel'))
  
  # 绘图
  p <- ggplot(df_plot, aes(x = Day, y = FPKM, fill = Genotype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    theme_bw(base_size = 13) +
    scale_fill_manual(values = c(Hetero = "#d95f02", HomoDel = "#7570b3", HomoIns = "#1b9e77"),
                      name = "Genotype") +
    labs(title = gene_name, x = "Days", y = "the experssion of fpkm level") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ Tissue, nrow = 1)

  # 保存图像
  ggsave(filename = paste0("GSVA_target/", gene_name, ".png"),
         plot = p,
         device = "png",
         width = 12,
         height = 8)

  # 打印进度
  print(paste0('已完成：',i,'/',length(target_gene_id)))
}

# 创建新的列，每组三列取平均
selected_cols <- names(fpkm_new)[grepl('_Sp_', names(fpkm_new))]
df <- fpkm_new |> dplyr::select(gene_id, all_of(selected_cols)) |>
  column_to_rownames('gene_id')|>as.matrix()


HomoIns_means <- data.frame(
  HomoIns_Sp_0d = rowMeans(df[, grepl("HomoIns_Sp_0d", colnames(df))]),
  HomoIns_Sp_1d = rowMeans(df[, grepl("HomoIns_Sp_1d", colnames(df))]),
  HomoIns_Sp_5d = rowMeans(df[, grepl("HomoIns_Sp_5d", colnames(df))])
) |> as.matrix()


Hetero_means <- data.frame(
  Hetero_Sp_0d = rowMeans(df[, grepl("Hetero_Sp_0d", colnames(df))]),
  Hetero_Sp_1d = rowMeans(df[, grepl("Hetero_Sp_1d", colnames(df))]),
  Hetero_Sp_5d = rowMeans(df[, grepl("Hetero_Sp_5d", colnames(df))])
) |> as.matrix()


HomoDel_means <- data.frame(
  HomoDel_Sp_0d = rowMeans(df[, grepl("HomoDel_Sp_0d", colnames(df))]),
  HomoDel_Sp_1d = rowMeans(df[, grepl("HomoDel_Sp_1d", colnames(df))]),
  HomoDel_Sp_5d = rowMeans(df[, grepl("HomoDel_Sp_5d", colnames(df))])
) |> as.matrix()

df <- cbind(HomoIns_means, Hetero_means, HomoDel_means) |> as.data.frame()

df <- df %>%
  mutate(
    Ratio_0d_Hetero_vs_HomoIns = Hetero_Sp_0d / HomoIns_Sp_0d,
    Ratio_0d_Hetero_vs_HomoDel = Hetero_Sp_0d / HomoDel_Sp_0d,

    Ratio_1d_Hetero_vs_HomoIns = Hetero_Sp_1d / HomoIns_Sp_1d,
    Ratio_1d_Hetero_vs_HomoDel = Hetero_Sp_1d / HomoDel_Sp_1d,

    Ratio_5d_Hetero_vs_HomoIns = Hetero_Sp_5d / HomoIns_Sp_5d,
    Ratio_5d_Hetero_vs_HomoDel = Hetero_Sp_5d / HomoDel_Sp_5d
  )
df = df |> mutate(gene_id = as.character(rownames(df)), .before = HomoIns_Sp_0d)
df_filter <- df |> 
  filter(gene_id %in% target_gene_id)
df_filter <- df_filter |>
  left_join(id_mapping_clean, by = c('gene_id' = 'GeneID'))
# 计算两个比值的乘积（或取最大值用于排序）
df_clean <- df_filter %>%
  mutate(combined_score = Ratio_5d_Hetero_vs_HomoIns + Ratio_5d_Hetero_vs_HomoDel)
