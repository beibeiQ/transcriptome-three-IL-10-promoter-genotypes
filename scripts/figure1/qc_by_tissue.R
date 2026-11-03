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

# ============================================================
# 样本质控：按组织分别绘制 FPKM 箱线图 + 样本聚类树状图
# ============================================================

setwd(project_root())

library(data.table)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)

# ---- 读入数据 ----
fpkm <- fread(project_path("fpkm"))

# ---- 样本重命名 ----
old_names <- names(fpkm)[2:55]

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

new_names <- sapply(old_names, rename_sample)

fpkm_new <- fpkm
names(fpkm_new)[2:55] <- new_names

# ---- 整理为长格式，方便后续按组织拆分 ----
pca_df <- fpkm_new |> dplyr::select(1:55) |> column_to_rownames('gene_id')
pca_df_filter <- pca_df[rowSums(pca_df) > 0, ]

# ---- 按组织拆分 ----
split_tissues <- list(
  I  = pca_df_filter[, grep("_I_", colnames(pca_df_filter))],
  Sp = pca_df_filter[, grep("_Sp_", colnames(pca_df_filter))]
)

tissue_labels <- c(I = "Intestine", Sp = "Spleen")

# ---- 创建输出目录 ----
dir.create("qc_by_tissue", showWarnings = FALSE)

# ============================================================
# 逐组织绘制：箱线图 + 树状图
# ============================================================
for (tissue in names(split_tissues)) {

  mat <- split_tissues[[tissue]]
  tissue_name <- tissue_labels[tissue]

  # --- 样本元信息 ---
  meta <- data.frame(
    Sample   = colnames(mat),
    Genotype = sub("_.*", "", colnames(mat)),
    Day      = sub(".*_([0-9]+)d_.*", "\\1", colnames(mat)),
    Rep      = sub(".*_Rep([0-9]+)", "\\1", colnames(mat)),
    stringsAsFactors = FALSE
  )
  meta$Day <- paste0(meta$Day, "d")
  # 按 基因型+时间点 排序
  meta <- meta |> arrange(Genotype, Day, Rep)
  mat  <- mat[, meta$Sample]

  # --- 构建长格式数据 ---
  long_df <- mat |>
    rownames_to_column("gene_id") |>
    pivot_longer(-gene_id, names_to = "Sample", values_to = "FPKM") |>
    left_join(meta, by = "Sample")
  # 按排序后的 Sample 水平，保证箱线图顺序
  long_df$Sample <- factor(long_df$Sample, levels = meta$Sample)

  # --- 按基因型+时间点着色 ---
  long_df$Group <- paste0(long_df$Genotype, "_", long_df$Day)

  # 箱线图颜色映射
  group_levels <- unique(long_df$Group)
  n_groups     <- length(group_levels)
  # 简单调色板：不同基因型用不同色系，时间点区分深浅
  palette_base <- c(HomoIns = "#E41A1C", HomoDel = "#377EB8", Hetero = "#4DAF4A")
  palette_map  <- setNames(
    rep(palette_base, each = 3)[1:n_groups],
    group_levels
  )

  # ---- 箱线图 ----
  p_box <- ggplot(long_df, aes(x = Sample, y = FPKM, fill = Group)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_y_log10() +
    scale_fill_manual(values = palette_map) +
    theme_bw() +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, size = 8, color = "black"),
      axis.text.y  = element_text(size = 9),
      axis.title    = element_text(size = 11),
      legend.title  = element_text(size = 10),
      legend.text   = element_text(size = 9),
      plot.title    = element_text(hjust = 0.5, size = 13)
    ) +
    labs(
      title = paste0("FPKM Expression Distribution - ", tissue_name),
      y     = "FPKM (log10)",
      x     = "Sample",
      fill  = "Group"
    )

  ggsave(
    filename = file.path("qc_by_tissue", paste0("Boxplot_", tissue, ".pdf")),
    plot = p_box, width = 12, height = 6
  )

  # ---- 树状聚类图 ----
  cor_mat  <- cor(mat)
  dist_mat <- as.dist(1 - cor_mat)
  hc       <- hclust(dist_mat, method = "average")

  # 按聚类顺序给叶节点标签着色
  hc_meta <- data.frame(
    Sample   = colnames(mat),
    Genotype = sub("_.*", "", colnames(mat)),
    stringsAsFactors = FALSE
  )
  group_color <- c(HomoIns = "#E41A1C", HomoDel = "#377EB8", Hetero = "#4DAF4A")

  # 简化标签：HomoIns_I_1d_Rep1 -> I_1d_R1
  short_labels <- colnames(mat) |>
    sub("^HomoIns_", "Ins_", x = _) |>
    sub("^HomoDel_", "Del_", x = _) |>
    sub("_Rep", "_R", x = _) |>
    sub("^Hetero_", "Het_", x = _)

  # 转为 dendrogram 对象，用 dendrapply 给每个叶子标签上色 + 缩短标签
  dend <- as.dendrogram(hc)
  leaf_order   <- labels(dend)
  leaf_genos   <- hc_meta$Genotype[match(leaf_order, hc_meta$Sample)]
  leaf_cols    <- group_color[leaf_genos]
  leaf_shorts  <- short_labels[match(leaf_order, colnames(mat))]

  i <- 0
  dend <- dendrapply(dend, function(n) {
    if (is.leaf(n)) {
      i <<- i + 1
      attr(n, "label")          <- leaf_shorts[i]
      attr(n, "nodePar")        <- list(lab.col = leaf_cols[i])
      attr(n, "nodePar")$cex    <- 0.75
    }
    n
  })

  pdf(file.path("qc_by_tissue", paste0("Dendrogram_", tissue, ".pdf")),
      width = 12, height = 5)
  par(mar = c(5, 4, 3, 1))
  plot(dend,
       main = paste0("Hierarchical Clustering - ", tissue_name),
       xlab = "", sub = "Correlation-based distance",
       ylab = "Height", cex.lab = 0.75)
  rect.hclust(hc, k = 3, border = "grey60")
  legend("topright",
         legend = names(group_color), col = group_color,
         pch = 15, pt.cex = 2, cex = 0.9, bty = "n", title = "Genotype")
  dev.off()

  message("Done: ", tissue_name)
}

message("\n所有文件输出至 qc_by_tissue/ 目录")
