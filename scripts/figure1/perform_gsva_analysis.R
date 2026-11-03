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

# GSVA富集分析函数
perform_gsva_analysis <- function(
    fpkm_data,                    # 输入的表达矩阵数据框
    gene_set,                     # 基因集向量
    gene_set_name = "M",          # 基因集名称
    gene_id_col = "gene_id",      # 基因ID列名
    exclude_tissue = "I",         # 需要排除的组织类型
    plot_title = "GSVA Analysis", # 图表标题
    start_col = 1,                # 表达数据开始列
    end_col = 55,                 # 表达数据结束列
    genotype_levels = c('HomoIns','Hetero','HomoDel'), # 基因型水平顺序
    color_palette = c(Hetero = "#d95f02", HomoDel = "#7570b3", HomoIns = "#1b9e77"), # 颜色设置
    base_size = 13,               # 图表基础字体大小
    return_data = FALSE           # 是否返回数据而不是图表
) {
  
  # 加载必要的包
  required_packages <- c("dplyr", "tibble", "tidyr", "ggplot2", "GSVA", "stringr")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }
  
  # 数据预处理和表达矩阵准备
  tryCatch({
    expr_matrix <- fpkm_data %>% 
      dplyr::select(all_of(c(gene_id_col, names(fpkm_data)[start_col:end_col]))) %>%
      column_to_rownames(gene_id_col) %>% 
      as.matrix()
    
    # Log2转换
    log_fpkm <- log2(expr_matrix + 1)
    
    # 准备基因集
    gene_list <- list(gene_set)
    names(gene_list)[1] <- gene_set_name
    
    # 检查基因集与表达矩阵的交集
    intersect_genes <- intersect(gene_set, rownames(log_fpkm))
    if (length(intersect_genes) == 0) {
      stop("基因集与表达矩阵中的基因没有交集！")
    }
    
    cat("基因集中", length(gene_set), "个基因，与表达矩阵交集为", length(intersect_genes), "个基因\n")
    
    # GSVA分析
    cat("开始GSVA分析...\n")
    gsvaPar <- gsvaParam(exprData = log_fpkm, geneSets = gene_list)
    gsva.es <- gsva(gsvaPar, verbose = TRUE)
    
    # 数据长格式转换
    gsva_long <- gsva.es %>%
      as.data.frame() %>%
      rownames_to_column(var = gene_set_name) %>%
      pivot_longer(-all_of(gene_set_name), names_to = "Sample", values_to = "GSVA_score")
    
    # 提取样本信息
    gsva_long$group <- gsub("_Rep[0-9]+$", "", gsva_long$Sample)
    
    # 使用更灵活的字符串分割
    sample_parts <- str_split(gsva_long$Sample, pattern = '_', simplify = TRUE)
    
    # 检查分割结果的列数
    if (ncol(sample_parts) >= 3) {
      gsva_long$Genetype <- sample_parts[, 1]
      gsva_long$Tissue <- sample_parts[, 2] 
      gsva_long$Days <- sample_parts[, 3]
    } else {
      stop("样本名格式不符合预期，无法正确分割")
    }
    
    # 设置因子水平
    gsva_long$Genetype <- factor(gsva_long$Genetype, levels = genotype_levels)
    
    # 过滤数据
    df_plot <- gsva_long %>% filter(Tissue != exclude_tissue)
    
    # 如果请求返回数据
    if (return_data) {
      return(list(
        gsva_scores = gsva.es,
        long_data = gsva_long,
        plot_data = df_plot
      ))
    }
    
    # 创建图表
    p <- ggplot(df_plot, aes(x = Days, y = GSVA_score, fill = Genetype)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8) +
      theme_bw(base_size = base_size) +
      scale_fill_manual(values = color_palette, name = "Genotype") +
      labs(title = plot_title, x = "Days", y = "GSVA enrichment score") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(~ Genetype, nrow = 1)
    
    return(p)
    
  }, error = function(e) {
    stop(paste("分析过程中出现错误:", e$message))
  })
}

# 使用示例
# 基本用法
# plot <- perform_gsva_analysis(
#   fpkm_data = fpkm_new,
#   gene_set = genset,
#   plot_title = "M03"
# )

# 高级用法 - 自定义参数
# plot <- perform_gsva_analysis(
#   fpkm_data = fpkm_new,
#   gene_set = genset,
#   gene_set_name = "MacrophageSignature",
#   plot_title = "Macrophage GSVA Analysis",
#   exclude_tissue = c("I", "B"),  # 排除多个组织
#   genotype_levels = c('HomoDel', 'Hetero', 'HomoIns'),  # 改变顺序
#   color_palette = c(HomoDel = "red", Hetero = "blue", HomoIns = "green")
# )

# 仅获取数据不绘图
# results <- perform_gsva_analysis(
#   fpkm_data = fpkm_new,
#   gene_set = genset,
#   return_data = TRUE
# )

# 使用不同的列范围
# plot <- perform_gsva_analysis(
#   fpkm_data = fpkm_new,
#   gene_set = genset,
#   start_col = 2,
#   end_col = 60
# )