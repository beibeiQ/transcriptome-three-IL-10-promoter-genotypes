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

read_deg_files <- function(path) {
  # 1. 获取所有对比名称
  deg_names <- list.files(path)
  
  # 2. 映射表：基因型和组织名
  genotype_map <- c(A = "HomoIns", B = "HomoDel", C = "Hetero")
  tissue_map   <- c(P = "Spleen", C = "Intestine")
  
  # 3. 批量读取并改名
  deg_list <- lapply(deg_names, function(name) {
    file_path <- file.path(paste0(paste0(project_path("deg_lists"), .Platform$file.sep),
                                  name))
    df <- read.delim(file_path, header = TRUE, sep = "\t")
    
    # 拆解对比名称，例如 AC0 vs BC0
    groups <- unlist(strsplit(name, "vs"))
    
    parse_group <- function(g) {
      geno   <- genotype_map[substr(g, 1, 1)]
      tissue <- tissue_map[substr(g, 2, 2)]   # 取第二位字符
      day    <- substr(g, 3, 3)
      paste0(geno, "_", tissue, "_", day, "d")
    }
    
    name_new <- paste(parse_group(groups[1]), "vs", parse_group(groups[2]), sep = "_")
    
    attr(df, "name") <- name_new  # 也可以存在 names 中
    df <- df |> mutate(description = if_else(
      log2FoldChange > 1 & padj < 0.05, 'Up', if_else(
        log2FoldChange < -1 & padj < 0.05, 'Down', 'Nosign'
      )), .before = gene_chr) 
    df <- df|> arrange(desc(log2FoldChange))
    df <- df |> dplyr::select(-gene_name)
    df <- df |>
      left_join(gene_id, by = 'gene_id') |>
      relocate(gene_name, .after = padj)
    return(df)
  })
  
  # 4. 添加命名
  names(deg_list) <- sapply(deg_list, function(df) attr(df, "name"))
  
  return(deg_list)
}

# 读入差异基因的文件
deg_path <- project_path("deg_lists")
deg_data <- read_deg_files(deg_path)

# 脾脏第一天差异基因数量韦恩图比较
{
  deg_HomoDel <- deg_data$HomoDel_Spleen_1d_vs_HomoDel_Spleen_0d |>
    filter(description != 'Nosign') |> pull(gene_id)
  deg_Hetero <- deg_data$Hetero_Spleen_1d_vs_Hetero_Spleen_0d |>
    filter(description != 'Nosign') |> pull(gene_id)
  deg_HomoIns <- deg_data$HomoIns_Spleen_1d_vs_HomoIns_Spleen_0d |>
    filter(description != 'Nosign') |> pull(gene_id)
  library(ggvenn)
  venn_list <- list(
    HomoDel = deg_HomoDel,
    Hetero  = deg_Hetero,
    HomoIns = deg_HomoIns
  )
  ggvenn(venn_list, 
         fill_color = c("#fc8d59", "#ffffbf", "#91bfdb"), 
         stroke_size = 0.5, 
         set_name_size = 3, 
         text_size = 2) +   # 圈内数字字体大小
    theme(
      text = element_text(family = il10_font_family),
      plot.title = element_text(family = il10_font_family),
      legend.text = element_text(family = il10_font_family)
    )-> p
  ggsave('02.DEGs/deg_sp_1d_veen.pdf',p,device = 'pdf',width = 50,height = 50,units = 'mm')
  
}

# 脾脏第五天差异基因数量韦恩图比较
{
  deg_HomoDel <- deg_data$HomoDel_Spleen_5d_vs_HomoDel_Spleen_1d |>
    filter(description != 'Nosign') |> pull(gene_id)
  deg_Hetero <- deg_data$Hetero_Spleen_5d_vs_Hetero_Spleen_1d |>
    filter(description != 'Nosign') |> pull(gene_id)
  deg_HomoIns <- deg_data$HomoIns_Spleen_5d_vs_HomoIns_Spleen_1d |>
    filter(description != 'Nosign') |> pull(gene_id)
  library(ggvenn)
  venn_list <- list(
    HomoDel = deg_HomoDel,
    Hetero  = deg_Hetero,
    HomoIns = deg_HomoIns
  )
  ggvenn(venn_list, 
         fill_color = c("#fc8d59", "#ffffbf", "#91bfdb"), 
         stroke_size = 0.5, 
         set_name_size = 3, 
         text_size = 2) +   # 圈内数字字体大小
    theme(
      text = element_text(family = il10_font_family),
      plot.title = element_text(family = il10_font_family),
      legend.text = element_text(family = il10_font_family)
    )-> p
  ggsave('02.DEGs/deg_sp_5d_veen.pdf',p,device = 'pdf',width = 50,height = 50,units = 'mm')
  
}

# 肠道第一天差异基因数量韦恩图比较
{
  deg_HomoDel <- deg_data$HomoDel_Intestine_1d_vs_HomoDel_Intestine_0d |>
    filter(description != 'Nosign') |> pull(gene_id)
  deg_Hetero <- deg_data$Hetero_Intestine_1d_vs_Hetero_Intestine_0d |>
    filter(description != 'Nosign') |> pull(gene_id)
  deg_HomoIns <- deg_data$HomoIns_Intestine_1d_vs_HomoIns_Intestine_0d |>
    filter(description != 'Nosign') |> pull(gene_id)
  library(ggvenn)
  venn_list <- list(
    HomoDel = deg_HomoDel,
    Hetero  = deg_Hetero,
    HomoIns = deg_HomoIns
  )
  ggvenn(venn_list, 
         fill_color = c("#fc8d59", "#ffffbf", "#91bfdb"), 
         stroke_size = 0.5, 
         set_name_size = 3, 
         text_size = 2) +   # 圈内数字字体大小
    theme(
      text = element_text(family = il10_font_family),
      plot.title = element_text(family = il10_font_family),
      legend.text = element_text(family = il10_font_family)
    )-> p
  ggsave('02.DEGs/deg_int_1d_veen.pdf',p,device = 'pdf',width = 50,height = 50,units = 'mm')
  
}


# 肠道第五天差异基因数量韦恩图比较
{
  deg_HomoDel <- deg_data$HomoDel_Intestine_5d_vs_HomoDel_Intestine_1d |>
    filter(description != 'Nosign') |> pull(gene_id)
  deg_Hetero <- deg_data$Hetero_Intestine_5d_vs_Hetero_Intestine_1d |>
    filter(description != 'Nosign') |> pull(gene_id)
  deg_HomoIns <- deg_data$HomoIns_Intestine_5d_vs_HomoIns_Intestine_1d |>
    filter(description != 'Nosign') |> pull(gene_id)
  library(ggvenn)
  venn_list <- list(
    HomoDel = deg_HomoDel,
    Hetero  = deg_Hetero,
    HomoIns = deg_HomoIns
  )
  ggvenn(venn_list, 
         fill_color = c("#fc8d59", "#ffffbf", "#91bfdb"), 
         stroke_size = 0.5, 
         set_name_size = 3, 
         text_size = 2) +   # 圈内数字字体大小
    theme(
      text = element_text(family = il10_font_family),
      plot.title = element_text(family = il10_font_family),
      legend.text = element_text(family = il10_font_family)
    )-> p
  ggsave('02.DEGs/deg_int_5d_veen.pdf',p,device = 'pdf',width = 50,height = 50,units = 'mm')
  
}
