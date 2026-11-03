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

# 读入数据
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
fpkm_plot <- fpkm_new |> dplyr::select(c(1:55)) |>
  pivot_longer(
    cols = -gene_id,
    names_to = c("Genotype", "Tissue", "Day", "Rep"),
    names_pattern = "(HomoIns|HomoDel|Hetero)_([A-Za-z]+)_([0-9]+)d_Rep([0-9]+)",
    values_to = "FPKM"
  )

# 提取字段
library(stringr)
library(data.table)

lines <- readLines(annotation_path("gene_name.txt"))
# 提取 NC 编号
nc_id <- str_extract(lines, "NC_\\d+\\.\\d+")
# 提取 GeneID
gene_id <- str_extract(lines, "GeneID:\\d+")
gene_id <- sub("GeneID:", "", gene_id)
parent_vals <- str_extract(lines, "Parent=[^;]+")
parent_clean <- sub("Parent=gene-", "", parent_vals)
id_mapping_raw <- data.frame(NC_ID = nc_id, GeneID = gene_id, GeneName = parent_clean)

id_mapping_clean <- id_mapping_raw |> distinct(GeneID,.keep_all = T)

