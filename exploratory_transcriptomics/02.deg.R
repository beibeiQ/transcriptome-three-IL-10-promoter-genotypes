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
fpkm_new_long <- fpkm_new |> dplyr::select(c(1:55)) |>
  pivot_longer(
  cols = -gene_id,
  names_to = c("Genotype", "Tissue", "Day", "Rep"),
  names_pattern = "(HomoIns|HomoDel|Hetero)_([A-Za-z]+)_([0-9]+)d_Rep([0-9]+)",
  values_to = "FPKM"
)

# geneid转化的数据
# 读取原始文件（每行为一个注释信息）
lines <- readLines("grass_carp_gennid.txt")

# 只保留包含有效字段的行
lines <- lines[grepl("^>lcl\\|NC_", lines)]

# 提取字段
library(stringr)

# 提取 NC 编号
nc_id <- str_extract(lines, "NC_\\d+\\.\\d+")

# 提取 GeneID
gene_id <- str_extract(lines, "GeneID:\\d+")

# 提取 protein_id
protein_id <- str_extract(lines, "protein_id=XP_\\d+\\.\\d+")

# 去掉前缀（只保留数字部分）
gene_id <- sub("GeneID:", "", gene_id)
protein_id <- sub("protein_id=", "", protein_id)

# 组成数据框
id_mapping <- data.frame(NC_ID = nc_id, GeneID = gene_id, ProteinID = protein_id)

library(data.table)
lines <- readLines(annotation_path("gene_name.txt"))
# 提取 NC 编号
nc_id <- str_extract(lines, "NC_\\d+\\.\\d+")
# 提取 GeneID
gene_id <- str_extract(lines, "GeneID:\\d+")
gene_id <- sub("GeneID:", "", gene_id)
parent_vals <- str_extract(lines, "Parent=[^;]+")
parent_clean <- sub("Parent=gene-", "", parent_vals)
id_mapping_1 <- data.frame(NC_ID = nc_id, GeneID = gene_id, GeneName = parent_clean)

emapper <- read.csv(annotation_path("GCF_019924925.1_HZGC01_protein_from_eggnog.faa.emapper.annotations"), skip = 4, sep = '\t') %>%
  dplyr::select(ProteinID=X.query,
                Gene_Symbol=Preferred_name)
emapper$Gene_Symbol <- str_to_lower(emapper$Gene_Symbol)

id_mapping_final <- id_mapping |>
  left_join(id_mapping_1 |> 
              dplyr::select(GeneID, GeneName), 
                                      by = 'GeneID')
id_mapping_final <- id_mapping_final |> distinct(GeneID, .keep_all = T)
grepl('^LOC', id_mapping_final$GeneName) |> table()

# 读入fpkm
library(data.table)
fpkm <- fread(project_path("fpkm"))
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

# 添加基因名
fpkm_new <- fpkm_new |>
  mutate(gene_name_new = str_to_lower(
    str_extract(
      fpkm_new$gene_description, "(?<=GN=)[^\\s&]+")),
         .before = gene_name)

fpkm_new$gene_name_new <- ifelse(is.na(fpkm_new$gene_name_new), 
                             fpkm_new$gene_id, 
                             fpkm_new$gene_name_new)
gene_id <- fpkm_new |> dplyr::select(gene_id, gene_name_new)
names(gene_id)[2] <- 'gene_name'

# 读入差异基因
library(tidyverse)
read_deg_files <- function(path) {
  # 1. 获取所有对比名称
  deg_names <- list.files(path)
  deg_names <- deg_names[-c(37:39)]
  # 2. 映射表：基因型和组织名
  genotype_map <- c(A = "HomoIns", B = "HomoDel", C = "Hetero")
  tissue_map   <- c(P = "Spleen", C = "Intestine")
  
  # 3. 批量读取并改名
  deg_list <- lapply(deg_names, function(name) {
    file_path <- file.path(path, name, paste0(name, "_deg.xls"))
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
names(deg_data)

# 读入原始数据
df_count <- fread(project_path("raw_counts"))
df_count_new <- df_count
names(df_count_new)[2:55] <- new_names
df <- tibble(
  raw = names(df_count),
  clean = names(df_count_new)
)

# 做差异分析
df_vs <- fread('Hetero_vs_Others_Comparisons.txt')
comparisons <- df_vs[, c("Group1", "Group2")] |> as.data.frame()

library(DESeq2)

# 假设 count_data 为原始表达矩阵，meta 为样本信息表
for (i in 1:nrow(comparisons)) {
  group1 <- comparisons[i, 1]
  group2 <- comparisons[i, 2]
  
  message("🔍 比较：", group1, " vs ", group2)

  selected_cols_1 <- grep(paste0("^",group1,"_.*Rep\\d$"), colnames(count_data), value = TRUE)
    selected_cols_2 <- grep(paste0("^",group2,"_.*Rep\\d$"), colnames(count_data), value = TRUE)
  selected_cols <- c(selected_cols_2, selected_cols_1)
  selected_counts <- df_count_new |> dplyr::select("gene_id", selected_cols)
  
  count_matrix <- selected_counts |> 
    as.data.frame() |>
    column_to_rownames("gene_id") |>
    dplyr::select(selected_cols)  # 确保列顺序和 colData 对齐
  
  
  condition <- c(rep('control',3),rep('treatment',3))
  
  colData <- data.frame(
    row.names = selected_cols,
    condition = factor(condition, levels = c('control', 'treatment'))
  )
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = colData, design = ~condition)
  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$pvalue), ]
  df <- as.data.frame(resOrdered)
  df_clean <- df |>
  filter(!is.na(padj))

  df_clean <- df_clean |> mutate(gene_id = rownames(df_clean))
  df_clean <- df_clean |>
  left_join(id_mapping_final |> dplyr::select(GeneID, GeneName), by = c("gene_id" = "GeneID")) |>
  mutate(GeneName = if_else(is.na(GeneName), paste0("LOC", gene_id), GeneName))
  names(df_clean)[8] <- 'gene_name'
  df_clean <- df_clean |>
    mutate(description = if_else(
      log2FoldChange > 1 & padj < 0.05, "Up",
      if_else(log2FoldChange < -1 & padj < 0.05, "Down", "Nosign")
    ))
  df_clean <- df_clean |> arrange(padj, desc(abs(log2FoldChange)))  
  
  # 导出数据
  dir.create("02.DEGs/raw/", showWarnings = FALSE, recursive = TRUE)
  dir.create("02.DEGs/clean/", showWarnings = FALSE, recursive = TRUE)
   write.csv(as.data.frame(df), file = paste0("02.DEGs/raw/DEGs_", group1, "_vs_", group2, ".csv"))
  write.csv(as.data.frame(df_clean), file = paste0("02.DEGs/clean/DEGs_", group1, "_vs_", group2, ".csv"))
}

# 读入文件
path <- "02.DEGs/clean"  # 替换为实际路径
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
deg_list <- lapply(files, read_csv)
names(deg_list) <- basename(files)

plot_gene_expression_boxplot <- function(data, target_gene_id, gene_name = NULL, tissue = NULL) {
  
  # 正确过滤基因
  gene_plot <- data %>%
    filter(gene_id == target_gene_id)

  # 如果传入 tissue，就进一步过滤
   if (!is.null(tissue)) {
     gene_plot <- gene_plot %>% filter(Tissue == tissue)
   }

  # 基因名设置
  if (!is.null(gene_name)) {
     gene_name = gene_name
  } else {
     gene_name = target_gene_id
  }
  
  # 设置因子水平顺序
  gene_plot <- gene_plot %>%
    mutate(Genotype_fac = factor(Genotype, levels = c("HomoIns", "Hetero", "HomoDel")))
    
  p <- ggplot(gene_plot, aes(x = Genotype_fac, y = FPKM, fill = Genotype_fac)) +
        geom_boxplot(outlier.size = 0.5, width = 0.6) +
        labs(
          title = paste("Expression of", gene_name),
          x = "Genotype",
          y = "FPKM"
        ) +
        scale_fill_manual(values = c("HomoIns" = "#619CFF", "HomoDel" = "#00BA38", "Hetero" = "#F8766D")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5))

  if (is.null(tissue)) {
    p1 <- p + facet_wrap(~ Tissue)
  } else {
    p1 <- p + facet_wrap(~ Day, ncol = 3)
  }
  
  return(p1)
}


plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127522978', 
                             gene_name = 'IL10',tissue = 'Sp')
plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127520202', 
                             gene_name = 'IL1b',tissue = 'Sp')
plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127501686', 
                             gene_name = 'IL6',tissue = 'Sp')
plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127500528', 
                             gene_name = 'Tnfa',tissue = 'Sp')

plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127504671', 
                             gene_name = 'stat6',tissue = 'Sp')

# M1细胞maker
plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127512608', 
                             gene_name = 'nos2a',tissue = 'Sp')
plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127494921', 
                             gene_name = 'nos2a',tissue = 'Sp')
plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127502124', 
                             gene_name = 'ptgs2b',tissue = 'Sp')
plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127502124', 
                             gene_name = 'ptgs2b',tissue = 'Sp')

plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127519975', 
                             gene_name = 'stat1a',tissue = 'Sp')
plot_gene_expression_boxplot(fpkm_plot, 
                             target_gene_id = '127519048', 
                             gene_name = 'stat1b',tissue = 'Sp')

# M1细胞maker



run_enrichment_analysis <- function(
  genes,
  name,
  save_path = "enrichment_analysis/",
  species = "mmu",
  keyType = "ENSEMBL",
  OrgDb = org.Mm.eg.db,
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "ENTREZID"),
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  topn_go = 10,
  topn_kegg = 15,
  readable_result = TRUE,
  plot_go = TRUE,
  plot_kegg = TRUE,
  plot_width = 9,
  plot_height = 7,
  do_id_conversion = TRUE
) {
  library(clusterProfiler)
  library(ggplot2)
  library(stringr)

  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
  dir.create(file.path(save_path, "GO"), showWarnings = FALSE)
  dir.create(file.path(save_path, "KEGG"), showWarnings = FALSE)

  # 1. 基因 ID 转换
  if (do_id_conversion) {
    geneID <- tryCatch({
      bitr(genes, fromType = fromType, toType = toType, OrgDb = OrgDb)
    }, error = function(e) {
      warning("ID conversion failed: ", e$message)
      return(NULL)
    })

    if (is.null(geneID) || nrow(geneID) == 0) {
      warning("No valid gene IDs found after conversion, skipping: ", name)
      return(NULL)
    }

    entrez_ids <- geneID$ENTREZID
  } else {
    entrez_ids <- genes
  }

  # 2. GO 富集分析
  message("🟡 [", name, "] 正在进行 GO 富集分析（enrichGO）...")
  ego <- tryCatch({
    enrichGO(
      gene = entrez_ids,
      OrgDb = OrgDb,
      ont = "ALL",
      keyType = keyType,
      pAdjustMethod = "fdr",
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      readable = readable_result
    )
  }, error = function(e) {
    warning("GO 富集出错：", e$message)
    return(NULL)
  })

  if (!is.null(ego) && inherits(ego, "enrichResult") && nrow(ego@result) > 0) {
    df.go.result <- ego@result

    write.table(df.go.result,
                file.path(save_path, "GO", paste0(name, ".go.enrichment.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    if (plot_go) {
      goBP <- subset(df.go.result, ONTOLOGY == "BP")[1:min(topn_go, nrow(subset(df.go.result, ONTOLOGY == "BP"))), ]
      goCC <- subset(df.go.result, ONTOLOGY == "CC")[1:min(topn_go, nrow(subset(df.go.result, ONTOLOGY == "CC"))), ]
      goMF <- subset(df.go.result, ONTOLOGY == "MF")[1:min(topn_go, nrow(subset(df.go.result, ONTOLOGY == "MF"))), ]

      go.df <- rbind(goBP, goCC, goMF)
      go.df <- na.omit(go.df)
      go.df$Description <- factor(go.df$Description, levels = rev(go.df$Description))

      go_bar <- ggplot(go.df, aes(x = Description, y = Count, fill = ONTOLOGY)) +
        geom_bar(stat = "identity", width = 0.9) +
        coord_flip() +
        theme_bw() +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
        labs(x = "GO terms", y = "Gene Number", title = paste0(name, " - Enriched GO Terms")) +
        theme(axis.title = element_text(size = 13),
              axis.text = element_text(size = 11),
              plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
              legend.title = element_text(size = 13),
              legend.text = element_text(size = 11))

      ggsave(go_bar,
             filename = file.path(save_path, "GO", paste0(name, ".go.enrichment.png")),
             width = plot_width, height = plot_height)
    }
  } else {
    message("⚠️ [", name, "] GO富集无结果，跳过保存与绘图。")
    ego <- NULL
  }

  # 3. KEGG 富集分析
  message("🟡 [", name, "] 正在进行 KEGG 富集分析（enrichKEGG）...")
  kegg <- tryCatch({
    enrichKEGG(
      gene = entrez_ids,
      organism = species,
      keyType = "kegg",
      pAdjustMethod = "fdr",
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff
    )
  }, error = function(e) {
    warning("KEGG 富集失败：", e$message)
    return(NULL)
  })

  # 3.1 readable 转换（有的物种 OrgDb 不支持）
  kk <- if (!is.null(kegg) && readable_result) {
    tryCatch({
      setReadable(kegg, OrgDb = OrgDb, keyType = "ENTREZID")
    }, error = function(e) {
      warning("⚠️ KEGG readable 转换失败：", e$message)
      kegg  # fallback
    })
  } else {
    kegg
  }

  if (!is.null(kk) && inherits(kk, "enrichResult") && nrow(kk@result) > 0) {
    df.kegg.result <- kk@result

    write.table(df.kegg.result,
                file.path(save_path, "KEGG", paste0(name, ".kegg.enrichment.txt")),
                sep = "\t", quote = FALSE, row.names = FALSE)

    if (plot_kegg) {
      if ("Description" %in% colnames(df.kegg.result) && all(!is.na(df.kegg.result$Description))) {
        tmp <- str_split_fixed(df.kegg.result$Description, pattern = " - ", n = 2)
        df.kegg.result$Description <- tmp[, 1]
      }

      df.kegg.result <- df.kegg.result[1:min(topn_kegg, nrow(df.kegg.result)), ]
      df.kegg.result <- na.omit(df.kegg.result)

      if ("GeneRatio" %in% colnames(df.kegg.result)) {
        tmp <- str_split(df.kegg.result$GeneRatio, "/", simplify = TRUE)
        generatio <- as.numeric(tmp[, 1]) / as.numeric(tmp[, 2])
        df.kegg.result$GeneRatio <- generatio
      }

      kegg_dot <- ggplot(df.kegg.result, aes(x = GeneRatio, y = reorder(Description, Count))) +
        geom_point(aes(size = Count, color = -log10(p.adjust))) +
        theme_bw() +
        scale_colour_gradient(low = "green", high = "red") +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
        labs(x = "GeneRatio", y = "", title = paste0(name, " - Enriched KEGG Pathways"),
             color = expression(-log10(p.adjust)), size = "Count") +
        theme(axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
              legend.title = element_text(size = 11),
              legend.text = element_text(size = 10))

      ggsave(kegg_dot,
             filename = file.path(save_path, "KEGG", paste0(name, ".kegg.enrichment.png")),
             width = plot_width, height = plot_height)
    }
  } else {
    message("⚠️ [", name, "] KEGG富集无结果，跳过保存与绘图。")
    kk <- NULL
  }

  message("✅ Enrichment analysis completed: ", name)
  return(list(GO = ego, KEGG = kk))
}


library(clusterProfiler)
library(org.CI.eg.db)
enrichment_go_list <- list()
enrichment_kegg_list <- list()

for (i in seq_along(deg_list)) {
  sample_id <- names(deg_list)[i]
  df <- deg_list[[sample_id]]
  df_degs_up <-  df |>
    filter(description == 'Up')
  gene_id <- df_degs_up$gene_id
  
  filename <- str_split_1(string = sample_id, pattern = '\\.')[1]
  res <- run_enrichment_analysis(genes = gene_id,
                          name = filename,
                          save_path = '02.DEGs/',
                          species = 'cide',
                          keyType = 'GID',
                          OrgDb = org.CI.eg.db,
                          readable_result = F,
                          do_id_conversion = F)
  enrichment_go_list[[i]] <- res$GO
  enrichment_kegg_list[[i]] <- res$KEGG
  names(enrichment_go_list)[i] <- filename
  names(enrichment_kegg_list)[i] <- filename
}



id1 <- DEGs_Hetero_Sp_5d_vs_HomoDel_Sp_5d |> 
  filter(description == 'Up') |>
  pull(gene_id)
id2 <- DEGs_Hetero_Sp_5d_vs_HomoIns_Sp_5d |> 
  filter(description == 'Up') |>
  pull(gene_id)
id <- intersect(id1, id2)
