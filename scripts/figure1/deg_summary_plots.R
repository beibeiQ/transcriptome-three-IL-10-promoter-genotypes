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

source('load_data.R')
library(showtext)
il10_font_family <- register_plot_fonts()
showtext_auto()

spl_file <- list.files(paste0(project_path("deg_lists"), .Platform$file.sep),"^.P[^v]*vs")
library(stringr)
library(xlsx)
res_df <- data.frame('Genotype'= NA, 
                 'compare'= NA,
                 'Up' = NA, 
                 'down' = NA)
for (i in seq_along(spl_file)) {
  file <- spl_file[i]
  Genotype <- if_else(str_sub(file,0,1) == 'A', 
                        'HomoIns', 
                        if_else(
                          str_sub(file,0,1) == 'B',
                          'HomoDel','Hetero'
                          ))
  treat <- str_sub(file,3,3)
  con <- str_sub(file,8,8)
  compare <- paste0(treat,'_vs_',con)
  # 读入数据
  df <- fread(paste0(paste0(project_path("deg_lists"), .Platform$file.sep),file))
  df <- df |> mutate(level = if_else(
    log2FoldChange > 1 & padj < 0.05,'Up',if_else(
      log2FoldChange < -1 & padj < 0.05,'Down','Nosign'
    )
  ))
  Up <- count(df, level) |> filter(level == 'Down') |> pull(n)
  Down <- count(df, level) |> filter(level == 'Up') |> pull(n)
  res_df[i,1] <- Genotype
  res_df[i,2] <- compare
  res_df[i,3] <- Up
  res_df[i,4] <- Down
  print(paste0('已完成：', i))
}

library(ggplot2)
library(dplyr)
library(tidyr)

# 构造标签列（推荐用于y轴）
df_deg <- res_df %>%
  mutate(group_label = paste(Genotype, compare, sep = "_"))

# 转为长格式，Down 为负数
df_long <- df_deg %>%
  dplyr::select(group_label, Up, down) %>%
  pivot_longer(cols = c("Up", "down"), names_to = "Direction", values_to = "Count") %>%
  mutate(Count = ifelse(Direction == "down", -Count, Count),
         Direction = ifelse(Direction == "down", "Down", "Up"))

df_long$com <- factor(str_split(df_long$group_label, pattern = '_',n = 2, simplify = T) |>
                        as.data.frame()|> pull(V2),
                      levels = c('1_vs_0','5_vs_1','5_vs_0'))
df_long$genotype <- factor(str_split(df_long$group_label, pattern = '_',n = 2, simplify = T) |>
                             as.data.frame()|> pull(V1),
                           levels = c('HomoDel','Hetero','HomoIns'))
df_long$group_label <- factor(df_long$group_label,
                              levels = c(
                                'HomoDel_1_vs_0',
                                'HomoDel_5_vs_1',
                                'HomoDel_5_vs_0',
                                'Hetero_1_vs_0',
                                'Hetero_5_vs_1',
                                'Hetero_5_vs_0',
                                'HomoIns_1_vs_0',
                                'HomoIns_5_vs_1',
                                'HomoIns_5_vs_0'
                              ))
df_long$group_label <- factor(df_long$group_label,
                              levels = c(
                                'HomoDel_1_vs_0',
                                'Hetero_1_vs_0',
                                'HomoIns_1_vs_0',
                                'HomoDel_5_vs_1',
                                'Hetero_5_vs_1',
                                'HomoIns_5_vs_1',
                                'HomoDel_5_vs_0',
                                'Hetero_5_vs_0',
                                'HomoIns_5_vs_0'
                              ))
df_long$genotype <- factor(df_long$genotype,
                              levels = c(
                                'HomoDel',
                                'Hetero',
                                'HomoIns'
                              ))
group_labels_vec <- c(
  'HomoDel_1_vs_0' = 'HomoDel',
  'Hetero_1_vs_0' = 'Hetero',
  'HomoIns_1_vs_0' = 'HomoIns',
  'HomoDel_5_vs_1' = 'HomoDel',
  'Hetero_5_vs_1' = 'Hetero',
  'HomoIns_5_vs_1' = 'HomoIns',
  'HomoDel_5_vs_0' = 'HomoDel',
  'Hetero_5_vs_0' = 'Hetero',
  'HomoIns_5_vs_0' = 'HomoIns'
)
group_labels_vec <- c(
  'HomoDel_1_vs_0' = 'HomoDel',
  'Hetero_1_vs_0' = 'Hetero',
  'HomoIns_1_vs_0' = 'HomoIns',
  'HomoDel_5_vs_1' = 'HomoDel',
  'Hetero_5_vs_1' = 'Hetero',
  'HomoIns_5_vs_1' = 'HomoIns'
)
max_x <- max(abs(df_long$Count)) * 1.1
# 绘图
df_long_new <- df_long |>
  filter(com != '5_vs_0')
  
ggplot(df_long_new, aes(x = Count, y = group_label, fill = com)) +
  
  # 背景色（先画，确保在底层）
  annotate("rect", xmin = -max_x, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#d6ecf3", alpha = 0.4) +  # 左淡蓝
  annotate("rect", xmin = 0, xmax = max_x, ymin = -Inf, ymax = Inf, fill = "#fde0dd", alpha = 0.4) +  # 右淡红
  
  # 柱状图
  geom_bar(stat = "identity", width = 0.7) +
  
  # 文字标签
  annotate("text", x = -max_x * 0.9, y = length(unique(df_long_new$group_label)) + 0.5,
           label = "Down", size = 1.5, fontface = "bold", color = "gray40") +
  annotate("text", x = max_x * 0.9, y = length(unique(df_long_new$group_label)) + 0.5,
           label = "Up", size = 1.5, fontface = "bold", color = "gray40") +
  
  # 配色
  scale_fill_manual(values = c("1_vs_0" = "#69ba8a", "5_vs_1" = "#e08b4a")) +
  
  # 主题美化
  theme_minimal(base_family = il10_font_family) +
  theme(
    axis.text.y       = element_text(size = 6),
    axis.text.x       = element_text(size = 6),
    legend.title      = element_blank(),
    legend.position   = 'none',
    panel.grid        = element_blank(),
    panel.background  = element_blank()
  ) +
  geom_text(
    aes(label = abs(Count)),
    position = position_stack(vjust = 0.5),
    hjust = ifelse(df_long_new$Direction == "Up", 0.4 , 0.4),  # Up 靠左，Down 靠右
    size = 1.5,
    color = "white"
  ) +
  scale_y_discrete(labels = group_labels_vec) +
  # 中线
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.3) +
  xlim(-max_x, max_x) +
  labs(x = NULL, y = NULL) -> p1
ggsave(filename = '02.DEGs/sp_degs_nums.pdf',
       plot = p1,
       device = 'pdf',
       width = 44.75,
       height = 33.86,
       units = 'mm')


 {
  int_file <- list.files(paste0(project_path("deg_lists"), .Platform$file.sep),"^.C[^v]*vs")
  library(stringr)
  library(xlsx)
  res_df <- data.frame('Genotype'= NA, 
                       'compare'= NA,
                       'Up' = NA, 
                       'down' = NA)
  for (i in seq_along(int_file)) {
    file <- int_file[i]
    Genotype <- if_else(str_sub(file,0,1) == 'A', 
                        'HomoIns', 
                        if_else(
                          str_sub(file,0,1) == 'B',
                          'HomoDel','Hetero'
                        ))
    treat <- str_sub(file,3,3)
    con <- str_sub(file,8,8)
    compare <- paste0(treat,'_vs_',con)
    # 读入数据
    df <- fread(paste0(paste0(project_path("deg_lists"), .Platform$file.sep),file))
    df <- df |> mutate(level = if_else(
      log2FoldChange > 1 & padj < 0.05,'Up',if_else(
        log2FoldChange < -1 & padj < 0.05,'Down','Nosign'
      )
    ))
    Up <- count(df, level) |> filter(level == 'Down') |> pull(n)
    Down <- count(df, level) |> filter(level == 'Up') |> pull(n)
    res_df[i,1] <- Genotype
    res_df[i,2] <- compare
    res_df[i,3] <- Up
    res_df[i,4] <- Down
    print(paste0('已完成：', i))
  }
  
  # 构造标签列（推荐用于y轴）
  df_deg <- res_df %>%
    mutate(group_label = paste(Genotype, compare, sep = "_"))
  
  # 转为长格式，Down 为负数
  df_long <- df_deg %>%
    dplyr::select(group_label, Up, down) %>%
    pivot_longer(cols = c("Up", "down"), names_to = "Direction", values_to = "Count") %>%
    mutate(Count = ifelse(Direction == "down", -Count, Count),
           Direction = ifelse(Direction == "down", "Down", "Up"))
  
  df_long$com <- factor(str_split(df_long$group_label, pattern = '_',n = 2, simplify = T) |>
                          as.data.frame()|> pull(V2),
                        levels = c('1_vs_0','5_vs_1','5_vs_0'))
  df_long$genotype <- factor(str_split(df_long$group_label, pattern = '_',n = 2, simplify = T) |>
                               as.data.frame()|> pull(V1),
                             levels = c('HomoDel','Hetero','HomoIns'))
  df_long$group_label <- factor(df_long$group_label,
                                levels = c(
                                  'HomoDel_1_vs_0',
                                  'HomoDel_5_vs_1',
                                  'HomoDel_5_vs_0',
                                  'Hetero_1_vs_0',
                                  'Hetero_5_vs_1',
                                  'Hetero_5_vs_0',
                                  'HomoIns_1_vs_0',
                                  'HomoIns_5_vs_1',
                                  'HomoIns_5_vs_0'
                                ))
  df_long$group_label <- factor(df_long$group_label,
                                levels = c(
                                  'HomoDel_1_vs_0',
                                  'Hetero_1_vs_0',
                                  'HomoIns_1_vs_0',
                                  'HomoDel_5_vs_1',
                                  'Hetero_5_vs_1',
                                  'HomoIns_5_vs_1',
                                  'HomoDel_5_vs_0',
                                  'Hetero_5_vs_0',
                                  'HomoIns_5_vs_0'
                                ))
  df_long$genotype <- factor(df_long$genotype,
                             levels = c(
                               'HomoDel',
                               'Hetero',
                               'HomoIns'
                             ))
  group_labels_vec <- c(
    'HomoDel_1_vs_0' = 'HomoDel',
    'Hetero_1_vs_0' = 'Hetero',
    'HomoIns_1_vs_0' = 'HomoIns',
    'HomoDel_5_vs_1' = 'HomoDel',
    'Hetero_5_vs_1' = 'Hetero',
    'HomoIns_5_vs_1' = 'HomoIns'
  )
  df_long_new <- df_long |>
    filter(com != '5_vs_0')
  max_x <- max(abs(df_long$Count)) * 1.1
  # 绘图
  ggplot(df_long_new, aes(x = Count, y = group_label, fill = com)) +
    
    # 背景色（先画，确保在底层）
    annotate("rect", xmin = -max_x, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#d6ecf3", alpha = 0.4) +  # 左淡蓝
    annotate("rect", xmin = 0, xmax = max_x, ymin = -Inf, ymax = Inf, fill = "#fde0dd", alpha = 0.4) +  # 右淡红
    
    # 柱状图
    geom_bar(stat = "identity", width = 0.7) +
    
    # 文字标签
    annotate("text", x = -max_x * 0.9, y = length(unique(df_long_new$group_label)) + 0.5,
             label = "Down", size = 1.5, fontface = "bold", color = "gray40") +
    annotate("text", x = max_x * 0.9, y = length(unique(df_long_new$group_label)) + 0.5,
             label = "Up", size = 1.5, fontface = "bold", color = "gray40") +
    
    # 配色
    scale_fill_manual(values = c("1_vs_0" = "#69ba8a", "5_vs_1" = "#e08b4a")) +
    
    # 主题美化
    theme_minimal(base_family = il10_font_family) +
    theme(
      axis.text.y       = element_text(size = 6),
      axis.text.x       = element_text(size = 6),
      legend.title      = element_blank(),
      legend.position   = 'none',
      panel.grid        = element_blank(),
      panel.background  = element_blank()
    ) +
    geom_text(
      data = df_long_new,
      aes(label = abs(Count)),   # 取绝对值显示，不显示负号
      hjust = ifelse(df_long_new$Direction == "Up", 1.5, -0.6),
      size = 1.5,
      color = "white"
    ) +
    scale_y_discrete(labels = group_labels_vec) +
    scale_x_continuous(limits = c(-3000, 3000), breaks = seq(-2000, 2000, by = 2000))+
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.3) +
    labs(x = NULL, y = NULL) -> p1
  ggsave(filename = '02.DEGs/int_degs_nums.pdf',
         plot = p1,
         device = 'pdf',
         width = 44.75,
         height = 33.86,
         units = 'mm')
  
  
}

