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

library(showtext)
il10_font_family <- register_plot_fonts()
showtext_auto()




# =============================================================================
# 1. 脾脏IL10
# =============================================================================

source('load_data.R')
df_sp <- fpkm_plot |>
  filter(gene_id == '127522978' & Tissue == 'Sp')
df_sp$Genotype <- factor(df_sp$Genotype,
                         levels = c('HomoDel','Hetero','HomoIns'))
df_sp$days <- paste0('Days ',df_sp$Day)

df_summary <- df_sp %>%
  group_by(Genotype, Day) %>%
  summarise(mean_FPKM = mean(FPKM), .groups = "drop")

# 将 Day 转换为因子，按时间顺序排列
df_summary$Day <- factor(df_summary$Day, levels = c(0, 1, 5))

# 画图
library(RColorBrewer)
library(cowplot)
# 假设 df_summary 已存在
# Step 1: 计算每个 Day 内的最大值
bg_data <- df_summary %>%
  group_by(Day) %>%
  summarise(max_FPKM = max(mean_FPKM), .groups = "drop")

# Step 2: 合并最大值回原数据
df_plot <- df_summary %>%
  left_join(bg_data, by = "Day")

# Step 3: 绘图
ggplot(df_plot, aes(x = mean_FPKM, y = Genotype, fill = Day)) +
  # 添加灰色背景：xmin=0, xmax=max_FPKM
  geom_rect(data = df_plot,
            aes(xmin = 0, xmax = max_FPKM, ymin = as.numeric(Genotype) - 0.4, ymax = as.numeric(Genotype) + 0.4),
            fill = "gray", alpha = 0.3, inherit.aes = FALSE) +
  
  # 添加柱子
  geom_bar(stat = "identity") +
  
  # 分面图
  facet_wrap(~Day, ncol = 3) +
  scale_fill_brewer(palette = 'Set2') +
  labs(x = "", y = "") +
  theme_void()+
  theme(
    text = element_text(family = il10_font_family),
    strip.text = element_text(size = 6),
    strip.background = element_blank(),
    legend.position = 'none',
        axis.title.x = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)) -> p
ggsave(filename = '01.PCA/il10 expression/sp_il10_v2.pdf',
       plot = p,
       device = 'pdf',
       width = 685,
       height = 228,
       units = 'px')


# =============================================================================
# 2. 肠道IL10
# =============================================================================

df_int <- fpkm_plot |>
  filter(gene_id == '127522978' & Tissue == 'I')
df_int$Genotype <- factor(df_int$Genotype,
                         levels = c('HomoDel','Hetero','HomoIns'))
df_int$days <- paste0('Days ',df_int$Day)

df_summary <- df_int %>%
  group_by(Genotype, Day) %>%
  summarise(mean_FPKM = mean(FPKM), .groups = "drop")

# 将 Day 转换为因子，按时间顺序排列
df_summary$Day <- factor(df_summary$Day, levels = c(0, 1, 5))

# 画图
# Step 1: 计算每个 Day 内的最大值
bg_data <- df_summary %>%
  group_by(Day) %>%
  summarise(max_FPKM = max(mean_FPKM), .groups = "drop")

# Step 2: 合并最大值回原数据
df_plot <- df_summary %>%
  left_join(bg_data, by = "Day")

# Step 3: 绘图
ggplot(df_plot, aes(x = mean_FPKM, y = Genotype, fill = Day)) +
  # 添加灰色背景：xmin=0, xmax=max_FPKM
  geom_rect(data = df_plot,
            aes(xmin = 0, xmax = max_FPKM, ymin = as.numeric(Genotype) - 0.4, ymax = as.numeric(Genotype) + 0.4),
            fill = "gray", alpha = 0.3, inherit.aes = FALSE) +
  
  # 添加柱子
  geom_bar(stat = "identity") +
  
  # 分面图
  facet_wrap(~Day, ncol = 3) +
  scale_fill_brewer(palette = 'Set2') +
  labs(x = "", y = "") +
  theme_void()+
  theme(
    text = element_text(family = il10_font_family),
    strip.text = element_text(size = 6),
    strip.background = element_blank(),
    legend.position = 'none',
    axis.title.x = element_text(size = 6),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)) -> p
ggsave(filename = '01.PCA/il10 expression/int_il10_v2.pdf',
       plot = p,
       device = 'pdf',
       width = 685,
       height = 228,
       units = 'px')
