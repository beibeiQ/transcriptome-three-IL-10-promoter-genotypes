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

library(readr)
gsva_es_kegg <- read_csv(result_path("figure1", "gsva_es_kegg.csv"))
View(gsva_es_kegg)
library(dplyr)

gsva_es_kegg <- gsva_es_kegg |> column_to_rownames('...1')
library(pheatmap)
pheatmap(gsva_es_kegg,show_rownames = F, show_colnames = F)

gsva_es_hallmark <- read_csv(result_path("figure1", "gsva.es.hallmark.csv"))
gsva_es_hallmark <- gsva_es_hallmark |> column_to_rownames('...1')
gsva_es_hallmark_sp  <- gsva_es_hallmark |> dplyr::select(contains('_Sp_'))


gsva_es_hallmark_sp_0_mean <- data.frame()
gsva_es_hallmark_sp_1_mean <- data.frame()
gsva_es_hallmark_sp_5_mean <- data.frame()

focus_hallmark <- c('Inflammatory_response',
                    'Interferon_alpha_response',
                    'Interferon_gamma_response',
                    'Tnfa_signaling_via_nfkb',
                    'Il2_stat5_signaling',
                    'Il6_jak_stat3_signaling',
                    'Reactive_oxygen_species_pathway',
                    'Complement',
                    'Apoptosis')

gsva_es_hallmark_sp_0  <- gsva_es_hallmark_sp |> 
  dplyr::select(contains('_0d_'))
gsva_es_hallmark_sp_0 <- gsva_es_hallmark_sp_0 |>
  mutate(
    HomoDel_mean = gsva_es_hallmark_sp_0 |>dplyr::select(starts_with('HomoDel')) |> apply(1,mean),
    Hetero_mean = gsva_es_hallmark_sp_0 |>dplyr::select(starts_with('Hetero')) |> apply(1,mean),
    HomoIns_mean= gsva_es_hallmark_sp_0 |>dplyr::select(starts_with('HomoIns')) |> apply(1,mean)
  )
gsva_es_hallmark_sp_0_mean <-  gsva_es_hallmark_sp_0[,10:12]


gsva_es_hallmark_sp_1  <- gsva_es_hallmark_sp |> 
  dplyr::select(contains('_1d_'))
gsva_es_hallmark_sp_1 <- gsva_es_hallmark_sp_1 |>
  mutate(
    HomoDel_mean = gsva_es_hallmark_sp_1 |>dplyr::select(starts_with('HomoDel')) |> apply(1,mean),
    Hetero_mean = gsva_es_hallmark_sp_1 |>dplyr::select(starts_with('Hetero')) |> apply(1,mean),
    HomoIns_mean= gsva_es_hallmark_sp_1 |>dplyr::select(starts_with('HomoIns')) |> apply(1,mean)
  )
gsva_es_hallmark_sp_1_mean <-  gsva_es_hallmark_sp_1[,10:12]


gsva_es_hallmark_sp_5  <- gsva_es_hallmark_sp |> 
  dplyr::select(contains('_5d_'))
gsva_es_hallmark_sp_5 <- gsva_es_hallmark_sp_5 |>
  mutate(
    HomoDel_mean = gsva_es_hallmark_sp_5 |>dplyr::select(starts_with('HomoDel')) |> apply(1,mean),
    Hetero_mean = gsva_es_hallmark_sp_5 |>dplyr::select(starts_with('Hetero')) |> apply(1,mean),
    HomoIns_mean= gsva_es_hallmark_sp_5 |>dplyr::select(starts_with('HomoIns')) |> apply(1,mean)
  )
gsva_es_hallmark_sp_5_mean <- gsva_es_hallmark_sp_5[,10:12]

names(gsva_es_hallmark_sp_0_mean) <- c('HomoDel', 'Hetero','HomoIns')
names(gsva_es_hallmark_sp_1_mean) <- c('HomoDel', 'Hetero','HomoIns')
names(gsva_es_hallmark_sp_5_mean) <- c('HomoDel', 'Hetero','HomoIns')
new_name <- c('Apoptosis',
              'Complement',
              'Il2_stat5',
              'Il6_jak',
              'Inflammatory',
              'Interferon_alpha',
              'Interferon_gamma',
              'Reactive_oxygen',
              'Tnfa_nfkb')
df_sp_0 <- gsva_es_hallmark_sp_0_mean[rownames(gsva_es_hallmark_sp_0_mean) %in% focus_hallmark,]
rownames(df_sp_0) <- new_name
df_sp_1 <- gsva_es_hallmark_sp_1_mean[rownames(gsva_es_hallmark_sp_1_mean) %in% focus_hallmark,]
rownames(df_sp_1) <- new_name
df_sp_5 <- gsva_es_hallmark_sp_5_mean[rownames(gsva_es_hallmark_sp_5_mean) %in% focus_hallmark,]
rownames(df_sp_5) <- new_name


library(colorspace)
hcl_palettes(plot = T)
breaks <- seq(-1.2, 1.2, length.out = 100)
colors <- colorRampPalette(c('white','#f4b1bd','#ea5a72'))(length(breaks))

colors <- colorRampPalette(c('#B5E0C0','#FFFFFF','#FDBF92'))(9)
colors <- colorRampPalette(c('#5BB770','#B5E0C0','#FFFFFF','#FDBF92','#E88A3D'))(9)
colors <- colorRampPalette(c('#F9F9F9','#3E6EB0'))(length(breaks))

df_scaled <- t(apply(df_sp_0, 1, scale))  # 按行标准化
colnames(df_scaled) <- colnames(df_sp_0)
rownames(df_scaled) <- rownames(df_sp_0)
df_long <- as.data.frame(df_scaled) %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(cols = -Pathway, names_to = "Group", values_to = "Zscore")
df_long$Group <- factor(df_long$Group, levels = c('HomoDel','Hetero','HomoIns'))
ggplot(df_long, aes(x = Group, y = Pathway, fill = Zscore)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = colors,
                       limits = c(-1.2, 1.2), name = "Z-score") +
  theme_minimal(base_size = 6) +
  theme(text = element_text(family = il10_font_family),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none') -> p1
ggsave(filename = '02.DEGs/heatmap_sp_0.pdf',
       plot = p1,
       device = 'pdf',
       width = 44.75,
       height = 33.86,
       units = 'mm')

df_scaled <- t(apply(df_sp_1, 1, scale))  # 按行标准化
colnames(df_scaled) <- colnames(df_sp_1)
rownames(df_scaled) <- rownames(df_sp_1)
df_long <- as.data.frame(df_scaled) %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(cols = -Pathway, names_to = "Group", values_to = "Zscore")
df_long$Group <- factor(df_long$Group, levels = c('HomoDel','Hetero','HomoIns'))

ggplot(df_long, aes(x = Group, y = Pathway, fill = Zscore)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = colors,
                       limits = c(-1.2, 1.2), name = "Z-score") +
  theme_minimal(base_size = 6) +
  theme(text = element_text(family = il10_font_family),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 6),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none') -> p1
ggsave(filename = '02.DEGs/heatmap_sp_1.pdf',
       plot = p1,
       device = 'pdf',
       width = 26.75,
       height = 33.86,
       units = 'mm')

df_scaled <- t(apply(df_sp_5, 1, scale))  # 按行标准化
colnames(df_scaled) <- colnames(df_sp_5)
rownames(df_scaled) <- rownames(df_sp_5)
df_long <- as.data.frame(df_scaled) %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(cols = -Pathway, names_to = "Group", values_to = "Zscore")
df_long$Group <- factor(df_long$Group, levels = c('HomoDel','Hetero','HomoIns'))

ggplot(df_long, aes(x = Group, y = Pathway, fill = Zscore)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = colors,
                       limits = c(-1.2, 1.2), name = "Z-score") +
  theme_minimal(base_size = 6) +
  theme(text = element_text(family = il10_font_family),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 6),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none') -> p1
ggsave(filename = '02.DEGs/heatmap_sp_5.pdf',
       plot = p1,
       device = 'pdf',
       width = 26.75,
       height = 33.86,
       units = 'mm')

{
  gsva_es_hallmark_int  <- gsva_es_hallmark |> dplyr::select(contains('_I_'))
  
  
  gsva_es_hallmark_int_0_mean <- data.frame()
  gsva_es_hallmark_int_1_mean <- data.frame()
  gsva_es_hallmark_int_5_mean <- data.frame()
  
  focus_hallmark <- c('Inflammatory_response',
                      'Interferon_alpha_response',
                      'Interferon_gamma_response',
                      'Tnfa_signaling_via_nfkb',
                      'Il2_stat5_signaling',
                      'Il6_jak_stat3_signaling',
                      'Reactive_oxygen_species_pathway',
                      'Complement',
                      'Apoptosis')
  
  gsva_es_hallmark_int_0  <- gsva_es_hallmark_int |> 
    dplyr::select(contains('_0d_'))
  gsva_es_hallmark_int_0 <- gsva_es_hallmark_int_0 |>
    mutate(
      HomoDel_mean = gsva_es_hallmark_int_0 |>dplyr::select(starts_with('HomoDel')) |> apply(1,mean),
      Hetero_mean = gsva_es_hallmark_int_0 |>dplyr::select(starts_with('Hetero')) |> apply(1,mean),
      HomoIns_mean= gsva_es_hallmark_int_0 |>dplyr::select(starts_with('HomoIns')) |> apply(1,mean)
    )
  gsva_es_hallmark_int_0_mean <-  gsva_es_hallmark_int_0[,10:12]
  
  
  gsva_es_hallmark_int_1  <- gsva_es_hallmark_int |> 
    dplyr::select(contains('_1d_'))
  gsva_es_hallmark_int_1 <- gsva_es_hallmark_int_1 |>
    mutate(
      HomoDel_mean = gsva_es_hallmark_int_1 |>dplyr::select(starts_with('HomoDel')) |> apply(1,mean),
      Hetero_mean = gsva_es_hallmark_int_1 |>dplyr::select(starts_with('Hetero')) |> apply(1,mean),
      HomoIns_mean= gsva_es_hallmark_int_1 |>dplyr::select(starts_with('HomoIns')) |> apply(1,mean)
    )
  gsva_es_hallmark_int_1_mean <-  gsva_es_hallmark_int_1[,10:12]
  
  
  gsva_es_hallmark_int_5  <- gsva_es_hallmark_int |> 
    dplyr::select(contains('_5d_'))
  gsva_es_hallmark_int_5 <- gsva_es_hallmark_int_5 |>
    mutate(
      HomoDel_mean = gsva_es_hallmark_int_5 |>dplyr::select(starts_with('HomoDel')) |> apply(1,mean),
      Hetero_mean = gsva_es_hallmark_int_5 |>dplyr::select(starts_with('Hetero')) |> apply(1,mean),
      HomoIns_mean= gsva_es_hallmark_int_5 |>dplyr::select(starts_with('HomoIns')) |> apply(1,mean)
    )
  gsva_es_hallmark_int_5_mean <- gsva_es_hallmark_int_5[,10:12]
  
  names(gsva_es_hallmark_int_0_mean) <- c('HomoDel', 'Hetero','HomoIns')
  names(gsva_es_hallmark_int_1_mean) <- c('HomoDel', 'Hetero','HomoIns')
  names(gsva_es_hallmark_int_5_mean) <- c('HomoDel', 'Hetero','HomoIns')
  new_name <- c('Apoptosis',
                'Complement',
                'Il2_stat5',
                'Il6_jak',
                'Inflammatory',
                'Interferon_alpha',
                'Interferon_gamma',
                'Reactive_oxygen',
                'Tnfa_nfkb')
  df_int_0 <- gsva_es_hallmark_int_0_mean[rownames(gsva_es_hallmark_int_0_mean) %in% focus_hallmark,]
  rownames(df_int_0) <- new_name
  df_int_1 <- gsva_es_hallmark_int_1_mean[rownames(gsva_es_hallmark_int_1_mean) %in% focus_hallmark,]
  rownames(df_int_1) <- new_name
  df_int_5 <- gsva_es_hallmark_int_5_mean[rownames(gsva_es_hallmark_int_5_mean) %in% focus_hallmark,]
  rownames(df_int_5) <- new_name
  
  

  colors <- colorRampPalette(c('#5BB770','#B5E0C0','#FFFFFF','#FDBF92','#E88A3D'))(9)
  
  
  df_scaled <- t(apply(df_int_0, 1, scale))  # 按行标准化
  colnames(df_scaled) <- colnames(df_int_0)
  rownames(df_scaled) <- rownames(df_int_0)
  df_long <- as.data.frame(df_scaled) %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(cols = -Pathway, names_to = "Group", values_to = "Zscore")
  df_long$Group <- factor(df_long$Group, levels = c('HomoDel','Hetero','HomoIns'))
  ggplot(df_long, aes(x = Group, y = Pathway, fill = Zscore)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = colors,
                         limits = c(-1.2, 1.2), name = "Z-score") +
    theme_minimal(base_size = 6) +
    theme(text = element_text(family = il10_font_family),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 6),
          axis.text.y = element_text(size = 6),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none') -> p1
  ggsave(filename = '02.DEGs/heatmap_int_0.pdf',
         plot = p1,
         device = 'pdf',
         width = 44.75,
         height = 33.86,
         units = 'mm')
  
  df_scaled <- t(apply(df_int_1, 1, scale))  # 按行标准化
  colnames(df_scaled) <- colnames(df_int_1)
  rownames(df_scaled) <- rownames(df_int_1)
  df_long <- as.data.frame(df_scaled) %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(cols = -Pathway, names_to = "Group", values_to = "Zscore")
  df_long$Group <- factor(df_long$Group, levels = c('HomoDel','Hetero','HomoIns'))
  
  ggplot(df_long, aes(x = Group, y = Pathway, fill = Zscore)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = colors,
                         limits = c(-1.2, 1.2), name = "Z-score") +
    theme_minimal(base_size = 6) +
    theme(text = element_text(family = il10_font_family),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 6),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none') -> p1
  ggsave(filename = '02.DEGs/heatmap_int_1.pdf',
         plot = p1,
         device = 'pdf',
         width = 26.75,
         height = 33.86,
         units = 'mm')
  
  df_scaled <- t(apply(df_int_5, 1, scale))  # 按行标准化
  colnames(df_scaled) <- colnames(df_int_5)
  rownames(df_scaled) <- rownames(df_int_5)
  df_long <- as.data.frame(df_scaled) %>%
    rownames_to_column("Pathway") %>%
    pivot_longer(cols = -Pathway, names_to = "Group", values_to = "Zscore")
  df_long$Group <- factor(df_long$Group, levels = c('HomoDel','Hetero','HomoIns'))
  
  ggplot(df_long, aes(x = Group, y = Pathway, fill = Zscore)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = colors,
                         limits = c(-1.2, 1.2), name = "Z-score") +
    theme_minimal(base_size = 6) +
    theme(text = element_text(family = il10_font_family),
          axis.text.x = element_text(angle = 45, hjust = 1,size = 6),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none') -> p1
  ggsave(filename = '02.DEGs/heatmap_int_5.pdf',
         plot = p1,
         device = 'pdf',
         width = 26.75,
         height = 33.86,
         units = 'mm')
}


