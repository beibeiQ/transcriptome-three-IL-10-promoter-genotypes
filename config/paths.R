# Project directories.
project_root <- function() {
  root <- Sys.getenv("IL10_PROJECT_ROOT", unset = .il10_root)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}
data_path <- function(...) file.path(Sys.getenv("IL10_DATA_DIR", unset = file.path(project_root(), "data")), ...)
annotation_path <- function(...) file.path(Sys.getenv("IL10_ANNOTATION_DIR", unset = data_path("annotation")), ...)
result_path <- function(...) file.path(Sys.getenv("IL10_RESULTS_DIR", unset = file.path(project_root(), "outputs")), ...)
project_path <- function(key) {
  switch(key,
    raw_counts = data_path("expression", "gene_count.xls"),
    fpkm = data_path("expression", "gene_fpkm.xls"),
    models = data_path("models"),
    figure2_analysis = result_path("figure2_analysis"),
    figure2_revision = result_path("figure2_revision"),
    contrast_inputs = result_path("contrast_inputs"),
    gene_contrasts = result_path("gene_contrasts"),
    deg_lists = data_path("figure1", "deg_lists"),
    source_tables = data_path("figure_source_tables"),
    stop("Unknown configured path key: ", key))
}
register_plot_fonts <- function() {
  folder <- Sys.getenv("IL10_FONT_DIR", unset = data_path("fonts"))
  files <- file.path(folder, c("ARIAL.TTF", "ARIALBD.TTF", "ARIALI.TTF", "ARIALBI.TTF"))
  if (all(file.exists(files)) && requireNamespace("sysfonts", quietly = TRUE)) {
    sysfonts::font_add("Arial", regular=files[1], bold=files[2], italic=files[3], bolditalic=files[4])
    return("Arial")
  }
  if (requireNamespace("sysfonts", quietly = TRUE) && "Arial" %in% sysfonts::font_families()) return("Arial")
  "sans"
}
il10_font_family <- register_plot_fonts()
for (folder in c(result_path(), project_path("figure2_analysis"), project_path("figure2_revision"),
                  project_path("contrast_inputs"), project_path("gene_contrasts"),
                  file.path(project_root(), "01.PCA", "il10 expression"),
                  file.path(project_root(), "01.PCA", "new_pca_plot"),
                  file.path(project_root(), "02.DEGs"), result_path("figure1")))
  dir.create(folder, recursive=TRUE, showWarnings=FALSE)
