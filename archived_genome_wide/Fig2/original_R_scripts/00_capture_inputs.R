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

suppressPackageStartupMessages(library(data.table))
args<-commandArgs(TRUE);base<-if(length(args))args[1] else project_root()
src<-project_path("figure2_analysis");out<-file.path(project_path("contrast_inputs"),'inputs');dir.create(out,recursive=TRUE,showWarnings=FALSE)
needed<-c('gene_annotations.csv.gz','gene_filter_audit.csv.gz','sample_metadata.csv','KEGG_annotation_input.tsv')
stopifnot(all(file.copy(file.path(src,needed),file.path(out,needed),overwrite=TRUE)))
raw<-file.path(src,'raw_counts_input.tsv');if(!file.exists(raw))raw<-project_path("raw_counts")
stopifnot(file.copy(raw,file.path(out,'raw_counts_input.tsv'),overwrite=TRUE))
fits<-file.path(project_path("models"),c('Spleen_dds.rds','Intestine_dds.rds'))
stopifnot(all(file.copy(fits,file.path(out,basename(fits)),overwrite=TRUE)))
co<-fread(cmd=paste('gzip -dc',shQuote(file.path(src,'all_contrasts.csv.gz'))),colClasses=c(gene_id='character'))
co<-co[family=='Response_pairwise'];co[,c('gene_name','gene_description'):=NULL]
fwrite(co,file.path(out,'all_contrasts.csv.gz'))
paths<-list.files(out,full.names=TRUE)
fwrite(data.table(file=basename(paths),md5=unname(tools::md5sum(paths))),file.path(dirname(out),'input_checksums.csv'))
cat('INPUT CAPTURE COMPLETE\n')
