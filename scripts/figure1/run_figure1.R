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


source(file.path(project_root(), "scripts/figure1/helpers.R"))
args <- commandArgs(TRUE)
panel <- if(length(args)) toupper(args[1]) else "CHECK"
if(panel=="CHECK") {
  needed <- c(project_path("fpkm"),annotation_path("grass_carp_hallmark_list.rds"))
  cat("Input presence:\n");print(data.frame(file=basename(needed),present=file.exists(needed)),row.names=FALSE)
  quit(status=if(all(file.exists(needed))) 0 else 2)
}
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
if (!file.exists(project_path("fpkm"))) stop("Missing input: data/expression/gene_fpkm.xls")
raw <- fread(project_path("fpkm")); ids <- grep("^[ABC][CP][015]_[123]$",names(raw),value=TRUE)
stopifnot(length(ids)==54, "gene_id" %in% names(raw),!anyDuplicated(raw$gene_id))
meta <- decode_samples(ids);stopifnot(all(table(meta$genotype,meta$tissue,meta$day)==3))
mat <- as.matrix(raw[,..ids]);rownames(mat)<-as.character(raw$gene_id)
stopifnot(all(is.finite(mat)),all(mat>=0))
colors <- c(HomoDel="#a8dad8",Hetero="#825ca6",HomoIns="#e7872b")
out <- result_path("figure1");dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(panel=="C") {
  stopifnot("127522978" %in% rownames(mat))
  df <- meta;df$FPKM <- as.numeric(mat["127522978",]);df$Day <- factor(df$day,levels=c(0,1,5))
  df$genotype <- factor(df$genotype,levels=c("HomoDel","Hetero","HomoIns"))
  for(tis in c("Spleen","Intestine")) {
    d <- df[df$tissue==tis,]
    med <- aggregate(FPKM~genotype+Day,d,median)
    p <- ggplot(d,aes(Day,FPKM,fill=genotype))+geom_boxplot(linewidth=.2)+
      geom_line(data=med,aes(group=genotype),color="black",linewidth=.3)+
      geom_point(data=med,color="black",size=.3)+facet_wrap(~genotype,nrow=1)+
      scale_fill_manual(values=colors)+theme_bw(base_family=il10_font_family)+
      theme(legend.position="none")+labs(x="Days post-challenge",y="il10 abundance (FPKM)")
    ggsave(file.path(out,paste0("C_il10_",tis,".pdf")),p,width=135,height=65,units="mm")
  }
} else if(panel %in% c("D","E")) {
  if(!requireNamespace("PCAtools",quietly=TRUE)) stop("Install PCAtools for tissue/time PCA panels")
  tis <- if(panel=="D") "Spleen" else "Intestine"
  for(day in c(0,1,5)) {
    md <- meta[meta$tissue==tis & meta$day==day,];rownames(md)<-md$sample
    x <- mat[,md$sample,drop=FALSE];x<-x[rowSums(x)>0,,drop=FALSE]
    fit <- PCAtools::pca(x,metadata=md,removeVar=.1,center=TRUE,scale=FALSE)
    p <- PCAtools::biplot(fit,colby="genotype",encircle=TRUE,legendPosition="none",title=paste(tis,day,"dpi"))
    ggsave(file.path(out,paste0(panel,"_PCA_",day,"dpi.pdf")),p,width=80,height=65,units="mm")
  }
} else if(panel %in% c("H","I")) {
  scorefile <- file.path(out,"gsva.es.hallmark.csv")
  if(!file.exists(scorefile)) stop("Missing input: outputs/figure1/gsva.es.hallmark.csv")
  score <- fread(scorefile);pathway<-score[[1]];score<-as.data.frame(score[,-1]);rownames(score)<-pathway
  md <- decode_samples(names(score));tis<-if(panel=="H") "Spleen" else "Intestine"
  focus<-c("Tnfa_signaling_via_nfkb","Reactive_oxygen_species_pathway","Interferon_gamma_response","Interferon_alpha_response","Inflammatory_response","Il6_jak_stat3_signaling","Il2_stat5_signaling","Complement","Apoptosis")
  stopifnot(all(focus %in% rownames(score)));sets<-list()
  for(day in c(0,1,5)) for(gene_set in focus) {
    means<-vapply(c("HomoDel","Hetero","HomoIns"),function(g) {
      cols<-md$sample[md$tissue==tis & md$day==day & md$genotype==g]
      stopifnot(length(cols)==3);mean(as.numeric(score[gene_set,cols]))
    },numeric(1))
    sets[[length(sets)+1]]<-data.frame(pathway=gene_set,day=day,genotype=names(means),mean_score=means,z=standardize_genotype_means(means))
  }
  d<-do.call(rbind,sets);d$pathway<-factor(d$pathway,levels=rev(focus));d$genotype<-factor(d$genotype,levels=c("HomoDel","Hetero","HomoIns"))
  fwrite(d,file.path(out,paste0(panel,"_heatmap_values.csv")))
  p<-ggplot(d,aes(genotype,pathway,fill=z))+geom_tile(color="white")+facet_wrap(~day,nrow=1)+
     scale_fill_gradient2(low="#5BB770",mid="white",high="#E88A3D",limits=c(-1.2,1.2),oob=scales::squish)+
     theme_minimal(base_family=il10_font_family)+theme(axis.text.x=element_text(angle=45,hjust=1))+
     labs(x=NULL,y=NULL,fill="Within-panel Z")
  ggsave(file.path(out,paste0(panel,"_immune_heatmap.pdf")),p,width=175,height=95,units="mm")
} else stop("Choose CHECK, C, D, E, H or I; DEG input directions are specified separately.")
