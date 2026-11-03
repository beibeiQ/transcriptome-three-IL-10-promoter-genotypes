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

suppressPackageStartupMessages({library(DESeq2);library(data.table)})
args<-commandArgs(TRUE);base<-if(length(args))args[1] else project_root()
out<-project_path("gene_contrasts");inp<-file.path(out,'inputs');dir.create(inp,recursive=TRUE,showWarnings=FALSE)
src<-file.path(project_path("contrast_inputs"),'inputs')
for(n in c('raw_counts_input.tsv','Spleen_dds.rds','Intestine_dds.rds','sample_metadata.csv'))if(!file.exists(file.path(inp,n)))stopifnot(file.copy(file.path(src,n),file.path(inp,n)))
raw<-fread(file.path(inp,'raw_counts_input.tsv'));raw[,gene_id:=as.character(gene_id)]
ann<-raw[,.(gene_id,gene_name,gene_description)];fwrite(ann,file.path(out,'gene_annotations.csv.gz'))
all<-list();qa<-list();filters<-list();definitions<-list();nc<-list()
for(tis in c('Spleen','Intestine')){
 dds<-readRDS(file.path(inp,paste0(tis,'_dds.rds')));md<-as.data.frame(colData(dds));md$sample<-colnames(dds)
 m<-as.matrix(raw[,md$sample,with=FALSE]);rownames(m)<-raw$gene_id;kk<-rowSums(m>=10)>=3
 stopifnot(identical(rownames(dds),raw$gene_id[kk]),all(counts(dds)==m[kk,]),all(table(md$genotype,md$time)==3),all(mcols(dds)$betaConv))
 X<-model.matrix(~genotype*time,md);stopifnot(qr(X)$rank==9)
 masks<-sapply(split(seq_len(nrow(md)),interaction(md$genotype,md$time)),function(j)rowSums(m[,j,drop=FALSE]>=10)>=2)
 filters[[tis]]<-data.table(tissue=tis,gene_id=raw$gene_id,kept=kk,strict=apply(masks,1,all))
 rd<-function(g,t)as.numeric(model.matrix(~genotype*time,data.frame(genotype=factor(g,levels=levels(md$genotype)),time=factor(t,levels=levels(md$time)))))
 for(g in c('HomoDel','HomoIns'))for(t in c('d0','d1','d5')){
  v<-rd('Hetero',t)-rd(g,t);stage<-if(t=='d0')'Baseline'else paste0('Response_',t)
  if(t!='d0')v<-v-rd('Hetero','d0')+rd(g,'d0')
  test<-paste(stage,'Hetero_vs',g,sep='_')
  z<-as.data.table(as.data.frame(results(dds,contrast=v,independentFiltering=FALSE)),keep.rownames='gene_id')
  setnames(z,c('log2FoldChange','pvalue','padj'),c('effect','pvalue','q_within'))
  z[,`:=`(tissue=tis,stage=stage,comparator=g,test=test,linear_effect=as.numeric(coef(dds)%*%v),zero_count_convention=as.logical(DESeq2:::contrastAllZeroNumeric(dds,v))&!as.logical(mcols(dds)$allZero),ci95_low=effect-1.96*lfcSE,ci95_high=effect+1.96*lfcSE)]
  z[,q_within:=p.adjust(pvalue,'BH')];z<-merge(z,filters[[tis]][kept==TRUE,.(gene_id,strict)],by='gene_id')
  stopifnot(all(abs(z[zero_count_convention==FALSE,effect-linear_effect])<1e-5,na.rm=TRUE))
  all[[paste(tis,test)]]<-z;definitions[[paste(tis,test)]]<-data.table(tissue=tis,test=test,coefficient=colnames(X),weight=v)
 }
 n<-as.data.table(counts(dds,normalized=TRUE),keep.rownames='gene_id');n<-melt(n,id.vars='gene_id',variable.name='sample',value.name='normalized_count')
 n<-merge(n,as.data.table(md)[,.(sample,genotype,time,rep,fish_id)],by='sample');n[,tissue:=tis];nc[[tis]]<-n
 qa[[tis]]<-data.table(tissue=tis,genes=nrow(dds),samples=ncol(dds),raw_identical=TRUE,rank=qr(X)$rank)
}
z<-rbindlist(all);z[,family:=ifelse(stage=='Baseline','Baseline','Response')];z[,q_pair_family:=p.adjust(pvalue,'BH'),by=family]
z<-merge(z,ann,by='gene_id',all.x=TRUE)
fwrite(z,file.path(out,'DESeq2_pairwise.csv.gz'));fwrite(rbindlist(nc),file.path(out,'normalized_counts_long.csv.gz'))
fwrite(rbindlist(filters),file.path(out,'gene_filter_audit.csv.gz'));fwrite(rbindlist(qa),file.path(out,'DESeq2_model_QA.csv'));fwrite(rbindlist(definitions),file.path(out,'DESeq2_contrast_definitions.csv'))
fwrite(data.table(file=list.files(inp),md5=unname(tools::md5sum(list.files(inp,full.names=TRUE)))),file.path(out,'input_checksums.csv'))
capture.output(sessionInfo(),file=file.path(out,'DESeq2_sessionInfo.txt'))
