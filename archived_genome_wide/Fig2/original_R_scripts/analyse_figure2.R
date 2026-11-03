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

suppressPackageStartupMessages({library(DESeq2); library(data.table)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)>=1) args[1] else project_path("models")
src <- if(length(args)>=2) args[2] else project_path("raw_counts")
out <- if(length(args)>=3) args[3] else project_path("figure2_analysis")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
raw <- fread(src); raw[,gene_id:=as.character(gene_id)]
annot <- unique(raw[,.(gene_id,gene_name,gene_description)])
stopifnot(nrow(annot)==nrow(raw), !anyDuplicated(raw$gene_id))
fwrite(annot,file.path(out,'gene_annotations.csv.gz'))
input_paths<-c(src,file.path(root,c('Spleen_dds.rds','Intestine_dds.rds')))
fwrite(data.table(file=basename(input_paths[file.exists(input_paths)]), md5=unname(tools::md5sum(input_paths[file.exists(input_paths)]))), file.path(out,'input_md5.tsv'), sep='\t')
omnibus <- list(); cons <- list(); ncounts <- list(); quality <- list(); metas <- list(); filters <- list(); sizes <- list(); coeffs <- list()
for(tis in c('Spleen','Intestine')) {
 cat('START',tis,format(Sys.time()),'\n'); flush.console()
 fit_path<-file.path(root,paste0(tis,'_dds.rds'))
 if(file.exists(fit_path)) {
   dds <- readRDS(fit_path)
 } else {
   ids<-grep('^[ABC][CP][015]_[123]$',names(raw),value=TRUE)
   ids<-ids[substr(ids,2,2)==if(tis=='Spleen')'P' else 'C']
   meta<-data.frame(sample=ids,genotype=factor(unname(c(A='HomoIns',B='HomoDel',C='Hetero')[substr(ids,1,1)]),levels=c('HomoDel','Hetero','HomoIns')),tissue=tis,time=factor(paste0('d',substr(ids,3,3)),levels=c('d0','d1','d5')),rep=substr(ids,5,5))
   meta$fish_id<-paste0(substr(ids,1,1),meta$time,'_',meta$rep);meta$tank<-NA_character_;rownames(meta)<-ids
   m<-as.matrix(raw[,ids,with=FALSE]);rownames(m)<-raw$gene_id
   kk<-rowSums(m>=10)>=3
   dds<-DESeq(DESeqDataSetFromMatrix(m[kk,,drop=FALSE],meta,~genotype*time),quiet=TRUE)
   saveRDS(dds,file.path(out,paste0(tis,'_dds.rds')))
 }
 md <- as.data.frame(colData(dds)); md$sample <- colnames(dds)
 stopifnot(nrow(md)==27,all(table(md$genotype,md$time)==3), !anyDuplicated(md$sample))
 stopifnot(identical(levels(md$genotype),c('HomoDel','Hetero','HomoIns')),identical(levels(md$time),c('d0','d1','d5')))
 expected_g <- unname(c(A='HomoIns',B='HomoDel',C='Hetero')[substr(md$sample,1,1)])
 stopifnot(all(as.character(md$genotype)==expected_g),all(as.character(md$time)==paste0('d',substr(md$sample,3,3))))
 mat <- as.matrix(raw[,md$sample,with=FALSE]); rownames(mat)<-raw$gene_id
 stopifnot(all(is.finite(mat)),all(mat>=0),all(mat==round(mat)))
 keep <- rowSums(mat>=10)>=3
 stopifnot(identical(rownames(dds),raw$gene_id[keep]),all(counts(dds)==mat[keep,,drop=FALSE]))
 X <- model.matrix(~genotype*time,md)
 stopifnot(qr(X)$rank==9,ncol(coef(dds))==9)
 # Dropping genotype main effects ONLY constrains equality at reference day 0.
 base_reduced <- X[,!colnames(X)%in%c('genotypeHetero','genotypeHomoIns'),drop=FALSE]
 response_reduced <- model.matrix(~genotype+time,md)
 stopifnot(ncol(base_reduced)==7,qr(base_reduced)$rank==7,qr(response_reduced)$rank==5)
 masks <- sapply(split(seq_len(nrow(md)),interaction(md$genotype,md$time,drop=TRUE)),function(j)rowSums(mat[,j,drop=FALSE]>=10)>=2)
 fil <- data.table(gene_id=raw$gene_id,tissue=tis,kept=keep,expressed_all_groups=apply(masks,1,all))
 filters[[tis]] <- fil
 for(kind in c('Baseline','Response')) {
   reduced <- if(kind=='Baseline')base_reduced else response_reduced
   lrt <- nbinomLRT(dds,full=X,reduced=reduced,quiet=TRUE)
   tab <- as.data.table(as.data.frame(results(lrt,independentFiltering=FALSE)),keep.rownames='gene_id')
   tab[,c('log2FoldChange','lfcSE'):=NULL]
   setnames(tab,c('stat','padj'),c('LRT_stat','q_within'))
   tab[,`:=`(tissue=tis,test=kind,df=ncol(X)-ncol(reduced),pvalue=pvalue)]
   tab[,q_within:=p.adjust(pvalue,'BH')]
   omnibus[[paste(tis,kind)]] <- merge(tab,fil[kept==TRUE,.(gene_id,expressed_all_groups)],by='gene_id')
   cat(kind,'finite',sum(is.finite(tab$pvalue)),'q05',sum(tab$q_within<.05,na.rm=TRUE),'\n');flush.console()
 }
 rowX <- function(g,t) as.numeric(model.matrix(~genotype*time,data.frame(genotype=factor(g,levels=levels(md$genotype)),time=factor(t,levels=levels(md$time)))))
 defs <- list()
 pairs <- list(c('Hetero','HomoDel'),c('HomoIns','HomoDel'),c('HomoIns','Hetero'))
 for(pa in pairs) {
   nm <- paste(pa,collapse='_vs_')
   defs[[paste0('Baseline_',nm)]] <- list(v=rowX(pa[1],'d0')-rowX(pa[2],'d0'),family='Baseline_pairwise',g1=pa[1],g0=pa[2],time='d0')
   for(t in c('d1','d5')) defs[[paste0('Response_',nm,'_',t)]] <- list(v=rowX(pa[1],t)-rowX(pa[1],'d0')-rowX(pa[2],t)+rowX(pa[2],'d0'),family='Response_pairwise',g1=pa[1],g0=pa[2],time=t)
 }
 for(g in levels(md$genotype))for(t in c('d1','d5'))defs[[paste0('Within_',g,'_',t)]]<-list(v=rowX(g,t)-rowX(g,'d0'),family='Within_genotype',g1=g,g0=g,time=t)
 for(nm in names(defs)) {
   a<-defs[[nm]]; tab<-as.data.table(as.data.frame(results(dds,contrast=a$v,independentFiltering=FALSE)),keep.rownames='gene_id')
   setnames(tab,'padj','q_within')
   tab[,`:=`(tissue=tis,test=nm,family=a$family,g1=a$g1,g0=a$g0,time=a$time,ci95_low=log2FoldChange-1.96*lfcSE,ci95_high=log2FoldChange+1.96*lfcSE)]
   tab[,linear_predictor_log2FC:=as.numeric(coef(dds)%*%a$v)]
   tab[,zero_count_convention:=as.logical(DESeq2:::contrastAllZeroNumeric(dds,a$v)) & !as.logical(mcols(dds)$allZero)]
   stopifnot(all(abs(tab$log2FoldChange[!tab$zero_count_convention]-tab$linear_predictor_log2FC[!tab$zero_count_convention])<1e-5,na.rm=TRUE))
   tab[,q_within:=p.adjust(pvalue,'BH')];cons[[paste(tis,nm)]]<-tab
 }
 nc<-as.data.table(counts(dds,normalized=TRUE),keep.rownames='gene_id')
 nc<-melt(nc,id.vars='gene_id',variable.name='sample',value.name='normalized_count')
 nc<-merge(nc,as.data.table(md)[,.(sample,genotype,time,rep,fish_id)],by='sample');nc[,tissue:=tis];ncounts[[tis]]<-nc
 quality[[tis]]<-data.table(tissue=tis,genes_input=nrow(raw),genes_kept=sum(keep),samples=ncol(dds),beta_nonconverged=sum(!mcols(dds)$betaConv,na.rm=TRUE),raw_counts_identical=TRUE,full_rank=qr(X)$rank,baseline_df=2,response_df=4)
 sizes[[tis]]<-data.table(tissue=tis,sample=colnames(dds),size_factor=sizeFactors(dds),library_total=colSums(mat),genes_detected=colSums(mat>0))
 metas[[tis]]<-as.data.table(md)
 co<-as.data.table(coef(dds),keep.rownames='gene_id');co[,tissue:=tis];coeffs[[tis]]<-co
 cat('DONE',tis,format(Sys.time()),'\n');flush.console()
}
om<-rbindlist(omnibus);om[,q_across_tissues:=p.adjust(pvalue,'BH'),by=test]
om[,q_strict:=NA_real_];om[expressed_all_groups==TRUE,q_strict:=p.adjust(pvalue,'BH'),by=.(tissue,test)]
om<-merge(om,annot,by='gene_id',all.x=TRUE)
co<-rbindlist(cons);co[,q_family:=p.adjust(pvalue,'BH'),by=.(tissue,family)]
co<-merge(co,annot,by='gene_id',all.x=TRUE)
fwrite(om,file.path(out,'omnibus_tests.csv.gz'));fwrite(co,file.path(out,'all_contrasts.csv.gz'))
fwrite(rbindlist(ncounts),file.path(out,'normalized_counts_long.csv.gz'));fwrite(rbindlist(metas),file.path(out,'sample_metadata.csv'))
fwrite(rbindlist(quality),file.path(out,'model_QA.csv'));fwrite(rbindlist(filters),file.path(out,'gene_filter_audit.csv.gz'));fwrite(rbindlist(sizes),file.path(out,'sample_QC.csv'))
fwrite(rbindlist(coeffs),file.path(out,'model_coefficients.csv.gz'))
fwrite(om[,.(genes=.N,finite_p=sum(is.finite(pvalue)),q05=sum(q_within<.05,na.rm=TRUE),q05_across_tissues=sum(q_across_tissues<.05,na.rm=TRUE),strict_genes=sum(expressed_all_groups),strict_q05=sum(q_strict<.05,na.rm=TRUE)),by=.(tissue,test)],file.path(out,'omnibus_summary.csv'))
fwrite(co[,.(genes=.N,finite_p=sum(is.finite(pvalue)),q05=sum(q_within<.05,na.rm=TRUE),q05_family=sum(q_family<.05,na.rm=TRUE),zero_conventions=sum(zero_count_convention)),by=.(tissue,test,family)],file.path(out,'contrast_summary.csv'))
fwrite(om[gene_id=='127522978'],file.path(out,'IL10_omnibus.csv'));fwrite(co[gene_id=='127522978'],file.path(out,'IL10_contrasts.csv'))
capture.output(sessionInfo(),file=file.path(out,'analysis_sessionInfo.txt'))
cat('ALL ANALYSES COMPLETE\n')
