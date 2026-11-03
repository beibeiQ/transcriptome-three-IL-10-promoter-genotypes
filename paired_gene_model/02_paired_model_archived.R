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

suppressPackageStartupMessages({library(data.table);library(edgeR);library(limma)})
args<-commandArgs(TRUE);base<-if(length(args))args[1]else project_root()
out<-project_path("gene_contrasts");inp<-file.path(out,'inputs')
raw<-fread(file.path(inp,'raw_counts_input.tsv'));raw[,gene_id:=as.character(gene_id)]
md<-as.data.frame(fread(file.path(inp,'sample_metadata.csv')))
md$genotype<-factor(md$genotype,levels=c('HomoDel','Hetero','HomoIns'));md$time<-factor(md$time,levels=c('d0','d1','d5'));md$tissue<-factor(md$tissue,levels=c('Intestine','Spleen'))
stopifnot(nrow(md)==54,uniqueN(md$fish_id)==27,all(table(md$fish_id,md$tissue)==1),all(table(md$genotype,md$time,md$tissue)==3))
for(i in seq_len(nrow(md)))stopifnot(md$fish_id[i]==paste0(substr(md$sample[i],1,1),md$time[i],'_',md$rep[i]))
m<-as.matrix(raw[,md$sample,with=FALSE]);rownames(m)<-raw$gene_id
keep<-rowSums(m[,md$tissue=='Spleen']>=10)>=3&rowSums(m[,md$tissue=='Intestine']>=10)>=3
X<-model.matrix(~tissue*genotype*time,md);stopifnot(qr(X)$rank==18)
y<-calcNormFactors(DGEList(m[keep,]),method='TMM');v<-voom(y,X,plot=FALSE)
cr1<-duplicateCorrelation(v,X,block=md$fish_id);v<-voom(y,X,plot=FALSE,block=md$fish_id,correlation=cr1$consensus)
cr2<-duplicateCorrelation(v,X,block=md$fish_id);fit<-lmFit(v,X,block=md$fish_id,correlation=cr2$consensus)
rd<-function(s,g,t)as.numeric(model.matrix(~tissue*genotype*time,data.frame(tissue=factor(s,levels=levels(md$tissue)),genotype=factor(g,levels=levels(md$genotype)),time=factor(t,levels=levels(md$time)))))
delta<-function(s,g,t)rd(s,'Hetero',t)-rd(s,g,t)-rd(s,'Hetero','d0')+rd(s,g,'d0')
cts<-list();ds<-list()
add<-function(nm,vec,tissue,stage,comparator,family){cts[[nm]]<<-vec;ds[[nm]]<<-data.table(test=nm,tissue=tissue,stage=stage,comparator=comparator,family=family)}
for(g in c('HomoDel','HomoIns')){
 for(s in c('Spleen','Intestine')){
  add(paste(s,'Baseline',g,sep='_'),rd(s,'Hetero','d0')-rd(s,g,'d0'),s,'Baseline',g,'Joint_baseline')
  for(t in c('d1','d5'))add(paste(s,'Response',g,t,sep='_'),delta(s,g,t),s,paste0('Response_',t),g,'Joint_response')
 }
 b<-(rd('Spleen','Hetero','d0')-rd('Spleen',g,'d0'))-(rd('Intestine','Hetero','d0')-rd('Intestine',g,'d0'))
 add(paste('Tissue_Baseline',g,sep='_'),b,'Spleen_minus_Intestine','Baseline',g,'Tissue_baseline')
 for(t in c('d1','d5'))add(paste('Tissue_Response',g,t,sep='_'),delta('Spleen',g,t)-delta('Intestine',g,t),'Spleen_minus_Intestine',paste0('Response_',t),g,'Tissue_response')
}
C<-do.call(cbind,cts);rownames(C)<-colnames(X);cf<-eBayes(contrasts.fit(fit,C),robust=TRUE)
res<-rbindlist(lapply(seq_len(ncol(C)),function(k){
 z<-as.data.table(topTable(cf,coef=k,number=Inf,sort.by='none',confint=.95),keep.rownames='gene_id')
 setnames(z,c('logFC','P.Value','adj.P.Val','CI.L','CI.R'),c('effect','pvalue','q_within','ci95_low','ci95_high'))
 z[,`:=`(test=colnames(C)[k],SE=cf$stdev.unscaled[,k]*sqrt(cf$s2.post),df_moderated=cf$df.total,tcrit=qt(.975,cf$df.total),linear_effect=as.numeric(fit$coefficients%*%C[,k]))]
 z
}))
res<-merge(res,rbindlist(ds),by='test');res[,q_family:=p.adjust(pvalue,'BH'),by=family]
stopifnot(all(abs(res$effect-res$linear_effect)<1e-8),all(abs(res$t-res$effect/res$SE)<1e-8))
fwrite(res,file.path(out,'paired_all_contrasts.csv.gz'))
defs<-rbindlist(lapply(names(cts),function(n)data.table(test=n,coefficient=colnames(X),weight=cts[[n]])))
fwrite(defs,file.path(out,'paired_contrast_definitions.csv'));fwrite(as.data.table(fit$coefficients,keep.rownames='gene_id'),file.path(out,'paired_coefficients.csv.gz'))
idx<-which(grepl('tissueSpleen:genotype.*:time',colnames(X)));stopifnot(length(idx)==4)
om<-as.data.table(topTable(eBayes(fit,robust=TRUE),coef=idx,number=Inf,sort.by='none'),keep.rownames='gene_id');fwrite(om,file.path(out,'threeway_omnibus.csv.gz'))
fwrite(as.data.table(v$E,keep.rownames='gene_id'),file.path(out,'paired_voom_logCPM.csv.gz'))
fwrite(data.table(sample=md$sample,tmm_factor=y$samples$norm.factors,library_size=y$samples$lib.size),file.path(out,'paired_normalization.csv'))
fwrite(as.data.table(md),file.path(out,'sample_metadata.csv'))
fwrite(data.table(genes_input=nrow(raw),genes_common=sum(keep),libraries=ncol(v$E),fish=uniqueN(md$fish_id),rank=qr(X)$rank,nominal_residual_df=54-qr(X)$rank,correlation_pass1=cr1$consensus,correlation_pass2=cr2$consensus,threeway_q05=sum(om$adj.P.Val<.05)),file.path(out,'paired_model_QA.csv'))
fwrite(res[,.(genes=.N,q05=sum(q_family<.05),min_q=min(q_family)),by=.(family,tissue,stage,comparator)],file.path(out,'paired_summary.csv'))
saveRDS(list(metadata=md,voom=v,fit=fit,contrast_fit=cf,contrast_matrix=C,correlation=cr2$consensus),file.path(out,'paired_model.rds'))
capture.output(sessionInfo(),file=file.path(out,'paired_sessionInfo.txt'))
cat('PAIRED COMPLETE',sum(keep),'genes, correlation',cr2$consensus,'\n');print(res[family=='Tissue_response',.(genes=.N,q05=sum(q_family<.05)),by=.(stage,comparator)])
