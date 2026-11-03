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

args<-commandArgs(trailingOnly=TRUE)
root<-normalizePath(args[1],winslash='/'); wd<-file.path(root,'work/macrophage_revision')
.libPaths(c(file.path(root,'work/Rlib'),.libPaths()))
suppressPackageStartupMessages({library(GSVA);library(data.table);library(limma);library(BiocParallel)})
meta<-as.data.frame(fread(file.path(wd,'sample_metadata.csv')));meta<-meta[meta$tissue=='Spleen',]
meta$genotype<-factor(meta$genotype,levels=c('HomoDel','Hetero','HomoIns'));meta$time<-factor(meta$time,levels=c('d0','d1','d5'))
raw<-fread(file.path(wd,'Spleen_FPKM.csv.gz'));mat<-as.matrix(raw[,meta$sample,with=FALSE]);rownames(mat)<-raw$gene_id
mat<-mat[apply(mat,1,var)>0,,drop=FALSE]
gmt<-strsplit(readLines(file.path(wd,'marker_sets.gmt')),'\t',fixed=TRUE)
sets<-setNames(lapply(gmt,function(v)v[-c(1,2)]),sapply(gmt,'[',1))
sets<-lapply(sets,function(x)intersect(x,rownames(mat)))
param<-gsvaParam(log2(mat+1),sets,kcdf='Gaussian',minSize=3,maxSize=Inf,tau=1,maxDiff=TRUE,absRanking=FALSE,sparse=FALSE,ondisk='no')
es<-gsva(param,verbose=FALSE,BPPARAM=SerialParam())
fwrite(data.table(signature=rownames(es),es),file.path(wd,'signature_scores.csv'))
oldfile<-file.path(wd,'archived_marker_scores.csv')
old<-fread(oldfile);oldm<-as.matrix(old[,meta$sample,with=FALSE]);rownames(oldm)<-old$signature
repro<-data.table(signature=rownames(oldm),max_abs_difference=sapply(rownames(oldm),function(n)max(abs(es[n,]-oldm[n,]))))
fwrite(repro,file.path(wd,'archived_score_reproduction.csv'))
design<-model.matrix(~genotype*time,meta)
rowdesign<-function(g,t)as.numeric(model.matrix(~genotype*time,data.frame(genotype=factor(g,levels=levels(meta$genotype)),time=factor(t,levels=levels(meta$time)))))
primary<-c('Macrophage','M1','M2');sensitivity<-c('M1_without_shared','M2_without_shared')
omnibus<-rbindlist(lapply(rownames(es),function(n){d<-meta;d$score<-as.numeric(es[n,]);a<-anova(lm(score~genotype+time,d),lm(score~genotype*time,d));data.table(signature=n,F=a$F[2],df1=a$Df[2],df2=a$Res.Df[2],P=a$`Pr(>F)`[2])}))
omnibus[,family:=ifelse(signature %in% primary,'primary','shared-gene sensitivity')]
omnibus[,Q:=p.adjust(P,'BH'),by=family]
fwrite(omnibus,file.path(wd,'signature_interaction_tests.csv'))
contr<-list();info<-list()
for(g in c('HomoDel','HomoIns')){
 for(t in c('d0','d1','d5')){
   nm<-paste('Endpoint',g,t,sep='|');contr[[nm]]<-rowdesign('Hetero',t)-rowdesign(g,t);info[[nm]]<-data.table(test=nm,question='Endpoint',comparator=g,time=t)
 }
 for(t in c('d1','d5')){
   nm<-paste('Response',g,t,sep='|');contr[[nm]]<-(rowdesign('Hetero',t)-rowdesign('Hetero','d0'))-(rowdesign(g,t)-rowdesign(g,'d0'));info[[nm]]<-data.table(test=nm,question='Response',comparator=g,time=t)
 }
}
cm<-do.call(cbind,contr);rownames(cm)<-colnames(design);fit<-contrasts.fit(lmFit(es,design),cm)
res<-rbindlist(lapply(seq_along(contr),function(k){est<-fit$coefficients[,k];se<-fit$stdev.unscaled[,k]*fit$sigma;df<-fit$df.residual;data.table(signature=rownames(es),test=names(contr)[k],estimate=est,SE=se,df=df,CI_low=est-qt(.975,df)*se,CI_high=est+qt(.975,df)*se,P=2*pt(-abs(est/se),df))}))
res<-merge(res,rbindlist(info),by='test');res[,membership_family:=ifelse(signature %in% primary,'primary','shared-gene sensitivity')]
res[,Q:=p.adjust(P,'BH'),by=.(question,membership_family)]
fwrite(res,file.path(wd,'signature_contrasts.csv'))
long<-melt(data.table(signature=rownames(es),es),id.vars='signature',variable.name='sample',value.name='score');long<-merge(long,meta,by='sample')
fwrite(long,file.path(wd,'signature_individual_scores.csv'))
summ<-long[,.(n=.N,mean=mean(score),SD=sd(score)),by=.(signature,genotype,time)]
fwrite(summ,file.path(wd,'signature_group_summary.csv'))
me<-fread(file.path(wd,'ME1_individual_counts.csv'));me<-me[tissue=='Spleen'];me<-me[match(meta$sample,me$sample),]
y<-log2(me$normalized_count+1);resy<-lm.fit(design,y)$residuals
corres<-rbindlist(lapply(primary,function(n){s<-as.numeric(es[n,]);ress<-lm.fit(design,s)$residuals;r<-cor(resy,ress);df<-nrow(meta)-qr(design)$rank-1;t<-r*sqrt(df/(1-r*r));data.table(signature=n,partial_r=r,df=df,P=2*pt(-abs(t),df))}))
corres[,Q:=p.adjust(P,'BH')];fwrite(corres,file.path(wd,'ME1_adjusted_score_correlations.csv'))
capture.output(sessionInfo(),file=file.path(wd,'R_sessionInfo.txt'))
print(repro);print(omnibus);print(summ[signature %in% primary]);print(res[membership_family=='primary' & time=='d1']);print(corres)
