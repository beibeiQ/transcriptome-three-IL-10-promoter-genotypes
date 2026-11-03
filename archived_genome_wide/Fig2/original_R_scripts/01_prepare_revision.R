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
args<-commandArgs(TRUE);base<-if(length(args))args[1]else project_root()
old<-project_path("figure2_analysis");out<-project_path("figure2_revision");inp<-file.path(out,'inputs');dir.create(inp,recursive=TRUE,showWarnings=FALSE)
rgz<-function(p)fread(cmd=paste('gzip -dc',shQuote(p)),colClasses=c(gene_id='character'))
for(n in c('omnibus_tests.csv.gz','model_coefficients.csv.gz','Figure_legends.txt'))if(!file.exists(file.path(inp,n)))stopifnot(file.copy(file.path(old,n),file.path(inp,n)))
for(n in c('Figure2A_counts.csv','Figure2C_points.csv','Figure2C_summary.csv','Figure2C_gene_selection.csv','Figure2C_model_effects.csv','Figure2D_selected_pathways.csv'))if(!file.exists(file.path(out,n)))stopifnot(file.copy(file.path(old,n),file.path(out,n)))
om<-rgz(file.path(inp,'omnibus_tests.csv.gz'));target<-om[test=='Response' & q_within<.05,.(tissue,gene_id,response_q=q_within)]
src<-file.path(inp,'selected_contrasts.csv.gz')
if(!file.exists(src)){
 co<-rgz(file.path(old,'all_contrasts.csv.gz'));co<-merge(co[family%in%c('Within_genotype','Response_pairwise')],target,by=c('tissue','gene_id'));fwrite(co,src)
}else co<-rgz(src)
geno<-c('HomoDel','Hetero','HomoIns');days<-c('d1','d5');keys<-as.vector(outer(geno,days,paste,sep='_'));stopifnot(identical(keys,c('HomoDel_d1','Hetero_d1','HomoIns_d1','HomoDel_d5','Hetero_d5','HomoIns_d5')))
main<-co[family=='Within_genotype'];main[,column:=paste(g1,time,sep='_')]
main[,display_available:=is.finite(log2FoldChange)&!zero_count_convention];main[,effect_display:=fifelse(display_available,log2FoldChange,NA_real_)]
heat<-list();orders<-list();qa<-list();centres<-list();invariance<-list()
set.seed(20261009)
for(tis in c('Spleen','Intestine')){
 z<-main[tissue==tis];wide<-dcast(z,gene_id~column,value.var='linear_predictor_log2FC');M<-as.matrix(wide[,..keys]);rownames(M)<-wide$gene_id
 stopifnot(nrow(M)==nrow(target[tissue==tis]),all(is.finite(M)))
 C<-M;for(j in list(1:3,4:6))C[,j]<-sweep(M[,j,drop=FALSE],1,rowMeans(M[,j,drop=FALSE]),'-')
 # Alternative genotype references give the exact same centred ordering features.
 diffs<-sapply(c(1,2,3),function(ref){Q<-M;for(j in list(1:3,4:6)){Q[,j]<-sweep(M[,j,drop=FALSE],1,M[,j[ref]],'-');Q[,j]<-sweep(Q[,j,drop=FALSE],1,rowMeans(Q[,j,drop=FALSE]),'-')};max(abs(Q-C))})
 stopifnot(max(diffs)<1e-10);invariance[[tis]]<-data.table(tissue=tis,reference=geno,max_abs_feature_difference=diffs)
 km<-kmeans(C,centers=4,nstart=30,iter.max=100);rank<-do.call(order,as.data.frame(km$centers));cl<-match(km$cluster,rank);ix<-order(cl,C[,1],C[,2],rownames(C))
 ord<-data.table(gene_id=rownames(M)[ix],row_index=seq_along(ix),cluster=cl[ix],tissue=tis);orders[[tis]]<-ord
 hh<-merge(z,ord,by=c('tissue','gene_id'));hh[,color_value:=pmax(-4,pmin(4,effect_display))];heat[[tis]]<-hh
 feature<-as.data.table(C,keep.rownames='gene_id');feature[,tissue:=tis];centres[[tis]]<-feature
 qa[[tis]]<-data.table(tissue=tis,selected_genes=nrow(M),shown_genes=uniqueN(hh$gene_id),cells=nrow(hh),grey_cells=sum(!hh$display_available),saturated_cells=sum(abs(hh$effect_display)>4,na.rm=TRUE),min_effect=min(hh$effect_display,na.rm=TRUE),max_effect=max(hh$effect_display,na.rm=TRUE))
}
heat<-rbindlist(heat);ord<-rbindlist(orders);fwrite(heat,file.path(out,'Figure2B_own_baseline_source.csv.gz'));fwrite(ord,file.path(out,'Figure2B_row_order.csv.gz'));fwrite(rbindlist(centres),file.path(out,'Figure2B_ordering_features.csv.gz'));fwrite(rbindlist(invariance),file.path(out,'reference_invariance_QA.csv'));fwrite(rbindlist(qa),file.path(out,'Figure2B_display_QA.csv'))
pa<-co[family=='Response_pairwise'];pa[,comparison:=paste(g1,g0,sep='_vs_')]
pa<-merge(pa,ord,by=c('tissue','gene_id'));pa[,display_available:=is.finite(log2FoldChange)&!zero_count_convention];pa[,effect_display:=fifelse(display_available,log2FoldChange,NA_real_)];pa[,color_value:=pmax(-4,pmin(4,effect_display))]
fwrite(pa,file.path(out,'FigureS2_all_pairwise_source.csv.gz'));fwrite(pa[,.(genes=uniqueN(gene_id),cells=.N,grey_cells=sum(!display_available),saturated_cells=sum(abs(effect_display)>4,na.rm=TRUE)),by=tissue],file.path(out,'FigureS2_pairwise_display_QA.csv'))
fwrite(data.table(file=list.files(inp),md5=unname(tools::md5sum(list.files(inp,full.names=TRUE)))),file.path(out,'input_checksums.csv'));capture.output(sessionInfo(),file=file.path(out,'preparation_sessionInfo.txt'));print(rbindlist(qa))
