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

# ============================================================
# Paired-tissue functional-context analysis
# Inputs: FPKM, KEGG membership and sample metadata.
# Outputs: complete score effects, diagnostics and reproducible RDS
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[1], winslash = '/')
wd <- file.path(root, 'work/tissue_function_revision')
.libPaths(c(file.path(root, 'work/Rlib'), .libPaths()))
suppressPackageStartupMessages({library(GSVA); library(data.table); library(BiocParallel)})
set.seed(42)
inp <- file.path(wd, 'inputs')
out <- file.path(wd, 'analysis')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
meta <- as.data.frame(fread(file.path(inp, 'sample_metadata.csv')))
meta$genotype <- factor(meta$genotype, levels = c('HomoDel', 'Hetero', 'HomoIns'))
meta$time <- factor(meta$time, levels = c('d0', 'd1', 'd5'))
raw <- fread(file.path(inp, 'paired_FPKM.csv.gz'))
mat <- as.matrix(raw[, meta$sample, with = FALSE]); rownames(mat) <- raw$gene_id
stopifnot(ncol(mat) == 54, !anyNA(mat), all(is.finite(mat)), all(mat >= 0))
gmt <- strsplit(readLines(file.path(inp, 'pathways.gmt')), '\t', fixed = TRUE)
sets <- setNames(lapply(gmt, function(v) v[-c(1,2)]), sapply(gmt, '[', 1))
score_fun <- function(x, sets) {
  param <- gsvaParam(log2(x + 1), sets, kcdf = 'Gaussian', minSize = 10,
                     maxSize = Inf, tau = 1, maxDiff = TRUE,
                     absRanking = FALSE, sparse = FALSE, ondisk = 'no')
  gsva(param, verbose = FALSE, BPPARAM = SerialParam())
}
scores <- score_fun(mat, sets)
stopifnot(identical(colnames(scores), meta$sample), nrow(scores) == 8)
fwrite(data.table(pathway_id = rownames(scores), scores), file.path(out, 'pooled_GSVA_scores.csv'))
long <- melt(data.table(pathway_id = rownames(scores), scores), id.vars = 'pathway_id', variable.name = 'sample', value.name = 'score')
long <- merge(long, as.data.table(meta), by = 'sample')
fwrite(long, file.path(out, 'individual_scores.csv'))
fwrite(long[, .(n = .N, mean = mean(score), SD = sd(score)), by = .(pathway_id, tissue, genotype, time)], file.path(out, 'group_score_summary.csv'))

make_design <- function(d) {
  d$group <- factor(paste(d$genotype, d$time, sep = '_'),
                    levels = as.vector(outer(c('HomoDel','Hetero','HomoIns'), c('d0','d1','d5'), paste, sep = '_')))
  x <- model.matrix(~0 + group, d)
  colnames(x) <- levels(d$group)
  x
}
make_contrasts <- function(x) {
  rv <- function(g,t) {v <- rep(0, ncol(x)); names(v) <- colnames(x); v[paste(g,t,sep='_')] <- 1; v}
  cc <- list(); info <- list()
  for (g in c('HomoDel','HomoIns')) {
    for (t in c('d0','d1','d5')) {
      nm <- paste('Endpoint', g, t, sep='|')
      cc[[nm]] <- rv('Hetero', t)-rv(g,t)
      info[[nm]] <- data.table(test=nm, question='Endpoint', comparator=g, time=t)
    }
    for (t in c('d1','d5')) {
      nm <- paste('Response', g, t, sep='|')
      cc[[nm]] <- rv('Hetero',t)-rv('Hetero','d0')-rv(g,t)+rv(g,'d0')
      info[[nm]] <- data.table(test=nm, question='Response', comparator=g, time=t)
    }
  }
  list(matrix = do.call(cbind,cc), info = rbindlist(info))
}
fit_scores <- function(es, d, label, loo=TRUE) {
  x <- make_design(d); cm <- make_contrasts(x); cmat <- cm$matrix
  y <- t(es); fit <- lm.fit(x, y); df <- nrow(x)-qr(x)$rank
  stopifnot(df == 18, qr(x)$rank==9)
  sigma <- sqrt(colSums(fit$residuals^2)/df)
  covx <- solve(crossprod(x)); se0 <- sqrt(diag(t(cmat)%*%covx%*%cmat))
  effects <- t(cmat)%*%fit$coefficients
  se <- outer(se0, sigma); tv <- effects/se
  effects_long <- rbindlist(lapply(seq_len(nrow(effects)),function(k) {
    data.table(pathway_id=rownames(es),test=rownames(effects)[k],estimate=effects[k,],SE=se[k,],df=df,
               CI_low=effects[k,]-qt(.975,df)*se[k,],CI_high=effects[k,]+qt(.975,df)*se[k,],P=2*pt(-abs(tv[k,]),df))
  }))
  effects_long <- merge(effects_long,cm$info,by='test'); effects_long[,tissue:=label]
  # Omnibus genotype x time vs additive genotype + time.
  xa <- model.matrix(~genotype+time,d); fa <- lm.fit(xa,y)
  rss_full <- colSums(fit$residuals^2); rss_add <- colSums(fa$residuals^2)
  fval <- pmax(0,(rss_add-rss_full)/4)/(rss_full/df)
  omnibus <- data.table(pathway_id=rownames(es),tissue=label,F=fval,df1=4,df2=df,P=pf(fval,4,df,lower.tail=FALSE))
  leverage <- rowSums((x%*%covx)*x)
  student <- sweep(fit$residuals,2,sigma,'/')/sqrt(1-leverage)
  cook <- student^2*leverage/(ncol(x)*(1-leverage))
  diag <- rbindlist(lapply(seq_len(ncol(y)),function(j) data.table(pathway_id=rownames(es)[j],tissue=label,sample=d$sample,fish_id=d$fish_id,
           genotype=d$genotype,time=d$time,fitted=fit$fitted.values[,j],residual=fit$residuals[,j],internally_studentized=student[,j],leverage=leverage,Cooks_distance=cook[,j])))
  influence <- NULL
  if (loo) {
    shifts <- lapply(seq_len(nrow(d)),function(i) t(cmat)%*%lm.fit(x[-i,,drop=FALSE], y[-i,,drop=FALSE])$coefficients)
    influence <- rbindlist(lapply(seq_len(nrow(effects)),function(k) {
      vals <- sapply(shifts,function(z) z[k,])
      data.table(pathway_id=rownames(es),test=rownames(effects)[k],tissue=label,
       loo_min=apply(vals,1,min),loo_max=apply(vals,1,max),
       maximum_absolute_shift=apply(abs(vals-effects[k,]),1,max),
       sign_retention=rowMeans(sign(vals)==sign(effects[k,])))
    }))
  }
  list(effects=effects_long,omnibus=omnibus,diagnostics=diag,influence=influence,
       model=list(design=x,contrast=cm,fit=fit,sigma=sigma,df=df))
}

fits <- lapply(c('Spleen','Intestine'),function(tis) {
  idx <- which(meta$tissue==tis); fit_scores(scores[,idx,drop=FALSE],meta[idx,],tis)
}); names(fits) <- c('Spleen','Intestine')
effects <- rbindlist(lapply(fits,'[[','effects')); effects[,Q:=p.adjust(P,'BH'),by=question]
omnibus <- rbindlist(lapply(fits,'[[','omnibus')); omnibus[,Q:=p.adjust(P,'BH')]
fwrite(effects,file.path(out,'program_contrasts.csv'))
fwrite(omnibus,file.path(out,'program_interactions.csv'))
fwrite(rbindlist(lapply(fits,'[[','diagnostics')),file.path(out,'score_model_diagnostics.csv'))
fwrite(rbindlist(lapply(fits,'[[','influence')),file.path(out,'score_leave_one_fish_ranges.csv'))

# Pair by verified metadata fish IDs; never infer pairing from sample position.
sp <- meta[meta$tissue=='Spleen',]; it <- meta[meta$tissue=='Intestine',]
it <- it[match(sp$fish_id,it$fish_id),]
stopifnot(identical(as.character(sp$fish_id),as.character(it$fish_id)),
          identical(as.character(sp$genotype),as.character(it$genotype)),
          identical(as.character(sp$time),as.character(it$time)))
paired <- scores[,sp$sample,drop=FALSE]-scores[,it$sample,drop=FALSE]
colnames(paired)<-sp$sample
direct <- fit_scores(paired,sp,'Spleen_minus_Intestine')
de <- direct$effects[question=='Response' | time=='d0']
de[,question:=ifelse(question=='Endpoint','Baseline','Response')]
de[,Q:=p.adjust(P,'BH'),by=question]
stopifnot(nrow(de[question=='Response'])==32,nrow(de[question=='Baseline'])==16)
fwrite(de,file.path(out,'direct_tissue_program_contrasts.csv'))
dp <- melt(data.table(pathway_id=rownames(paired),paired),id.vars='pathway_id',variable.name='sample',value.name='paired_difference')
dp <- merge(dp,as.data.table(sp),by='sample')
fwrite(dp,file.path(out,'paired_program_differences.csv'))
fwrite(direct$diagnostics,file.path(out,'direct_score_model_diagnostics.csv'))
fwrite(direct$influence,file.path(out,'direct_score_leave_one_fish_ranges.csv'))

# Sensitivity: separate tissue scoring, keeping fixed pathway definitions.
sens_fits<-list(); coverage<-list(); sensitivity_scores<-list()
for (tis in c('Spleen','Intestine')) {
  idx<-which(meta$tissue==tis); x<-mat[,idx,drop=FALSE]; keep<-apply(x,1,var)>0
  sm<-score_fun(x[keep,,drop=FALSE],sets)
  sensitivity_scores[[tis]]<-sm
  sens_fits[[tis]]<-fit_scores(sm,meta[idx,],tis,loo=FALSE)
  coverage[[tis]]<-data.table(tissue=tis,pathway_id=names(sets),scored_genes=sapply(sets,function(v)sum(v%in%rownames(x)[keep])),gene_universe=sum(keep))
}
sens<-rbindlist(lapply(sens_fits,'[[','effects'))[question=='Response'];sens[,Q:=p.adjust(P,'BH')]
fwrite(sens,file.path(out,'within_tissue_scoring_sensitivity.csv'))
fwrite(rbindlist(coverage),file.path(out,'sensitivity_mapping_coverage.csv'))
joined<-merge(effects[question=='Response',.(pathway_id,tissue,comparator,time,pooled_effect=estimate,pooled_Q=Q)],
              sens[,.(pathway_id,tissue,comparator,time,separate_effect=estimate,separate_Q=Q)],by=c('pathway_id','tissue','comparator','time'))
joined[,same_direction:=sign(pooled_effect)==sign(separate_effect)]
fwrite(joined,file.path(out,'scoring_sensitivity_comparison.csv'))
saveRDS(list(meta=meta,expression=mat,sets=sets,scores=scores,within_tissue=fits,
             paired_scores=paired,direct=direct,sensitivity=sens_fits,
             sensitivity_scores=sensitivity_scores,complete_effects=effects),file.path(out,'tissue_function_analysis.rds'))
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
print(omnibus); print(effects[question=='Response' & Q<.05]); print(de[question=='Response' & Q<.05])
print(joined[,.(n=.N,same_direction=sum(same_direction),primary_detected=sum(pooled_Q<.05),both_detected=sum(pooled_Q<.05 & separate_Q<.05))])
