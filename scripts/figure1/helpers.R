# Base-R helpers shared by Figure 1 entry points and unit tests.
decode_samples <- function(ids) {
  ids <- as.character(ids)
  raw <- grepl("^[ABC][CP][015]_[123]$", ids)
  if (all(raw)) {
    genotype <- unname(c(A="HomoIns", B="HomoDel", C="Hetero")[substr(ids,1,1)])
    tissue <- unname(c(P="Spleen", C="Intestine")[substr(ids,2,2)])
    day <- as.integer(substr(ids,3,3)); replicate <- as.integer(substr(ids,5,5))
  } else {
    fields <- strsplit(ids, "_", fixed=TRUE)
    if (!all(lengths(fields)==4)) stop("Unrecognized sample identifiers")
    genotype <- vapply(fields, `[`, "", 1)
    tissue <- vapply(fields, `[`, "", 2)
    tissue <- unname(c(Sp="Spleen", I="Intestine", Spleen="Spleen", Intestine="Intestine")[tissue])
    day <- as.integer(sub("d$", "", vapply(fields, `[`, "", 3)))
    replicate <- as.integer(sub("^Rep", "", vapply(fields, `[`, "", 4)))
  }
  stopifnot(!anyNA(genotype),!anyNA(tissue),all(day %in% c(0,1,5)),all(replicate %in% 1:3))
  data.frame(sample=ids, genotype=genotype, tissue=tissue, day=day, replicate=replicate)
}
standardize_genotype_means <- function(values) {
  stopifnot(length(values)==3,all(is.finite(values)))
  s <- sd(values)
  if (s==0) return(rep(0,3))
  (values-mean(values))/s
}
count_directional_degs <- function(log2fc, padj, numerator, denominator, target_num, target_den) {
  if (numerator==target_num && denominator==target_den) effect <- log2fc
  else if (numerator==target_den && denominator==target_num) effect <- -log2fc
  else stop("Contrast direction does not match the requested comparison")
  usable <- is.finite(effect) & is.finite(padj)
  c(Up=sum(usable & effect>1 & padj<0.05), Down=sum(usable & effect< -1 & padj<0.05))
}
