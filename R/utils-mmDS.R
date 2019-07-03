#' @describeIn mmDS
#'
#' see details.
#'
#' @param dup_corr logical; whether to use
#'   \code{\link[limma]{duplicateCorrelation}}.
#' @param trended logical; whether to use expression-dependent variance priors
#'  in \code{\link[limma]{eBayes}}.
#' @param ddf character string specifying the method for estimating
#'  the effective degrees of freedom. For \code{method = "dream"},
#'  either \code{"Satterthwaite"} (faster) or \code{"Kenward-Roger"}
#'  (more accurate); see \code{\link[variancePartition]{dream}}.
#'  For \code{method = "vst"}, method \code{"lme4"} is also valid;
#'  see \code{\link[lmerTest]{contest.lmerModLmerTest}}.
#'
#' @details
#' \code{.mm_dream} and \code{.mm_vst} expect cells from a single cluster,
#' and do not perform filtering or handle incorrect parameters well.
#' Meant to be called by \code{mmDS} with \code{method = c("dream", "vst")} and
#' \code{vst = c("sctransform", "DESeq2")} to be applied across all clusters.
#' \describe{
#' \item{\code{method = "dream"}}{
#'   voom-lme4-implementation \code{\link[variancePartition]{dream}}
#'   of mixed models for RNAseq data.}
#' \item{\code{method = "vst"}}{
#'   \describe{
#'   \item{\code{vst = "sctransform"}}{
#'     \code{lmer} or \code{blmer} mixed models on
#'     \code{\link[sctransform]{vst}} transformed counts.}
#'   \item{\code{vst = "DESeq2"}}{
#'     \code{\link[DESeq2]{varianceStabilizingTransformation}}
#'     followed by \code{lme4} mixed models.}}}}
#'
#' @importFrom edgeR DGEList
#' @importFrom doParallel registerDoParallel
#' @importFrom dplyr %>% last mutate_at rename
#' @importFrom limma duplicateCorrelation eBayes topTable voom
#' @importFrom magrittr set_rownames
#' @importFrom matrixStats rowSds
#' @importFrom parallel makeCluster stopCluster
#' @importFrom scran computeSumFactors
#' @importFrom SingleCellExperiment counts sizeFactors
#' @importFrom stats as.formula model.matrix
#' @importFrom variancePartition dream getContrast
.mm_dream <- function(x,
    coef, covs, n_threads, verbose,
    dup_corr = FALSE, trended = FALSE,
    ddf = c("Satterthwaite", "Kenward-Roger")) {
    
    if (is.null(sizeFactors(x)))
        x <- computeSumFactors(x)

    ddf <- match.arg(ddf)
    x <- x[rowSds(as.matrix(counts(x))) > 0, ]
    y <- DGEList(counts(x), norm.factors = 1 / sizeFactors(x))

    cd <- .prep_cd(x, covs)

    formula <- paste0("~", paste(c(covs, "group_id"), collapse = "+"))
    mm <- model.matrix(as.formula(formula), data = cd)
    v <- voom(y, mm)

    if (dup_corr) {
        dup_corr <- duplicateCorrelation(v, mm, block = x$sample_id)
        v <- voom(y, mm, block = x$sample_id, correlation = dup_corr$consensus)
    }

    if (n_threads > 1) {
        cl <- makeCluster(n_threads)
        registerDoParallel(cl)
    }

    formula <- paste0(formula, "+(1|sample_id)")
    if (verbose) print(formula)

    if (is.null(coef)) {
        coef <- last(colnames(mm))
        if (verbose)
            message("Argument 'coef' not specified; ",
                sprintf("testing for %s.", dQuote(coef)))
    }

    contrast <- getContrast(v, as.formula(formula), cd, coef)
    fit <- dream(v, formula, cd, contrast, ddf = ddf, suppressWarnings = !verbose)
    fit <- eBayes(fit, trend = trended, robust = TRUE)
    if (n_threads > 1) stopCluster(cl)

    topTable(fit, number = Inf, sort.by = "none") %>%
        rename(p_val = "P.Value", p_adj.loc = "adj.P.Val")
}

#' @describeIn mmDS
#'
#' see details.
#'
#' @param vst method to use as variance-stabilizing transformations.
#'   \code{"sctransform"} for \code{\link[sctransform]{vst}};
#'   \code{"DESeq2"} for \code{\link[DESeq2]{varianceStabilizingTransformation}}.
#' @param bayesian logical; whether to use bayesian mixed models.
#' @param blind logical; whether to ignore experimental design for the vst.
#' @param REML logical; whether to maximize REML instead of log-likelihood.
#'
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateDispersions
#'   sizeFactors varianceStabilizingTransformation
#' @importFrom dplyr last
#' @importFrom purrr set_names
#' @importFrom sctransform vst
#' @importFrom scran computeSumFactors
#' @importFrom SingleCellExperiment counts sizeFactors sizeFactors<-
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble add_column
.mm_vst <- function(x,
    coef, covs, n_threads, verbose,
    vst = c("sctransform", "DESeq2"),
    bayesian = FALSE, blind = TRUE, REML = TRUE,
    ddf = c("Satterthwaite", "Kenward-Roger", "lme4")) {

    vst <- match.arg(vst)
    ddf <- match.arg(ddf)
    cd <- .prep_cd(x, covs)

    # variance-stabilizing transformation
    fun_call <- switch(vst,
        sctransform = {
            # assure correct vst() function is used
            fun <- getFromNamespace("vst", vst)
            expression(fun(assay(x), show_progress = verbose)$y)
        },
        DESeq2 = {
            if (is.null(sizeFactors(x)))
                x <- computeSumFactors(x)
            formula <- paste("~", paste(c(covs, "sample_id"), collapse="+"))
            formula <- as.formula(formula)
            y <- as.matrix(counts(x))
            y <- suppressMessages(DESeqDataSetFromMatrix(y, cd, formula))
            sizeFactors(y) <- sizeFactors(x)
            if (!blind) y <- estimateDispersions(y)
            expression(assay(varianceStabilizingTransformation(y, blind)))
        })
    if (verbose) y <- eval(fun_call) else y <- suppressMessages(eval(fun_call))

    # get formula
    formula <- paste(c("~(1|sample_id)", covs, "group_id"), collapse = "+")
    if (verbose) print(formula)
    formula <- as.formula(paste("u", formula))

    # get coefficient to test
    if (is.null(coef)) {
        coef <- paste0("group_id", last(levels(x$group_id)))
        if (verbose)
            message("Argument 'coef' not specified; ",
                sprintf("testing for %s.", dQuote(coef)))
    }

    # fit mixed models for ea. gene
    fits <- bplapply(seq_len(nrow(y)), function(i)
        .fit_lmer(cbind(u = y[i, ], cd), formula, coef, bayesian, REML, ddf),
        BPPARAM = MulticoreParam(n_threads, progressbar=verbose)) %>% 
        set_names(rownames(y))

    if (verbose) message("Applying empirical Bayes moderation..")
    .mm_eBayes(fits, coef) %>% 
        add_column(.after = "p_val", p_adj.loc = p.adjust(.$p_val))
}

# helper to prepare colData for .mm_dream/vst
#' @importFrom dplyr %>% mutate_at mutate_if
#' @importFrom magrittr set_rownames
#' @importFrom SummarizedExperiment colData
.prep_cd <- function(x, covs) {
    colData(x)[c("sample_id", "group_id", covs)] %>% 
        data.frame(check.names = FALSE) %>%
        mutate_if(is.factor, droplevels) %>% 
        mutate_at(covs, function(u) if (is.numeric(u)) scale(u)) %>%
        set_rownames(colnames(x))
}


#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr %>% bind_rows last
#' @importFrom purrr set_names
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble add_column
.mm_glmm <- function(x, coef, covs, n_threads, verbose=TRUE, moderate=FALSE, family=c("nbinom", "poisson")){
    family <- match.arg(family)
    cd <- .prep_cd(x, covs)
    y <- counts(x)
    if(is.null(sizeFactors(x))){
        cd$ls <- log(colSums(y))
    }else{
        cd$ls <- sizeFactors(x)
    }
    
    # get formula
    formula <- paste(c("~(1|sample_id)+offset(ls)", covs, "group_id"), collapse = "+")
    if (verbose) print(formula)
    formula <- as.formula(paste("u", formula))
    
    # get coefficient to test
    if (is.null(coef)) {
        coef <- paste0("group_id", last(levels(x$group_id)))
        if (verbose) 
            message("Argument 'coef' not specified; ", 
                    sprintf("testing for %s.", dQuote(coef)))
    }
    
    # fit mixed model for ea. gene
    fits <- bplapply(seq_len(nrow(y)), function(i) {
        if(family=="nbinom"){
            return(.fit_nbinom(df=data.frame(u=y[i, ], cd), formula, coef, prepForModeration=moderate))
        }
        return(.fit_bglmer(df=data.frame(u=y[i, ], cd), formula, coef, prepForModeration=moderate))
    }, BPPARAM = MulticoreParam(n_threads, progressbar=verbose)) %>% 
        set_names(rownames(y))
    
    if(moderate){
        if (verbose) message("Applying empirical Bayes moderation..")
        fits <- .mm_eBayes(fits, coef)
    }else{
        fits <- as.data.frame(t(bind_rows(fits)))
        colnames(fits) <- c("beta", "SE", "stat", "p_val")
    }
    fits %>% add_column(.after = "p_val", p_adj.loc = p.adjust(.$p_val))
}

#' @import blme
.mm_poisson <- function(x, coef, covs, n_threads, verbose=TRUE, ...){
    .mm_glmm(x, coef, covs, n_threads, verbose=verbose, family="poisson", ...)
}

#' @import glmmTMB
.mm_nbinom <- function(x, coef, covs, n_threads, verbose=TRUE, ...){
    .mm_glmm(x, coef, covs, n_threads, verbose=verbose, family="nbinom")
}



#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom dplyr %>% bind_rows last rename
#' @importFrom purrr set_names
#' @importFrom SingleCellExperiment counts SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble add_column
#' @importFrom Matrix.utils aggregate.Matrix
#' @import glmmTMB
#' @import blme
.mm_hybrid <- function(x, coef, covs, n_threads, verbose=TRUE, fam=c("nbinom","poisson"), pthres4mm=0.1){
    fam <- match.arg(fam)
    x$cluster_id <- droplevels(x$cluster_id)
    cd <- .prep_cd(x, covs)
    y <- counts(x)
    if(is.null(sizeFactors(x))){
        cd$ls <- log(colSums(y))
    }else{
        cd$ls <- sizeFactors(x)
    }
    
    pb <- SingleCellExperiment( list(t(aggregate.Matrix(t(counts(x)), x$sample_id))) %>% 
                                    set_names(as.character(x$cluster_id[1])) )
    pb.cd <- as.data.frame(getPBcolData(x)[colnames(pb),])
    
    # get sample-level formula
    f <- paste(c("~group_id", intersect(covs, colnames(pb.cd))), collapse="+")
    mm <- model.matrix(as.formula(f), data=pb.cd)
    
    # get cell-level formula
    formula <- paste(c("~(1|sample_id)+offset(ls)", covs, "group_id"), collapse = "+")
    if (verbose) print(formula)
    formula <- as.formula(paste("u", formula))
    
    # get coefficient to test
    if (is.null(coef)) {
        coef <- paste0("group_id", last(levels(x$group_id)))
        if (verbose) 
            message("Argument 'coef' not specified; ", 
                    sprintf("testing for %s.", dQuote(coef)))
    }
    
    # pseudo-bulk analysis:
    res <- pbDS(x, pb, mm, coef=which(colnames(mm)==coef), method="edgeR")
    res <- res$table[[1]][[1]]
    res <- res[,setdiff(colnames(res), c("F","p_adj.loc","coef","p_adj.glb"))] %>% 
        rename(pb.p_val="p_val")
    row.names(res) <- res$gene
    names(wg) <- wg <- as.character(res$gene)[which(res$pb.p_val < pthres4mm)]
    
    # fit mixed model for ea. gene which is below threshold in pseudobulk analysis
    fits <- bplapply(wg, function(i) {
        if(fam=="nbinom"){
            return(.fit_nbinom(data.frame(u=y[i, ], cd), formula, coef, prepForModeration=FALSE))
        }
        .fit_bglmer(data.frame(u=y[i, ], cd), formula, coef, prepForModeration=FALSE)
    }, BPPARAM = MulticoreParam(n_threads, progressbar=verbose))
    
    res$glmm.estimate <- res$glmm.p_val <- NA_real_
    res[wg, c("glmm.estimate", "glmm.p_val")] <- t(bind_rows(fits))[,c(1,4)]
    
    res$geomean.p_val <- res$mean.p_val <- res$p_val <- res$pb.p_val
    
    res[wg,"geomean.p_val"] <- 10^-rowMeans(-log10(as.matrix(res[wg,c("pb.p_val","glmm.p_val")])))
    res[wg,"mean.p_val"] <- rowMeans(as.matrix(res[wg,c("pb.p_val","glmm.p_val")]))
    res[wg,"p_val"] <- res[wg,"glmm.p_val"]
    
    res[,-1:-2] %>% add_column(.after = "p_val", p_adj.loc = p.adjust(.$p_val))
}


# fits mixed models and returns fit information required for eBayes
#' @importFrom blme blmer
#' @importFrom lmerTest lmer contest
#' @importFrom lme4 .makeCC lmerControl
#' @importFrom purrr map set_names
#' @importFrom stats residuals
#' @importFrom utils getFromNamespace
.fit_lmer <- function(df, formula, coef, bayesian, REML, ddf) {
    # here we should do some handling of convergence/singularity
    fun <- ifelse(bayesian, blmer, getFromNamespace("lmer", "lmerTest"))
    mod <- tryCatch(fun(formula, df, REML, control = lmerControl(
        check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))),
        error = function(e) e)
    if (inherits(mod, "error")) mod
    tryCatch(.prep_mmFit(mod, coef, ddf), error=function(e) NULL)
}

# fits negative binomial mixed models and returns fit information required for eBayes
#' @import glmmTMB
#' @importFrom purrr map set_names
#' @importFrom stats residuals
.fit_nbinom <- function(df, formula, coef, prepForModeration=FALSE){
    mod <- tryCatch({
        glmmTMB(formula, family=nbinom1, data=df)
    }, error=function(e){ print(e); return(NULL)})
    if(prepForModeration){
        return(tryCatch(.prep_mmFit(mod, coef), error=function(e) NULL))
    }
    tryCatch(coef(summary(mod))[[1]][coef,], error=function(e) rep(NA_real_, 4))
}

# fits poisson mixed models and returns fit information required for eBayes
#' @import blme
#' @importFrom purrr map set_names
#' @importFrom stats residuals
.fit_bglmer <- function(df, formula, coef, prepForModeration=FALSE){
    mod <- tryCatch({
        bglmer(formula, family="poisson", data=df)
    }, error=function(e){ print(e); return(NULL)})
    if(prepForModeration){
        return(tryCatch(.prep_mmFit(mod, coef), error=function(e) NULL))
    }
    tryCatch(coef(summary(mod))[coef,], error=function(e) rep(NA_real_, 4))
}

# prepares a model fit for eBayes moderation
.prep_mmFit <- function(mod, coef, ddf=NULL){
    if(is(mod, "glmmTMB")){
        co <- coef(mod)[[1]][[1]]
    }else{
        co <- coef(mod)[[1]]
    }
    beta <- as.matrix(co[1,coef])
    
    L <- matrix(as.numeric(colnames(co)==coef), ncol=1)
    row.names(L) <- colnames(co)
    
    # lmer method and SE modified from variancePartition:::.eval_lmm  
    if(is(mod, "lmerModLmerTest")){
        test <- lmerTest::contest(mod, t(L), ddf = ddf)
        ll <- list( df.residual=test[,"NumDF"], 
                    stat=as.matrix(test[,"F value"]), p_val0=test[,"Pr(>F)"])
    }else{
        test <- coef(summary(mod))
        if(is(mod, "glmmTMB")) test <- test[[1]]
        ll <- list( df.residual=df.residual(mod),
                    stat=as.matrix(test[coef,"z value"]), p_val0=test[coef,"Pr(>|z|)"])
    }
    if(is(mod, "lmerModLmerTest") && ddf=="Kenward-Roger") {
        V = pbkrtest::vcovAdj.lmerMod(mod, 0)
    }else{
        V = vcov(mod)
    }
    SE <- as.matrix(sapply(1:ncol(L), FUN=function(j){
        as.matrix(sqrt(sum(L[, j] * (V %*% L[, j]))), ncol = 1)
    }))
    colnames(SE) <- colnames(beta) <- "logFC"
    ll$coefficients <- beta
    ll$stdev.unscaled <- SE/sigma(mod)
    ll$Amean <- mean(co[,1])
    ll$sigma <- sigma(mod)
    ll
}


# formats a list of .prep_mmFit results into
# an eBayes compatible list & performs moderation
#' @importFrom dplyr %>% bind_cols pull
#' @importFrom limma eBayes
#' @importFrom magrittr set_colnames
#' @importFrom purrr map set_names
.mm_eBayes <- function(fits, coef, trended = FALSE) {
    rmv <- vapply(fits, inherits, what = "error", logical(1))
    f <- fits[!rmv]
    if (length(f) > 0) {
        res <- as.data.frame(t(sapply(f, FUN=function(x){ sapply(x, as.numeric) })))
        saveRDS(res, file="TMP.rds")
        res <- eBayes(res, trend = trended, robust = TRUE)
        res <- res[,c("coefficients", "stat", "p_val0", "t", "p.value")]
    } else {
        res <- matrix(NA, nrow = 0, ncol = 4) %>% 
            as.data.frame
    }
    if (any(rmv)) {
        res[names(which(rmv)), "error"] <- 
            vapply(fits[rmv], as.character, character(1))
    }
    res[names(fits), ] %>% set_colnames(c("beta", "stat", "p_val0", "t", "p_val"))
}

getPBcolData <- function(sce){
    CD <- colData(sce)
    nsamps <- length(unique(CD$sample_id))
    # get rid of most columns that would slow down agg:
    CD <- CD[,which(sapply(CD, FUN=function(x) length(unique(x)))<=nsamps),drop=FALSE]
    pb <- aggregate(CD, by=list(CD$sample_id), FUN=function(x){
        ifelse(length(unique(x))>1, NA, x[1])
    })
    # re-factor:
    for(f in intersect(colnames(pb), colnames(CD))){
        if(is.factor(CD[[f]])) pb[[f]] <- factor(levels(CD[[f]])[pb[[f]]], levels=levels(CD[[f]]))
    }
    pb <- pb[,colSums(is.na(as.matrix(pb)))<nrow(pb),drop=FALSE]
    # save the number of cells:
    pb$nbCells <- as.numeric(table(CD$sample_id)[pb$sample_id])
    row.names(pb) <- pb[,1]
    pb <- pb[,-1]
    pb$sample_id <- NULL
    pb
}

