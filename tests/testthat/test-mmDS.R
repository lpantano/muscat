context("DS analysis using mixed-models")

# load packages
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(SingleCellExperiment)
})

# generate toy dataset
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)
sce <- toySCE()
cs <- sce$cluster_id %in% sample(sce$cluster_id, 2)
sce <- sce[sample(nrow(sce), 50), cs]
sce$cluster_id <- droplevels(sce$cluster_id)

kids <- sce$cluster_id
sids <- sce$sample_id
gids <- sce$group_id

# randomly select n_de DE genes & multiply counts by 10 for group 2
n_de <- 5
de_gs <- sample(rownames(sce), n_de)
g3 <- gids == "g3"
assay(sce[de_gs, g3]) <- assay(sce[de_gs, g3]) + 100

res <- mmDS(sce, method = "poisson", verbose = FALSE)

tbl <- map(res, dplyr::filter, p_adj.loc < 1e-12)
de_gs_res <- map(tbl, "gene")
n_de_res <- vapply(de_gs_res, length, numeric(1))
n_de_res
expect_true(all(n_de_res == n_de))


for (vst in c("DESeq2", "sctransform")) {
    test_that(paste0("mmDS-vst-", vst), {
        res <- mmDS(sce, method = "vst", vst = vst, verbose = TRUE)
        
        expect_is(res, "list")
        expect_identical(names(res), levels(kids))
        
        tbl <- map(res, dplyr::filter, p_adj.loc < 1e-3)
        de_gs_res <- map(tbl, "gene")
        n_de_res <- vapply(de_gs_res, length, numeric(1))
        expect_true(all(n_de_res == n_de))
        
        os <- map(p_adj, order)
        de_gs_res <- map(os, function(o) rownames(sce)[o][seq_len(n_de)])
        expect_true(all(unlist(map(de_gs_res, setequal, de_gs))))
    })
}



