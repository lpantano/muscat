% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{data}
\alias{data}
\alias{sce}
\title{Example datasets}
\value{
a \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
}
\description{
A \code{\link[SingleCellExperiment]{SingleCellExperiment}} containing 
10x droplet-based scRNA-seq PBCM data from 8 Lupus patients befor and after 
6h-treatment with INF-beta (16 samples in total).

The original data has been minimally filtered to remove undetected genes, 
cell multiplets, and unassigned cells. 
Assay \code{logcounts} corresponds to log-normalized values 
obtained from \code{\link[scater]{normalize}} with default parameters.
  
The complete raw, gene, and cell metadata is available through the NCBI GEO,
accession number GSE96583.
}
\references{
Kang et al. (2018). 
Multiplexed droplet single-cell RNA-sequencing 
using natural genetic variation. 
\emph{Nature Biotechnology}, 
\bold{36}(1): 89-94. DOI: 10.1038/nbt.4042.
}
\author{
Helena L. Crowell \email{helena.crowell@uzh.ch}.
}
