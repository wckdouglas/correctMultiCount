# This file was generated by Rcpp::compileAttributes
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @useDynLib correctMultiCount
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' @title correctCounts
#' @description Assigning multiply-mapped reads to most-abundantly-mapped locus
#' @param basecount: should be a dataframe for storing gene count from uniquely-mapped reads
#'                   **Columns named as 'id' for gene_id and 'count' for counts
#' @param multiCount should be a dataframe for storing multiply-mapped reads and their mapped locus.
#'                 **Columns named as 'fragment_id' for storing read id  and 'gene_id' for storing gene_id
#' @export
correctCounts <- function(baseCount, multiCount) {
    .Call('correctMultiCount_correctCounts', PACKAGE = 'correctMultiCount', baseCount, multiCount)
}

