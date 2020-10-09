#' Create an object representing tomo-seq data
#'
#' This is a generic function to create an object representing tomo-seq data. The input object can either be a \code{matrix} or a \code{SummarizeExperiment}.
#'
#' @param object Either a raw read count matrix or a SummarizedExperiment object.
#' @param matrix.normalized (Optional) A numeric matrix of normalized read count.
#' @param min.section Integer. Genes expressed in less than \code{min.section} sections will be filtered out.
#' @param normalize Logical, whether to perform normalization when creating the object. Default is TRUE.
#' @param normalize.method Character, must be one of \code{"median"}, or \code{"cpm"}.
#' @param scale Logical, whether to perform scaling when creating the object. Default is TRUE.
#' @param ... Additional parameters to pass to S4 methods.
#'
#' @details This is the generic function to create a \code{SummarizedExperiment} object for representing tomo-seq data. Either \code{matrix} or \code{SummarizedExperiment} object can be used for input.
#'
#' When using \code{matrix} for input, at least one of raw read count matrix and normalized read count matrix (like FPKM and TPM) must be used for input. If normalized matrix is available, input it with argument \code{matrix.normalized}. Matrices should have genes as rows and sections as columns. Columns should be sorted according to the order of sections.
#'
#' When using \code{SummarizedExperiment} object for input, it must contain at least one of 'count' assay and 'normalized' assay. Besides, the row data and column data of the input object will be retained in the output object.
#'
#' By default, all library sizes are normalized to the median library size across sections. Set \code{normalize.method = "cpm"} will make library sizes normalized to 1 million counts.
#' Scaling and centering is performed for all genes across sections.
#'
#' @return A \code{SummarizedExperiment} object. Raw read count matrix, normalized read count matrix and scaled read count matrix are saved in 'count', 'normalized' and 'scale' assays of the object.
#'
#' @seealso \itemize{
#'  \item{\code{\link{createTomo.matrix}}} : creating an object from \code{matrix}.
#'  \item{\code{\link{createTomo.SummarizedExperiment}}} : creating an object from \code{SummarizedExperiment}.
#'  \item{\code{\link{normalizeTomo}}} : normalization.
#'  \item{\code{\link{scaleTomo}}} : scaling.
#'  \item{\code{\link[SummarizedExperiment]{SummarizedExperiment-class}}} : operations on \code{SummarizedExperiment}.
#' }
#'
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo(zh.data)
#'
#' data(zh.data)
#' se <- SummarizedExperiment::SummarizedExperiment(assays=list(count=zh.data))
#' zh <- createTomo(se)
setGeneric('createTomo', function(object, ...) standardGeneric('createTomo'))
