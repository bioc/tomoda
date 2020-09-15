#' Create a \code{SummarizedExperiment} object
#'
#' Create a \code{SummarizedExperiment} object from raw read counts or normalized read counts.
#'
#' @param matrix.count A numeric matrix or matrix-like data stucture that can be coverted to matrix, with genes with rows, sections as columns and values as raw read counts.
#' Columns should be sorted according to section numbers.
#' @param matrix.normalized A numeric matrix or matrix-like data stucture that can be coverted to matrix, with genes with rows, sections as columns and values as normalized read counts.
#' Columns should be sorted according to section numbers.
#' @param min.section An integer. Genes expressed in less than \code{min.section} sections will be filtered out.
#'
#' @details This is the method to create a \code{Tomo} object, which is an instance of class \code{SummarizedExperiment}. At least one of raw read count matrix and normalized raw count matrix (like FPKM, TPM, ...) must be used for input.
#'
#' @return A \code{SummarizedExperiment} object
#'
#' @importFrom methods new
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
CreateTomo <- function(matrix.count=NULL, matrix.normalized=NULL, min.section=3)
{
    available.count <- !is.null(matrix.count)
    available.normalized <- !is.null(matrix.normalized)

    if(!available.count & !available.normalized)
    {
        print('Please input at least one of raw read count matrix and normalized count matrix!\n')
        return()
    }
    else if(available.count & available.normalized)
    {
        good_gene <- apply(matrix.normalized>0, 1, sum) >= min.section
        matrix.count <- as.matrix(matrix.count[good_gene,])
        matrix.normalized <- as.matrix(matrix.normalized[good_gene,])
        section <- colnames(matrix.count)
        if(is.null(section))
            section <- as.character(seq_len(ncol(matrix.count)))

        col_data <- data.frame(section=section, stringsAsFactors=FALSE)
        row_data <- data.frame(gene=rownames(matrix.count), stringsAsFactors=FALSE)
        tomo_object <- SummarizedExperiment(assays=list(count=matrix.count, normalized=matrix.normalized), colData=col_data, rowData=row_data)
    }
    else if(available.count)
    {
        good_gene <- apply(matrix.count>0, 1, sum) >= min.section
        matrix.count <- as.matrix(matrix.count[good_gene,])
        section <- colnames(matrix.count)
        if(is.null(section))
            section <- as.character(seq_len(ncol(matrix.count)))
        col_data <- data.frame(section=section, stringsAsFactors=FALSE)
        row_data <- data.frame(gene=rownames(matrix.count), stringsAsFactors=FALSE)
        tomo_object <- SummarizedExperiment(assays=list(count=matrix.count), colData=col_data, rowData=row_data)
    }
    else
    {
        good_gene <- apply(matrix.normalized>0, 1, sum) >= min.section
        matrix.normalized <- as.matrix(matrix.normalized[good_gene,])
        section <- colnames(matrix.count)
        if(is.null(section))
            section <- as.character(seq_len(ncol(matrix.normalized)))
        col_data <- data.frame(section=section, stringsAsFactors=FALSE)
        row_data <- data.frame(gene=rownames(matrix.normalized), stringsAsFactors=FALSE)
        tomo_object <- SummarizedExperiment(assays=list(normalized=matrix.normalized), colData=col_data, rowData=row_data)
    }

    return(tomo_object)
}

setGeneric("Normalize",function(object, ...) standardGeneric("Normalize"))
setMethod("Normalize", signature(object="SummarizedExperiment"), Normalize)

setGeneric("Scale",function(object, ...) standardGeneric("Scale"))
setMethod("Scale", signature(object="SummarizedExperiment"), Scale)

setGeneric("HClust",function(object, ...) standardGeneric("HClust"))
setMethod("HClust", signature(object="SummarizedExperiment"), HClust)

setGeneric("KMeans",function(object, ...) standardGeneric("KMeans"))
setMethod("KMeans", signature(object="SummarizedExperiment"), KMeans)

setGeneric("FindPeakGenes",function(object, ...) standardGeneric("FindPeakGenes"))
setMethod("FindPeakGenes", signature(object="SummarizedExperiment"), FindPeakGenes)

setGeneric("PCA",function(object, ...) standardGeneric("PCA"))
setMethod("PCA", signature(object="SummarizedExperiment"), PCA)

setGeneric("TSNE",function(object, ...) standardGeneric("TSNE"))
setMethod("TSNE", signature(object="SummarizedExperiment"), TSNE)

setGeneric("UMAP",function(object, ...) standardGeneric("UMAP"))
setMethod("UMAP", signature(object="SummarizedExperiment"), UMAP)

setGeneric("LinePlot",function(object, ...) standardGeneric("LinePlot"))
setMethod("LinePlot", signature(object="SummarizedExperiment"), LinePlot)

setGeneric("EmbedPlot",function(object, ...) standardGeneric("EmbedPlot"))
setMethod("EmbedPlot", signature(object="SummarizedExperiment"), EmbedPlot)

setGeneric("GeneEmbedPlot",function(object, ...) standardGeneric("GeneEmbedPlot"))
setMethod("GeneEmbedPlot", signature(object="SummarizedExperiment"), GeneEmbedPlot)

setGeneric("ExpHeatmap",function(object, ...) standardGeneric("ExpHeatmap"))
setMethod("ExpHeatmap", signature(object="SummarizedExperiment"), ExpHeatmap)

setGeneric("CorHeatmap",function(object, ...) standardGeneric("CorHeatmap"))
setMethod("CorHeatmap", signature(object="SummarizedExperiment"), CorHeatmap)

setGeneric("GeneCorHeatmap",function(object, ...) standardGeneric("GeneCorHeatmap"))
setMethod("GeneCorHeatmap", signature(object="SummarizedExperiment"), GeneCorHeatmap)
