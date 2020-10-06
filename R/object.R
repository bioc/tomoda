#' Create an object from matrix
#'
#' \code{createTomo.matrix} creates an object from raw read count matrix or normalized read count matrix.
#'
#' @param matrix.count A numeric matrix or matrix-like data stucture that can be coverted to matrix, with genes with rows, sections as columns and values as raw read counts. Columns should be sorted according to section numbers.
#' @param matrix.normalized A numeric matrix or matrix-like data stucture that can be coverted to matrix, with genes as rows, sections as columns and values as normalized read counts. Columns should be sorted according to order of sections.
#' @param min.section Integer. Genes expressed in less than \code{min.section} sections will be filtered out.
#' @param normalize Logical, whether to perform normalization when creating the object. Default is TRUE.
#' @param normalize.method Character, must be one of \code{"median"}, or \code{"cpm"}.
#' @param scale Logical, whether to perform scaling when creating the object. Default is TRUE.
#'
#' @return A \code{SummarizedExperiment} object
#'
#' @seealso \code{\link{createTomo}} for the generic function.
#'
#' @importFrom methods new
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo.matrix(zh.data)
#'
createTomo.matrix <- function(matrix.count=NULL, matrix.normalized=NULL, min.section=3, normalize=TRUE, normalize.method='median', scale=TRUE)
{
    available.count <- !is.null(matrix.count)
    available.normalized <- !is.null(matrix.normalized)

    if(!available.count & !available.normalized)
    {
        message('Please input at least one of raw read count matrix and normalized count matrix!\n')
        return()
    }
    else if(available.count & available.normalized)
    {
        good_gene <- apply(matrix.normalized>0, 1, sum) >= min.section
        matrix.count <- as.matrix(matrix.count[good_gene,,drop=FALSE])
        matrix.normalized <- as.matrix(matrix.normalized[good_gene,,drop=FALSE])
        section <- colnames(matrix.count)
        if(is.null(section))
            section <- as.character(seq_len(ncol(matrix.count)))

        col_data <- data.frame(section=section, stringsAsFactors=FALSE)
        row_data <- data.frame(gene=rownames(matrix.count), stringsAsFactors=FALSE)
        tomo_object <- SummarizedExperiment(assays=list(
            count=matrix.count,
            normalized=matrix.normalized),
            colData=col_data,
            rowData=row_data)
    }
    else if(available.count)
    {
        good_gene <- apply(matrix.count>0, 1, sum) >= min.section
        matrix.count <- as.matrix(matrix.count[good_gene,,drop=FALSE])
        section <- colnames(matrix.count)
        if(is.null(section))
            section <- as.character(seq_len(ncol(matrix.count)))

        col_data <- data.frame(section=section, stringsAsFactors=FALSE)
        row_data <- data.frame(gene=rownames(matrix.count), stringsAsFactors=FALSE)
        tomo_object <- SummarizedExperiment(assays=list(
            count=matrix.count),
            colData=col_data,
            rowData=row_data)
    }
    else
    {
        good_gene <- apply(matrix.normalized>0, 1, sum) >= min.section
        matrix.normalized <- as.matrix(matrix.normalized[good_gene,,drop=FALSE])
        section <- colnames(matrix.normalized)
        if(is.null(section))
            section <- as.character(seq_len(ncol(matrix.normalized)))

        col_data <- data.frame(section=section, stringsAsFactors=FALSE)
        row_data <- data.frame(gene=rownames(matrix.normalized), stringsAsFactors=FALSE)
        tomo_object <- SummarizedExperiment(assays=list(
            normalized=matrix.normalized),
            colData=col_data,
            rowData=row_data)
    }

    if(normalize)
        tomo_object <- normalizeTomo(tomo_object, normalize.method)
    if(scale)
        tomo_object <- scaleTomo(tomo_object)

    return(tomo_object)
}

#' Create an object from SummarizedExperiment
#'
#' \code{createTomo.SummarizedExperiment} creates an object from a SummarizedExperiment object.
#'
#' @param se A SummarizedExperiment object, it must contain at least one of 'count' assay and 'normalized' assay.
#' @param min.section Integer. Genes expressed in less than \code{min.section} sections will be filtered out.
#' @param normalize Logical, whether to perform normalization when creating the object. Default is TRUE.
#' @param normalize.method Character, must be one of \code{"median"}, or \code{"cpm"}.
#' @param scale Logical, whether to perform scaling when creating the object. Default is TRUE.
#'
#' @return A \code{SummarizedExperiment} object
#'
#' @seealso \code{\link{createTomo}} for the generic function.
#'
#' @importFrom methods new
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
#'
#' @examples
#' data(zh.data)
#' se <- SummarizedExperiment::SummarizedExperiment(assays=list(count=zh.data))
#' zh <- createTomo.SummarizedExperiment(se)
createTomo.SummarizedExperiment <- function(se, min.section=3, normalize=TRUE, normalize.method='median', scale=TRUE)
{
    available.count <- 'count' %in% assayNames(se)
    available.normalized <- 'normalized' %in% assayNames(se)

    if(!available.count & !available.normalized)
    {
        message("The input SummarizedExperiment object must have at least one of 'count' assay and 'normalized' assay!\n")
        return()
    }
    else if(available.count & available.normalized)
    {
        good_gene <- apply(assay(se, 'normalized')>0, 1, sum) >= min.section
        matrix.count <- assay(se, 'count')[good_gene,,drop=FALSE]
        matrix.normalized <- assay(se, 'normalized')[good_gene,,drop=FALSE]
        section <- colnames(matrix.count)
        if(is.null(section))
            section <- as.character(seq_len(ncol(matrix.count)))

        col_data <- data.frame(section=section, colData(se), stringsAsFactors=FALSE)
        row_data <- data.frame(gene=rownames(matrix.count),
                               rowData(se)[good_gene,,drop=FALSE], stringsAsFactors=FALSE)
        tomo_object <- SummarizedExperiment(assays=list(
            count=matrix.count,
            normalized=matrix.normalized),
            colData=col_data,
            rowData=row_data)
    }
    else if(available.count)
    {
        good_gene <- apply(assay(se, 'count')>0, 1, sum) >= min.section
        matrix.count <- assay(se, 'count')[good_gene,,drop=FALSE]
        section <- colnames(matrix.count)
        if(is.null(section))
            section <- as.character(seq_len(ncol(matrix.count)))

        col_data <- data.frame(section=section, colData(se), stringsAsFactors=FALSE)
        row_data <- data.frame(gene=rownames(matrix.count),
                               rowData(se)[good_gene,,drop=FALSE], stringsAsFactors=FALSE)
        tomo_object <- SummarizedExperiment(assays=list(
            count=matrix.count),
            colData=col_data,
            rowData=row_data)
    }
    else
    {
        good_gene <- apply(assay(se, 'normalized')>0, 1, sum) >= min.section
        matrix.normalized <- assay(se, 'normalized')[good_gene,,drop=FALSE]
        section <- colnames(matrix.normalized)
        if(is.null(section))
            section <- as.character(seq_len(ncol(matrix.normalized)))

        col_data <- data.frame(section=section,
                               colData(se), stringsAsFactors=FALSE)
        row_data <- data.frame(gene=rownames(matrix.normalized),
                               rowData(se)[good_gene,,drop=FALSE], stringsAsFactors=FALSE)
        tomo_object <- SummarizedExperiment(assays=list(
            normalized=matrix.normalized),
            colData=col_data,
            rowData=row_data)
    }

    if(normalize)
        tomo_object <- normalizeTomo(tomo_object, normalize.method)
    if(scale)
        tomo_object <- scaleTomo(tomo_object)

    return(tomo_object)
}

#' @rdname createTomo
#' @aliases createTomo,SummarizedExperiment-method
setMethod('createTomo', signature(object='SummarizedExperiment'),
          function(object, min.section=3, normalize=TRUE, normalize.method='median', scale=TRUE)
              createTomo.SummarizedExperiment(se=object,
                                              min.section=min.section,
                                              normalize=normalize,
                                              normalize.method=normalize.method,
                                              scale=scale)
)

#' @rdname createTomo
#' @aliases createTomo,matrix-method
setMethod('createTomo', signature(object='matrix'),
          function(object, matrix.normalized=NULL, min.section=3, normalize=TRUE, normalize.method='median', scale=TRUE)
              createTomo.matrix(matrix.count=object,
                                matrix.normalized=matrix.normalized,
                                min.section=min.section,
                                normalize=normalize,
                                normalize.method=normalize.method,
                                scale=scale)
)

#' @rdname createTomo
#' @aliases createTomo,missing-method
setMethod('createTomo', signature(object='missing'),
          function(matrix.normalized=NULL, min.section=3, normalize=TRUE, normalize.method='median', scale=TRUE, ...)
              createTomo.matrix(matrix.normalized=matrix.normalized,
                                min.section=min.section,
                                normalize=normalize,
                                normalize.method=normalize.method,
                                scale=scale)
)
