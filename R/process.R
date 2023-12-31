#' Normalize data
#'
#' Normalize the raw read count in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param method Character, must be one of \code{"median"}, or \code{"cpm"}.
#'
#' @details This function should be run for \code{SummarizedExperiment} object created from raw read count matrix.
#' If the \code{SummarizedExperiment} object already has a normalized count matrix. The function simply return the original object.
#' Library sizes of all sections are normalized to the median library size (method='median') or one million (method='cpm').
#'
#' @return A \code{SummarizedExperiment} object with normalized read count matrix saved in assay \code{'normalized'}.
#'
#' @importFrom stats median
#' @importFrom SummarizedExperiment assayNames assay assay<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo(zh.data, normalize=FALSE)
#' zh <- normalizeTomo(zh)
normalizeTomo <- function(object, method='median')
{
    if("normalized" %in% assayNames(object))
        message('Normalized data already exist!')
    else
    {
        library_size <- apply(assay(object, 'count'), 2, sum)
        target_library_size <- ifelse(method=='cpm', 1e6, median(library_size))
        if(!method %in% c('median', 'cpm'))
            message('Unknown normalization method. Using median library size.')
        assay(object, 'normalized') <- apply(assay(object, 'count'), 2,
                                           function(x) x * target_library_size / sum(x))

        message("Normalized count matrix is saved in assay 'normalized'.\n")
    }
    return(object)
}

#' Scale data
#'
#' Scale the normalized read count in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#'
#' @details This function should be run for \code{SummarizedExperiment} object with normalized read count matrix.
#' The normalized read counts of each gene are subjected to Z score transformation across sections.
#'
#' @return A \code{SummarizedExperiment} object with scaled read count matrix saved in assay \code{'scaled'}.
#'
#' @importFrom stats sd
#' @importFrom SummarizedExperiment assayNames assay assay<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo(zh.data, scale=FALSE)
#' zh <- scaleTomo(zh)
scaleTomo <- function(object)
{
    if(!"normalized" %in% assayNames(object))
        message("Normalized data does not exist! Please run function 'normalizeTomo' before scaling data.\n")
    else
    {
        mean_normalized <- apply(assay(object, 'normalized'), 1, mean)
        sd_normalized <- apply(assay(object, 'normalized'), 1, sd)
        assay(object, 'scaled') <- (assay(object, 'normalized') - mean_normalized) / sd_normalized
        message("Scaled count matrix is saved in assay 'scaled'.\n")
    }
    return(object)
}

#' Hierarchical clustering across sections
#'
#' Performs hierarchical clustering across sections in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param matrix Character, must be one of \code{"count"}, \code{"normalized"}, or \code{"scaled"}.
#' @param measure Character, must be one of \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"} or \code{"minkowski"}.
#' @param p Numeric, the power of the Minkowski distance.
#' @param agglomeration Character, must be one of \code{"ward.D"}, \code{"ward.D2"}, \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"median"} or \code{"centroid"}.
#'
#' @return A \code{hclust} object.
#'
#' @seealso \code{\link[stats]{dist}} for measuring distance and \code{\link[stats]{hclust}} for performing hierarchical clustering on a matrix.
#'
#' @importFrom stats dist hclust
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo(zh.data)
#' hclust_zh <- hierarchClust(zh)
#' plot(hclust_zh)
#'
#' # Use other agglomeration method
#' hclust_zh <- hierarchClust(zh, agglomeration="average")
#'
#' # (Not recommended) Use scaled read counts to calculate distance
#' zh <- scaleTomo(zh)
#' hclust_zh <- hierarchClust(zh, matrix="scaled")
hierarchClust <- function(object, matrix='normalized', measure='euclidean', p=2, agglomeration='complete')
{
    exp_matrix <- t(assay(object, matrix))
    dist_section <- dist(exp_matrix, method=measure, p)
    hclust_section <- hclust(dist_section, method=agglomeration)
    return(hclust_section)
}

#' K-Means clustering across sections
#'
#' Performs K-Means clustering across sections in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param centers Integer, number of clusters, namely \emph{k}.
#' @param matrix Character, must be one of \code{"count"}, \code{"normalized"}, or \code{"scaled"}.
#' @param ... other parameters passed to \code{kmeans}.
#'
#' @return A \code{SummarizedExperiment} object. The obtained cluster labels are saved in slot \code{meta}.
#'
#' @seealso \code{\link[stats]{kmeans}} for performing K-Means clustering on a matrix.
#'
#' @importFrom stats kmeans
#' @importFrom SummarizedExperiment assay colData colData<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo(zh.data)
#' zh <- kmeansClust(zh, 3)
#'
#' # Use scaled read counts to calculate distance
#' zh <- scaleTomo(zh)
#' zh <- kmeansClust(zh, 3, matrix="scaled")
kmeansClust <- function(object, centers, matrix='normalized', ...)
{
    exp_matrix <- t(assay(object, matrix))
    kmeans_section <- kmeans(exp_matrix,centers=centers, ...)
    percent_between <- kmeans_section$betweenss / kmeans_section$totss
    colData(object)$kmeans_cluster <- kmeans_section$cluster

    cluster_list <- list()
    for(i in seq_len(centers))
        cluster_list[[i]] <- as.character(colData(object)$section)[kmeans_section$cluster == i]

    message("K-Means clustering labels are saved in colData.")
    message("between_SS / total_SS =", percent_between,'\n')
    return(object)
}

#' Find peak in a vector
#'
#' Find the position of peak in a vector.
#'
#' @param x A numeric vector.
#' @param threshold Integer, only values bigger than \code{threshold} are recognized as part of the peak.
#' @param length Integer, minimum \code{length} of consecutive values bigger than \code{threshold} are recognized as a peak.
#'
#' @return A numeric vector. The first element is the start index and the second element is the end index of the peak.
#' If multiple peaks exist, only output the start and end index of the one with maximun length.
#' If no peak exist, return \code{c(0, 0)}.
#'
#' @export
#'
#' @examples
#' # return c(3, 10)
#' findPeak(c(0:5, 5:0), threshold=1, length=4)
#'
#' # Most likely return c(0, 0)
#' findPeak(rnorm(10), threshold=3, length=3)
findPeak <- function(x, threshold=1, length=4)
{
    gt_threshold <- x > threshold
    consecutive <- rep(0, length(x) + 1)

    for(i in seq_len(length(x)))
    {
        if(gt_threshold[i])
        {
            consecutive[i+1] <- consecutive[i] + 1
        }
    }
    if(all(consecutive < length))
    {
        return(c(0, 0))
    }
    else
    {
        peak_end <- which.max(consecutive)
        peak_start <- peak_end - consecutive[peak_end]
        return(c(peak_start, peak_end - 1) )
    }
}

#' Find peak genes
#'
#' Find peak genes (spatially upregulated genes) in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param threshold Integer, only scaled read counts bigger than \code{threshold} are recognized as part of the peak.
#' @param length Integer, scaled read counts bigger than \code{threshold} in minimum \code{length} of consecutive sections are recognized as a peak.
#' @param matrix Character, must be one of \code{"count"}, \code{"normalized"}, or \code{"scaled"}.
#' @param nperm Integer, number of random permutations to calculate p values. Set it to 0 if you are not interested in p values.
#' @param method Character, the method to adjust p values for multiple comparisons, must be one of \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}.
#'
#' @details Peak genes are selected based on scaled read counts. As scaled read counts are Z scores, suggested \code{threshold} are \eqn{[1,3]}.
#' Smaller \code{threshold} and \code{length} makes the function detect more peak genes, and vice versa.
#' P values are calculated by approximate permutation tests. For a given \code{threshold} and \code{length}, the scaled read counts of each gene is randomly permutated for \code{nperm} times. The p value is defined as the ratio of permutations containing peaks.
#' In order to speed up permutation process, genes whose expression exceeds threshold in same number of sections have same p values.
#' To be specific, only one of these genes will be used to calculate a p value by permutation, and other genes are assigned this p value.
#'
#' @return A data.frame with peak genes as rows. It has following columns:
#' \itemize{
#'  \item{\code{gene}} : Character, peak gene names.
#'  \item{\code{start}} : Numeric, the start index of peak.
#'  \item{\code{end}} : Numeric, the end index of peak.
#'  \item{\code{center}} : Numeric, the middle index of peak. If the length of the peak is even, \code{center} is defined as the left-middle index.
#'  \item{\code{p}} : Numeric, p values.
#'  \item{\code{p.adj}} : Numeric, adjusted p values.
#' }
#'
#' @importFrom SummarizedExperiment assay assayNames rowData
#' @importFrom stats p.adjust
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo(zh.data)
#' peak_genes <- findPeakGene(zh)
#' head(peak_genes)
#'
#' # Increase threshold so that less peak genes will be found.
#' peak_genes <- findPeakGene(zh, threshold=1.5)
#'
#' # Increase peak length so that less peak genes will be found.
#' peak_genes <- findPeakGene(zh, length=5)
#'
#' # Set nperm to 0 so that p values will not be calculated. This will save running time.
#' peak_genes <- findPeakGene(zh, nperm=0)
findPeakGene <- function(object, threshold=1, length=4, matrix='scaled', nperm=1e5, method='BH')
{
    if(matrix=='scaled' & !"scaled" %in% assayNames(object))
    {
        message("Scaled data does not exist! Please run function 'scaleTomo' before finding peak genes.\n")
        return()
    }
    else
    {
        exp_matrix <- assay(object, matrix)
        peak_position <- apply(exp_matrix, 1, findPeak, threshold, length)
        peak_exist <- peak_position[1,] != 0
        if(!any(peak_exist))
        {
            message("No peak gene is found!\n")
            return()
        }

        peak_genes <- rowData(object)$gene[peak_exist]
        peak_gene_df <- data.frame(gene=peak_genes,
                                   start=peak_position[1,peak_exist],
                                   end=peak_position[2,peak_exist],
                                   center=floor(apply(peak_position[,peak_exist], 2, mean)),
                                   stringsAsFactors=FALSE)

        # Using permutation to calculate p-values
        n_peak_genes <- length(peak_genes)
        if(nperm > 0)
        {
            pvals <- rep(NA, n_peak_genes)
            n_section <- nrow(object)
            saved_pvals <- rep(NA, n_section)

            for(i in seq_len(n_peak_genes))
            {
                exp_gene <- exp_matrix[peak_genes[i], ]
                n_gt_threshold <- sum(exp_gene > threshold)
                if(is.na(saved_pvals[n_gt_threshold]))
                {
                    n_peak <- 0
                    for(perm in seq_len(nperm))
                        n_peak <- n_peak + (findPeak(sample(exp_gene), threshold, length)[1] > 0)
                    pval <- n_peak / nperm
                    saved_pvals[n_gt_threshold] <- pval
                }
                else
                {
                    pval <- saved_pvals[n_gt_threshold]
                }
                pvals[i] <- pval
            }

            peak_gene_df$p=pvals
            peak_gene_df$p.adj=p.adjust(pvals, method=method)
        }

        sorted_df <- peak_gene_df[order(peak_gene_df$start, peak_gene_df$end), ]
        message(nrow(sorted_df), "peak genes (spatially upregulated genes) are found!\n")
        return(sorted_df)
    }
}

#' Perform PCA
#'
#' Perform PCA on sections or genes in a \code{SummarizedExperiment} object for dimensionality reduction.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param genes \code{NA} or a vector of character. Perform PCA on sections if it is \code{NA}, or on given genes if it is a vector of gene names.
#' @param matrix Character, must be one of \code{"auto"}, \code{"count"}, \code{"normalized"}, or \code{"scaled"}. If \code{"auto"}, normalized matrix is used for sections and scaled matrix is used for genes.
#' @param scree Logical, plot the scree plot for PCs if it is \code{TRUE}.
#' @param ... Other parameters passed to \code{prcomp}.
#'
#' @return A \code{SummarizedExperiment} object. The PC embeddings are saved in slot \code{meta} if PCA is performed on sections, or saved in slot \code{gene_embedding} if PCA is performed on genes.
#'
#' @seealso \code{\link[stats]{prcomp}} for performing PCA on a matrix.
#'
#' @importFrom stats prcomp
#' @importFrom SummarizedExperiment assay rowData<- colData<-
#' @importFrom ggplot2 ggplot aes_string geom_point labs theme_bw
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo(zh.data)
#'
#' # Perform PCA on sections.
#' zh <- runPCA(zh)
#'
#' # Plot the scree plot.
#' zh <- runPCA(zh, scree=TRUE)
#'
#' # Perform PCA on some genes.
#' zh <- runPCA(zh, genes=rownames(zh)[1:100])
runPCA <- function(object, genes=NA, matrix='auto', scree=FALSE, ...)
{
    if(all(is.na(genes)))
    {
        if(matrix == 'auto')
            matrix <- 'normalized'
        pca_result <- prcomp(assay(object, matrix), ...)
        colData(object)$PC1 <- pca_result$rotation[,1]
        colData(object)$PC2 <- pca_result$rotation[,2]
        message("PC embeddings for sections are saved in column data.\n")
    }
    else
    {
        if(matrix == 'auto')
            matrix <- 'scaled'
        exp_matrix <- assay(object, matrix)[genes, ]
        pca_result <- prcomp(t(exp_matrix))
        pc1 <- pca_result$rotation[,1]
        pc2 <- pca_result$rotation[,2]
        rowData(object)$PC1 <- NA
        rowData(object)$PC2 <- NA
        rowData(object)[genes, 'PC1'] <- pc1
        rowData(object)[genes, 'PC2'] <- pc2
        message("PC embeddings for genes are saved in row data.\n")
    }

    if(scree)
    {
        pca_sd <- data.frame(pc=seq_len(length(pca_result$sdev)), sd=pca_result$sdev)
        g <- ggplot(pca_sd, aes_string(x='pc', y='sd')) +
            geom_point() +
            labs(x='PC', y='Standard deviation') +
            theme_bw()
        print(g)
    }
    return(object)
}

#' Perform TSNE
#'
#' Perform TSNE on sections or genes in a \code{SummarizedExperiment} object for dimensionality reduction.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param genes \code{NA} or a vector of character. Perform TSNE on sections if it is \code{NA}, or on given genes if it is a vector of gene names.
#' @param matrix Character, must be one of \code{"auto"}, \code{"count"}, \code{"normalized"}, or \code{"scaled"}. If \code{"auto"}, normalized matrix is used for sections and scaled matrix is used for genes.
#' @param perplexity Numeric, perplexity parameter for Rtsne (default: 0.25 *(number of observations - 1)).
#' @param ... Other parameters passed to \code{Rtsne}.
#'
#' @return A \code{SummarizedExperiment} object. The TSNE embeddings are saved in slot \code{meta} if TSNE is performed on sections, or saved in slot \code{gene_embedding} if TSNE is performed on genes.
#'
#' @seealso \code{\link[Rtsne]{Rtsne}} for performing TSNE on a matrix.
#'
#' @importFrom Rtsne Rtsne
#' @importFrom SummarizedExperiment assay rowData<- colData<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo(zh.data)
#'
#' # Perform TSNE on sections.
#' zh <- runTSNE(zh)
#'
#' # Perform TSNE on sections with other perplexity.
#' zh <- runTSNE(zh, perplexity=10)
#'
#' # Perform TSNE on some genes.
#' zh <- runTSNE(zh, genes=rownames(zh)[1:100])
runTSNE <- function(object, genes=NA, matrix='auto', perplexity=NA, ...)
{
    if(all(is.na(genes)))
    {
        if(matrix == 'auto')
            matrix <- 'normalized'
        if(is.na(perplexity))
            perplexity <- (ncol(object) - 1) / 4
        tsne_result <- Rtsne(t(assay(object, matrix)), perplexity=perplexity, ...)
        colData(object)$TSNE1 <- tsne_result$Y[,1]
        colData(object)$TSNE2 <- tsne_result$Y[,2]
        message("TSNE embeddings for sections are saved in column data.\n")
    }
    else
    {
        if(matrix == 'auto')
            matrix <- 'scaled'
        if(is.na(perplexity))
            perplexity <- (length(genes) - 1) / 4
        exp_matrix <- assay(object, matrix)[genes, ]
        tsne_result <- Rtsne(exp_matrix, perplexity=perplexity, ...)
        tsne1 <- tsne_result$Y[,1]
        tsne2 <- tsne_result$Y[,2]
        rowData(object)$TSNE1 <- NA
        rowData(object)$TSNE2 <- NA
        rowData(object)[genes, 'TSNE1'] <- tsne1
        rowData(object)[genes, 'TSNE2'] <- tsne2
        message("TSNE embeddings for genes are saved in row data.\n")
    }
    return(object)
}

#' Perform UMAP
#'
#' Perform UMAP on sections or genes in a \code{SummarizedExperiment} object for dimensionality reduction.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param genes \code{NA} or a vector of character. Perform UMAP on sections if it is \code{NA}, or on given genes if it is a vector of gene names.
#' @param matrix Character, must be one of \code{"auto"}, \code{"count"}, \code{"normalized"}, or \code{"scaled"}. If \code{"auto"}, normalized matrix is used for sections and scaled matrix is used for genes.
#' @param ... Other parameters passed to \code{umap}.
#'
#' @return A \code{SummarizedExperiment} object. The UMAP embeddings are saved in slot \code{meta} if UMAP is performed on sections, or saved in slot \code{gene_embedding} if UMAP is performed on genes.
#'
#' @seealso \code{\link[umap]{umap}} for performing UMAP on a matrix.
#'
#' @importFrom umap umap
#' @importFrom SummarizedExperiment assay rowData<- colData<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- createTomo(zh.data)
#'
#' # Perform UMAP on sections.
#' zh <- runUMAP(zh)
#'
#' # Perform UMAP on some genes.
#' zh <- runUMAP(zh, genes=rownames(zh)[1:100])
runUMAP <- function(object, genes=NA, matrix='auto', ...)
{
    if(all(is.na(genes)))
    {
        if(matrix == 'auto')
            matrix <- 'normalized'
        umap_result <- umap(t(assay(object, matrix)), ...)
        colData(object)$UMAP1 <- umap_result$layout[,1]
        colData(object)$UMAP2 <- umap_result$layout[,2]
        message("UMAP embeddings for sections are saved in column data.\n")
    }
    else
    {
        if(matrix == 'auto')
            matrix <- 'scaled'
        exp_matrix <- assay(object, matrix)[genes, ]
        umap_result <- umap(exp_matrix, ...)
        umap1 <- umap_result$layout[,1]
        umap2 <- umap_result$layout[,2]
        rowData(object)$UMAP1 <- NA
        rowData(object)$UMAP2 <- NA
        rowData(object)[genes, 'UMAP1'] <- umap1
        rowData(object)[genes, 'UMAP2'] <- umap2
        message("UMAP embeddings for genes are saved in row data.\n")
    }
    return(object)
}
