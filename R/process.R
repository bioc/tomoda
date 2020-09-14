#' Normalize data
#'
#' Normalize the raw read count in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#'
#' @details This function should be run for \code{SummarizedExperiment} object created from raw read count matrix.
#' If the \code{SummarizedExperiment} object already has a normalized count matrix. The function simply return the original object.
#' Library sizes of all sections are normalized to the median library size.
#'
#' @return A \code{SummarizedExperiment} object with normalized read count matrix saved in assay \code{'normalized'}.
#'
#' @importFrom stats median
#' @importFrom SummarizedExperiment assays assays<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
Normalize <- function(object)
{
  if("normalized" %in% names(assays(object)))
    cat('Normalized data already exists!')
  else
  {
    library_size <- apply(assays(object)$count, 2, sum)
    assays(object)$normalized <- apply(assays(object)$count, 2, function(x) x * median(library_size) / sum(x))
    cat("Normalized count matrix is saved in assay 'normalized'.\n")
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
#' @importFrom SummarizedExperiment assays assays<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- Scale(zh)
Scale <- function(object)
{
  if(!"normalized" %in% names(assays(object)))
    cat("Normalized data does not exist! Please run function 'Normalize' before scaling data.\n")
  else
  {
    mean_normalized <- apply(assays(object)$normalized, 1, mean)
    sd_normalized <- apply(assays(object)$normalized, 1, sd)
    assays(object)$scaled <- (assays(object)$normalized - mean_normalized) / sd_normalized
    cat("Scaled count matrix is saved in assay 'scaled'.\n")
  }
  return(object)
}

#' Get a matrix from \code{SummarizedExperiment} object
#'
#' Get one of three stored matrics, \code{count}, \code{normalized}, or \code{scaled}, in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param matrix Character, must be one of \code{"count"}, \code{"normalized"}, or \code{"scaled"}.
#'
#' @return A numeric matrix. Returns \code{NA} if such matrix does not exist.
#'
#' @importFrom SummarizedExperiment assays
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#'
#' # Get the raw read count matrix.
#' matrix_count <- GetMatrix(zh, "count")
#'
#' # Get the normalized read count matrix.
#' zh <- Normalize(zh)
#' matrix_normalized <- GetMatrix(zh, "normalized")
#'
#' # Get the scaled read count matrix.
#' zh <- Scale(zh)
#' matrix_scaled <- GetMatrix(zh, "scaled")
GetMatrix <- function(object, matrix)
{
  if(matrix == 'count')
  {
    if("count" %in% names(assays(object)))
      return(assays(object)$count)
    else
    {
      cat('Count matrix does not exist!\n')
      return(NULL)
    }
  }

  else if(matrix == 'normalized')
  {
    if("normalized" %in% names(assays(object)))
      return(assays(object)$normalized)
    else
    {
      cat('Normalized count matrix does not exist!\n')
      return(NULL)
    }
  }

  else if(matrix == 'scaled')
  {
    if("scaled" %in% names(assays(object)))
      return(assays(object)$scaled)
    else
    {
      cat('Scaled count matrix does not exist!\n')
      return(NULL)
    }
  }

  else
  {
    cat('Unknown matrix', matrix,'. Please input correct name.\n')
    return(NULL)
  }
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
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' hclust_zh <- HClust(zh)
#' plot(hclust_zh)
#'
#' # Use other agglomeration method
#' hclust_zh <- HClust(zh, agglomeration="average")
#'
#' # (Not recommended) Use scaled read counts to calculate distance
#' zh <- Scale(zh)
#' hclust_zh <- HClust(zh, matrix="scaled")
HClust <-  function(object, matrix='normalized', measure='euclidean', p=2, agglomeration='complete')
{
  exp_matrix <- t(GetMatrix(object, matrix))
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
#' @importFrom SummarizedExperiment colData colData<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- KMeans(zh, 3)
#'
#' # Use scaled read counts to calculate distance
#' zh <- Scale(zh)
#' zh <- KMeans(zh, 3, matrix="scaled")
KMeans <- function(object, centers, matrix='normalized', ...)
{
  exp_matrix <- t(GetMatrix(object, matrix))
  kmeans_section <- kmeans(exp_matrix,centers=centers, ...)
  percent_between <- kmeans_section$betweenss / kmeans_section$totss
  colData(object)$kmeans_cluster <- kmeans_section$cluster

  cluster_list <- list()
  for(i in 1:centers)
    cluster_list[[i]] <- as.character(colData(object)$section)[kmeans_section$cluster == i]

  cat("KMeans results:\n")
  print(cluster_list)
  cat("between_SS / total_SS =", percent_between,'\n')
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
#' FindPeak(c(0:5, 5:0), threshold=1, length=4)
#'
#' # Most likely return c(0, 0)
#' FindPeak(rnorm(10), threshold=3, length=3)
FindPeak <- function(x, threshold=1, length=4)
{
  gt_threshold <- x > threshold
  consecutive <- rep(0, length(x) + 1)

  for(i in 1:length(x) )
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
#' @importFrom SummarizedExperiment assays rowData
#' @importFrom stats p.adjust
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- Scale(zh)
#' peak_genes <- FindPeakGene(zh)
#' head(peak_genes)
#'
#' # Increase threshold so that less peak genes will be found.
#' peak_genes <- FindPeakGene(zh, threshold=1.5)
#'
#' # Increase peak length so that less peak genes will be found.
#' peak_genes <- FindPeakGene(zh, length=5)
#'
#' # Set nperm to 0 so that p values will not be calculated. This will save running time.
#' peak_genes <- FindPeakGene(zh, nperm=0)
FindPeakGene <- function(object, threshold=1, length=4, nperm=1e5, method='BH')
{
  if(!"scaled" %in% names(assays(object)))
  {
    cat("Scaled data does not exist! Please run function 'Scale' before finding peak genes.\n")
    return()
  }
  else
  {
    scaled <- assays(object)$scaled
    peak_position <- apply(scaled, 1, FindPeak, threshold, length)
    peak_exist <- peak_position[1,] != 0
    if(!any(peak_exist))
    {
      cat("No peak gene is found!\n")
      return()
    }

    peak_genes <- rowData(object)$gene[peak_exist]
    peak_gene_df <- data.frame(gene=peak_genes,
                               start=peak_position[1,peak_exist],
                               end=peak_position[2,peak_exist],
                               center=floor(apply(peak_position[,peak_exist], 2, mean)),
                               stringsAsFactors=FALSE)

    # Using permutation to calculate p-values
    if(nperm > 0)
    {
      pvals <- NULL
      n_section <- nrow(object)
      saved_pvals <- rep(NA, n_section)
      #saved_pvals[length] <- (n_section - length + 1) / choose(n_section, length)
      for(gene in peak_genes)
      {
        exp_gene <- scaled[gene, ]
        n_gt_threshold <- sum(exp_gene > threshold)
        if(is.na(saved_pvals[n_gt_threshold]))
        {
          n_peak <- 0
          for(i in 1:nperm)
            n_peak <- n_peak + (FindPeak(sample(exp_gene))[1] > 0)
          pval <- n_peak / nperm
          saved_pvals[n_gt_threshold] <- pval
        }
        else
        {
          pval <- saved_pvals[n_gt_threshold]
        }
        pvals <- c(pvals, pval)
      }

      peak_gene_df$p=pvals
      peak_gene_df$p.adj=p.adjust(pvals, method=method)
    }

    sorted_df <- peak_gene_df[order(peak_gene_df$start, peak_gene_df$end), ]
    cat(nrow(sorted_df), "peak genes (spatially upregulated genes) are found!\n")
    return(sorted_df)
  }
}

#' Perform PCA
#'
#' Perform PCA on sections or genes in a \code{SummarizedExperiment} object for dimensionality reduction.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param genes \code{NA} or a vector of character. Perform PCA on sections if it is \code{NA}, or on given genes if it is a vector of gene names.
#' @param scree Logical, plot the scree plot for PCs if it is \code{TRUE}.
#' @param ... Other parameters passed to \code{prcomp}.
#'
#' @return A \code{SummarizedExperiment} object. The PC embeddings are saved in slot \code{meta} if PCA is performed on sections, or saved in slot \code{gene_embedding} if PCA is performed on genes.
#'
#' @seealso \code{\link[stats]{prcomp}} for performing PCA on a matrix.
#'
#' @importFrom stats prcomp
#' @importFrom SummarizedExperiment assays rowData<- colData<-
#' @importFrom ggplot2 ggplot aes_string geom_point labs theme_bw
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- Scale(zh)
#'
#' # Perform PCA on sections.
#' zh <- PCA(zh)
#'
#' # Plot the scree plot.
#' zh <- PCA(zh, scree=TRUE)
#'
#' # Perform PCA on some genes.
#' zh <- PCA(zh, genes=rownames(zh)[1:100])
PCA <- function(object, genes=NA, scree=FALSE, ...)
{
  if(all(is.na(genes)))
  {
    pca_result <- prcomp(assays(object)$normalized, ...)
    colData(object)$PC1 <- pca_result$rotation[,1]
    colData(object)$PC2 <- pca_result$rotation[,2]
    cat("PC embeddings for sections are saved in column data.\n")
  }
  else
  {
    exp_matrix <- assays(object)$scaled[genes, ]
    pca_result <- prcomp(t(exp_matrix))
    pc1 <- pca_result$rotation[,1]
    pc2 <- pca_result$rotation[,2]
    rowData(object)$PC1 <- NA
    rowData(object)$PC2 <- NA
    rowData(object)[genes, 'PC1'] <- pc1
    rowData(object)[genes, 'PC2'] <- pc2
    cat("PC embeddings for genes are saved in row data.\n")
  }

  if(scree)
  {
    pca_sd <- data.frame(pc=1:length(pca_result$sdev), sd=pca_result$sdev)
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
#' @param perplexity Numeric, perplexity parameter for Rtsne (default: 0.25 *(number of observations - 1)).
#' @param ... Other parameters passed to \code{Rtsne}.
#'
#' @return A \code{SummarizedExperiment} object. The TSNE embeddings are saved in slot \code{meta} if TSNE is performed on sections, or saved in slot \code{gene_embedding} if TSNE is performed on genes.
#'
#' @seealso \code{\link[Rtsne]{Rtsne}} for performing PCA on a matrix.
#'
#' @importFrom Rtsne Rtsne
#' @importFrom SummarizedExperiment assays rowData<- colData<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- Scale(zh)
#'
#' # Perform TSNE on sections.
#' zh <- TSNE(zh)
#'
#' # Perform TSNE on sections with other perplexity.
#' zh <- TSNE(zh, perplexity=10)
#'
#' # Perform TSNE on some genes.
#' zh <- TSNE(zh, genes=rownames(zh)[1:100])
TSNE <- function(object, genes=NA, perplexity=NA, ...)
{
  if(all(is.na(genes)))
  {
    if(is.na(perplexity))
      perplexity <- (ncol(object) - 1) / 4
    tsne_result <- Rtsne(t(assays(object)$normalized), perplexity=perplexity, ...)
    colData(object)$TSNE1 <- tsne_result$Y[,1]
    colData(object)$TSNE2 <- tsne_result$Y[,2]
    cat("TSNE embeddings for sections are saved in column data.\n")
  }
  else
  {
    if(is.na(perplexity))
      perplexity <- (length(genes) - 1) / 4
    exp_matrix <- assays(object)$scaled[genes, ]
    tsne_result <- Rtsne(exp_matrix, perplexity=perplexity, ...)
    tsne1 <- tsne_result$Y[,1]
    tsne2 <- tsne_result$Y[,2]
    #names(tsne1) <- rownames(exp_matrix)
    #names(tsne2) <- rownames(exp_matrix)
    rowData(object)$TSNE1 <- NA
    rowData(object)$TSNE2 <- NA
    rowData(object)[genes, 'TSNE1'] <- tsne1
    rowData(object)[genes, 'TSNE2'] <- tsne2
    cat("TSNE embeddings for genes are saved in row data.\n")
  }
  return(object)
}

#' Perform UMAP
#'
#' Perform UMAP on sections or genes in a \code{SummarizedExperiment} object for dimensionality reduction.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param genes \code{NA} or a vector of character. Perform UMAP on sections if it is \code{NA}, or on given genes if it is a vector of gene names.
#' @param ... Other parameters passed to \code{umap}.
#'
#' @return A \code{SummarizedExperiment} object. The UMAP embeddings are saved in slot \code{meta} if UMAP is performed on sections, or saved in slot \code{gene_embedding} if UMAP is performed on genes.
#'
#' @seealso \code{\link[umap]{umap}} for performing UMAP on a matrix.
#'
#' @importFrom umap umap
#' @importFrom SummarizedExperiment assays rowData<- colData<-
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- Scale(zh)
#'
#' # Perform TSNE on sections.
#' zh <- UMAP(zh)
#'
#' # Perform TSNE on some genes.
#' zh <- UMAP(zh, genes=rownames(zh)[1:100])
UMAP <- function(object, genes=NA, ...)
{
  if(all(is.na(genes)))
  {
    umap_result <- umap(t(assays(object)$normalized), ...)
    colData(object)$UMAP1 <- umap_result$layout[,1]
    colData(object)$UMAP2 <- umap_result$layout[,2]
    cat("UMAP embeddings for sections are saved in column data.\n")
  }
  else
  {
    exp_matrix <- assays(object)$scaled[genes, ]
    umap_result <- umap(exp_matrix, ...)
    umap1 <- umap_result$layout[,1]
    umap2 <- umap_result$layout[,2]
    #names(umap1) <- rownames(exp_matrix)
    #names(umap2) <- rownames(exp_matrix)
    rowData(object)$UMAP1 <- NA
    rowData(object)$UMAP2 <- NA
    rowData(object)[genes, 'UMAP1'] <- umap1
    rowData(object)[genes, 'UMAP2'] <- umap2
    cat("UMAP embeddings for genes are saved in row data.\n")
  }
  return(object)
}
