#' Line plot for expression traces
#'
#' Plot expression traces for genes across sections in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param genes A character vector of gene names for plotting expression traces.
#' @param matrix Character, must be one of \code{"count"}, \code{"normalized"}, or \code{"scaled"}.
#' @param facet Logical. Plot the expression trace of each gene in a facet if it is \code{TRUE}.
#' @param span Numeric, the amount of smoothing for the default loess smoother. Smaller numbers produce wigglier lines, larger numbers produce smoother lines. Set it to 0 for non-smoothing lines.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{geom_smooth} for plotting smooth lines, \code{facet_wrap} for faceting genes.
#'
#' @importFrom ggplot2 ggplot aes_string theme theme_bw geom_smooth geom_line labs facet_wrap
#' @importFrom SummarizedExperiment rowData colData
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' LinePlot(zh,
#'  c("ENSDARG00000002131", "ENSDARG00000003061", "ENSDARG00000076075", "ENSDARG00000076850"))
#'
#' # Do not smooth lines.
#' LinePlot(zh,
#'  c("ENSDARG00000002131", "ENSDARG00000003061", "ENSDARG00000076075", "ENSDARG00000076850"), span=0)
#'
#' # Plot genes in different facets.
#' LinePlot(zh,
#'  c("ENSDARG00000002131", "ENSDARG00000003061", "ENSDARG00000076075", "ENSDARG00000076850"),
#'  facet=TRUE)
LinePlot <- function(object, genes, matrix='normalized', facet=FALSE, span=0.3)
{
  exp_matrix <- GetMatrix(object, matrix)
  exp_genes <- subset(exp_matrix, rowData(object)$gene %in% genes)
  exp_genes_df <- NULL
  sections <- factor(colData(object)$section, levels=colData(object)$section)
  for(i in 1:nrow(exp_genes))
  {
    exp_gene_df <- data.frame(gene=genes[i], section=sections, value=exp_genes[i,])
    exp_genes_df <- rbind(exp_genes_df, exp_gene_df)
  }
  ylab <- c('Count', 'Normalized count', 'Scaled expression')[c('count', 'normalized', 'scaled') == matrix]

  g <- ggplot(exp_genes_df, aes_string(x='section', y='value', group='gene', color='gene')) +
    labs(x='', y=ylab) +
    theme_bw() +
    theme(legend.position='top', axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

  if(span > 0)
    g <- g + geom_smooth(method='loess', span = span, se=F)
  else
    g <- g + geom_line(size=1)

  if(facet)
    g <- g + facet_wrap(~gene, scales = 'free') + theme(legend.position='none')

  return(g)
}

#' Embedding plot for sections
#'
#' Scatter plot for sections with two-dimenstional embeddings in a \code{SummarizedExperiment} object. Each point stands for a section.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param group Character, a variable in slot \code{meta} defining the groups of sections. Sections in the same group have same colors.
#' @param method Character, the embeddings for scatter plot. Must be one of \code{"TSNE"}, \code{"UMAP"}, or \code{"PCA"}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes_string geom_point theme theme_bw
#' @importFrom ggrepel geom_text_repel
#'
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- TSNE(zh)
#' # Plot TSNE embeddings.
#' EmbedPlot(zh)
#'
#' # Plot UMAP embeddings.
#' zh <- UMAP(zh)
#' EmbedPlot(zh, method="UMAP")
#'
#' # Color sections by kmeans cluster labels.
#' zh <- KMeans(zh, 3)
#' EmbedPlot(zh, group="kmeans_cluster")
EmbedPlot <- function(object, group='section', method='TSNE')
{
  if(!group %in% names(colData(object)))
  {
    group <- 'section'
    cat("Unknown parameter 'group' for plotting sections!")
  }
  colData(object)[[group]] <- as.character(colData(object)[[group]])

  if(method == 'TSNE')
  {
    if(all(c('TSNE1','TSNE2') %in% names(colData(object))))
    {
      g <- ggplot(data.frame(colData(object)), aes_string(x='TSNE1', y='TSNE2', color=group)) +
        geom_point() +
        geom_text_repel(aes_string(label='section')) +
        theme_bw() +
        theme(legend.position='none')
      return(g)
    }
    else
    {
      cat("Function 'TSNE' must be run before plotting TSNE embeddings.\n")
    }
  }
  else if(method == 'UMAP')
  {
    if(all(c('UMAP1','UMAP2') %in% names(colData(object))))
    {
      g <- ggplot(data.frame(colData(object)), aes_string(x='UMAP1', y='UMAP2', color=group)) +
        geom_point(aes_string()) +
        geom_text_repel(aes_string(label='section')) +
        theme_bw() +
        theme(legend.position='none')
      return(g)
    }
    else
    {
      cat("Function 'UMAP' must be run before plotting UMAP embeddings.\n")
    }
  }
  else if(method == 'PCA')
  {
    if(all(c('PC1','PC2') %in% names(colData(object))))
    {
      g <- ggplot(data.frame(colData(object)), aes_string(x='PC1', y='PC2', color=group)) +
        geom_point(aes_string()) +
        geom_text_repel(aes_string(label='section')) +
        theme_bw() +
        theme(legend.position='none')
      return(g)
    }
    else
    {
      cat("Function 'PCA' must be run before plotting PC embeddings.\n")
    }
  }
  else
  {
    cat("Unknown embeddings!\n")
  }
}

#' Embedding plot for genes
#'
#' Scatter plot for genes with two-dimenstional embeddings in a \code{SummarizedExperiment} object. Each point stands for a gene.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param gene.df Data.frame. The first column must be a vector of gene names, and has the name \code{"gene"}. Additional columns in \code{gene.df} can be used to set the colors of genes.
#' @param group Character, a column name in \code{gene.df} defining the groups of genes. Genes in the same group have same colors.
#' @param method Character, the embeddings for scatter plot. Must be one of \code{"TSNE"}, \code{"UMAP"}, or \code{"PCA"}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom SummarizedExperiment assays rowData
#' @importFrom ggplot2 ggplot aes_string geom_point theme theme_bw
#'
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' peak_genes <- FindPeakGene(zh)
#' zh <- TSNE(zh, peak_genes$gene)
#' # Color genes by peak centers.
#' GeneEmbedPlot(zh, peak_genes)
#'
#' # Color genes by peak starts.
#' GeneEmbedPlot(zh, peak_genes, group="start")
#'
#' # Do not color genes.
#' GeneEmbedPlot(zh, peak_genes["gene"])
GeneEmbedPlot <- function(object, gene.df, group='center', method='TSNE')
{
  if(group %in% names(gene.df))
    gene.df[[group]] <- as.factor(gene.df[[group]])
  if(method == 'TSNE')
  {
    if(all(c('TSNE1','TSNE2') %in% names(rowData(object))))
    {
      gene.df$TSNE1 <- rowData(object)[as.character(gene.df$gene), 'TSNE1']
      gene.df$TSNE2 <- rowData(object)[as.character(gene.df$gene), 'TSNE2']
      g <- ggplot(gene.df, aes_string(x='TSNE1', y='TSNE2')) +
        theme_bw()

      if(group %in% names(gene.df))
        g <- g + geom_point(aes_string(color=group))
      else
        g <- g + geom_point()

      return(g)
    }
    else
    {
      cat("Function 'TSNE' must be run for input genes before plotting TSNE embeddings.\n")
    }
  }
  else if(method == 'UMAP')
  {
    if(all(c('UMAP1','UMAP2') %in% names(rowData(object))))
    {
      gene.df$UMAP1 <- rowData(object)[as.character(gene.df$gene), 'UMAP1']
      gene.df$UMAP2 <- rowData(object)[as.character(gene.df$gene), 'UMAP2']
      g <- ggplot(gene.df, aes_string(x='UMAP1', y='UMAP2')) +
        theme_bw()

      if(group %in% names(gene.df))
        g <- g + geom_point(aes_string(color=group))
      else
        g <- g + geom_point()

      return(g)
    }
    else
    {
      cat("Function 'UMAP' must be run for input genes before plotting UMAP embeddings.\n")
    }
  }
  else if(method == 'PCA')
  {
    if(all(c('PC1','PC2') %in% names(rowData(object))))
    {
      gene.df$PC1 <- rowData(object)[as.character(gene.df$gene), 'PC1']
      gene.df$PC2 <- rowData(object)[as.character(gene.df$gene), 'PC2']
      g <- ggplot(gene.df, aes_string(x='PC1', y='PC2')) +
        theme_bw()

      if(group %in% names(gene.df))
        g <- g + geom_point(aes_string(color=group))
      else
        g <- g + geom_point()

      return(g)
    }
    else
    {
      cat("Function 'PCA' must be run for input genes before plotting PC embeddings.\n")
    }
  }
  else
  {
    cat("Unknown embeddings!\n")
  }
}

#' Expression heatmap
#'
#' Heatmap for gene expression across sections in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param genes A vector of character, the name of genes to plot heatmap.
#' @param matrix Character, must be one of \code{"count"}, \code{"normalized"}, or \code{"scaled"}.
#' @param size Character, the size of gene names. Set it to 0 if you do not want to show gene names.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradient2 theme element_blank element_text
#'
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- Scale(zh)
#'
#' # Plot some genes.
#' ExpHeatmap(zh,
#'  c("ENSDARG00000002131", "ENSDARG00000003061", "ENSDARG00000076075", "ENSDARG00000076850"))
#'
#' # Plot peak genes.
#' peak_genes <- FindPeakGene(zh)
#' ExpHeatmap(zh, peak_genes$gene)
#'
#' # Remove gene names if too many genes are in the heatmap.
#' ExpHeatmap(zh, peak_genes$gene, size=0)
ExpHeatmap <- function(object, genes, matrix='scaled', size=5)
{
  genes <- rev(as.character(genes))
  exp_matrix <- GetMatrix(object, matrix)[genes, ]
  if(matrix == 'normalized')
    exp_matrix <- log10(exp_matrix + 1)
  genes_df <- melt(exp_matrix, varnames=c('gene','section'))
  exp_name <-  c('Log10(Count+1)', 'Log10(Normalized count+1)', 'Scaled expression\n(Z score)')[c('count', 'normalized', 'scaled') == matrix]

  g <- ggplot(genes_df, aes_string(x='section', y='gene', fill='value')) +
    geom_raster() +
    scale_fill_gradient2(low='magenta', mid='black', high='yellow', name=exp_name) +
    theme(axis.line=element_blank(),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          axis.ticks=element_blank(),
          panel.background=element_blank())
  if(size == 0)
    g <- g + theme(axis.text.y=element_blank())
  else
    g <- g + theme(axis.text.y=element_text(size=size))

  return(g)
}

#' Correlation heatmap of sections
#'
#' Heatmap pf correlation coefficients between any two sections in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param matrix Character, must be one of \code{"count"}, \code{"normalized"}, or \code{"scaled"}.
#' @param max.cor Numeric, correlation coefficients bigger than \code{max.cor} are set to \code{max.cor}. It is used to clearly show small correlation coefficients.
#' @param cor.method Character, the method to calculate correlation coefficients. must be one of \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom stats cor
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradientn labs theme element_blank element_text
#' @importFrom colorRamps rgb.tables
#'
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- Scale(zh)
#' CorHeatmap(zh)
#'
#' # Use Spearman correlation coefficients.
#' CorHeatmap(zh, cor.method='spearman')
#'
#' # Set max correlation coefficients to 0.3.
#' CorHeatmap(zh, max.cor=0.3)
CorHeatmap <- function(object, matrix='scaled', max.cor=0.5, cor.method='pearson')
{
  exp_matrix <- GetMatrix(object, matrix=matrix)
  cor_matrix <- cor(exp_matrix, method=cor.method)
  cor_matrix[cor_matrix > max.cor] <- max.cor
  cor_df <- melt(cor_matrix, varnames=c('section1','section2'))
  # c('Pearson r', 'Kendall τ', 'Spearman ρ')
  cor_name <- expression("Pearson"~r,"Kendall"~tau,"Spearman"~rho)[c('pearson', 'kendall', 'spearman') == cor.method]

  g <- ggplot(cor_df, aes_string(x='section1', y='section2', fill='value')) +
    geom_raster() +
    scale_fill_gradientn(colors=rgb.tables(10), name=cor_name) +
    labs(x='', y='') +
    theme(axis.line=element_blank(),
          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          axis.ticks=element_blank(),
          panel.background=element_blank())

  return(g)
}

#' Correlation heatmap of genes
#'
#' Heatmap of correlation coefficients between any two queried genes in a \code{SummarizedExperiment} object.
#'
#' @param object A \code{SummarizedExperiment} object.
#' @param gene.df Data.frame. The first column must be a vector of gene names, and has the name \code{"gene"}. Additional columns in \code{gene.df} can be used to set the colors of genes.
#' @param group Character, a column name in \code{gene.df} defining the groups of genes. Genes in the same group have same colors on the side bar.
#' @param matrix Character, must be one of \code{"count"}, \code{"normalized"}, or \code{"scaled"}.
#' @param size Numeric, the size of gene names. Set it to 0 if you do not want to show gene names.
#' @param cor.method Character, the method to calculate correlation coefficients. must be one of \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}.
#'
#' @details This method can create a pure heatmap or a heatmap with side bar. If you prefer a pure heatmap, input a \code{gene.df} with a single column of gene names.
#' However, you may want to show additional information of genes with a side bar, and the grouping information should be saved as additional column(s) of \code{gene.df}, and declared as \code{group}.
#' By default, you can use the output by \code{FindPeakGene} as input \code{gene.df}. Peak genes will be grouped by their centers on the side bar.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom stats cor
#' @importFrom grDevices rainbow
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradientn labs theme element_blank element_text ggplot_build scale_color_manual guides annotation_raster guide_legend coord_cartesian
#'
#' @export
#'
#' @examples
#' data(zh.data)
#' zh <- CreateTomo(zh.data)
#' zh <- Normalize(zh)
#' zh <- Scale(zh)
#'
#' # Correlation heatmap for all peak genes.
#' peak_genes <- FindPeakGene(zh)
#' GeneCorHeatmap(zh, peak_genes)
#'
#' # Use Spearman correlation coefficients.
#' GeneCorHeatmap(zh, peak_genes, cor.method="spearman")
#'
#' # Group genes by peak start.
#' GeneCorHeatmap(zh, peak_genes, group="start")
#'
#' # Plot without side bar.
#' GeneCorHeatmap(zh, data.frame(
#'  gene=c("ENSDARG00000002131", "ENSDARG00000003061", "ENSDARG00000076075", "ENSDARG00000076850")))
GeneCorHeatmap <- function(object, gene.df, group='center', matrix='scaled', size=5, cor.method='pearson')
{
  exp_matrix <- GetMatrix(object, matrix=matrix)
  genes <- as.character(gene.df$gene)
  cor_matrix <- cor(t(exp_matrix[genes, ]), method=cor.method)
  cor_df <- melt(cor_matrix, varnames=c('gene1', 'gene2'))
  cor_name <- expression("Pearson"~r,"Kendall"~tau,"Spearman"~rho)[c('pearson', 'kendall', 'spearman') == cor.method]
  # plot sidebar if genes are grouped
  plot_sidebar <- group %in% names(gene.df)

  if(plot_sidebar)
  {
    group_genes <- gene.df[[group]]
    names(group_genes) <- gene.df$gene
    cor_df$group <- as.factor(group_genes[cor_df$gene1])
    # Colors for column side bar
    n_groups <- length(unique(group_genes))
    color_legend <- rainbow(n_groups)
    names(color_legend) <- unique(group_genes)
  }

  g <- ggplot(cor_df, aes_string(x='gene1', y='gene2', fill='value')) +
    geom_raster() +
    scale_fill_gradientn(colors=rgb.tables(10), name=cor_name)  +
    labs(x='', y='') +
    theme(axis.line=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_blank(),
          legend.key=element_blank(),
          legend.box='horizontal')
  if(size == 0)
    g <- g + theme(axis.text.x=element_blank())
  else
    g <- g + theme(axis.text.x=element_text(angle=90, hjust=1, size=size))

  if(plot_sidebar)
  {
    gbuild <- ggplot_build(plot = g)
    y_range <- diff(x = gbuild$layout$panel_params[[1]]$y.range)
    y_min <- max(gbuild$layout$panel_params[[1]]$y.range) + 0.01 * y_range
    y_max <- y_min + 0.03 * y_range

    g <- g + geom_point(aes_string(color='group'), alpha=0, shape=15, size=5) +
      scale_color_manual(name=group, values=color_legend)  +
      guides(color=guide_legend(override.aes=list(alpha=1))) +
      annotation_raster(t(color_legend[as.character(group_genes)]), -Inf, Inf, ymin=y_min, ymax=y_max) +
      coord_cartesian(ylim = c(0, y_max), clip='off')
  }

  return(g)
}
