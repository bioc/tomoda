context("Test visualization functions.")

data(zh.data)
zh <- CreateTomo(zh.data)
zh <- Normalize(zh)
zh <- Scale(zh)
zh <- PCA(zh)
zh <- TSNE(zh)
zh <- UMAP(zh)

peak_genes <- FindPeakGene(zh, nperm=1e4)
zh <- PCA(zh, peak_genes$gene)
zh <- TSNE(zh, peak_genes$gene)
zh <- UMAP(zh, peak_genes$gene)

test_that("embedding plot for sections",
          {
            expect_is(EmbedPlot(zh), 'ggplot')
            expect_is(EmbedPlot(zh, method='PCA'), 'ggplot')
            expect_is(EmbedPlot(zh, method='UMAP'), 'ggplot')
          })

test_that("embedding plot for genes",
          {
            expect_is(GeneEmbedPlot(zh, peak_genes), 'ggplot')
            expect_is(GeneEmbedPlot(zh, peak_genes, method='PCA'), 'ggplot')
            expect_is(GeneEmbedPlot(zh, peak_genes, method='UMAP'), 'ggplot')
          })

test_that("expression heatmap",
          {
            expect_is(ExpHeatmap(zh, peak_genes$genes), 'ggplot')
          })

test_that("correlation heatmap for sections",
          {
            expect_is(CorHeatmap(zh), 'ggplot')
          })

test_that("correlation heatmap for genes",
          {
            expect_is(GeneCorHeatmap(zh, peak_genes), 'ggplot')
          })
