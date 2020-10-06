context("Test visualization functions.")

data(zh.data)
zh <- createTomo(zh.data)
zh <- runPCA(zh)
zh <- runTSNE(zh)
zh <- runUMAP(zh)

peak_genes <- findPeakGene(zh, nperm=1e4)
zh <- runPCA(zh, peak_genes$gene)
zh <- runTSNE(zh, peak_genes$gene)
zh <- runUMAP(zh, peak_genes$gene)

test_that("embedding plot for sections",
          {
            expect_is(embedPlot(zh), 'ggplot')
            expect_is(embedPlot(zh, method='PCA'), 'ggplot')
            expect_is(embedPlot(zh, method='UMAP'), 'ggplot')
          })

test_that("embedding plot for genes",
          {
            expect_is(geneEmbedPlot(zh, peak_genes), 'ggplot')
            expect_is(geneEmbedPlot(zh, peak_genes, method='PCA'), 'ggplot')
            expect_is(geneEmbedPlot(zh, peak_genes, method='UMAP'), 'ggplot')
          })

test_that("expression heatmap",
          {
            expect_is(expHeatmap(zh, peak_genes$genes), 'ggplot')
          })

test_that("correlation heatmap for sections",
          {
            expect_is(corHeatmap(zh), 'ggplot')
          })

test_that("correlation heatmap for genes",
          {
            expect_is(geneCorHeatmap(zh, peak_genes), 'ggplot')
          })
