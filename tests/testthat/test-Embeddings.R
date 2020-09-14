context("Test dimensionality reduction methods.")

library(SummarizedExperiment)
data(zh.data)
zh <- CreateTomo(zh.data)
zh <- Normalize(zh)
zh <- Scale(zh)
peak_genes <- FindPeakGene(zh, nperm=1e4)

test_that("PCA for sections",
          {
            zh <- PCA(zh)
            expect_is(zh, 'SummarizedExperiment')
            expect_true(all(c('PC1','PC2') %in% names(colData(zh))))
          })

test_that("PCA for genes",
          {
            zh <- PCA(zh, peak_genes$gene)
            expect_true(all(c('PC1','PC2') %in% names(rowData(zh))))
            expect_equal(sum(!is.na(rowData(zh)$PC1)), nrow(peak_genes))
          })

test_that("TSNE for sections",
          {
            zh <- TSNE(zh)
            expect_is(zh, 'SummarizedExperiment')
            expect_true(all(c('TSNE1','TSNE2') %in% names(colData(zh))))
          })

test_that("TSNE for genes",
          {
            zh <- TSNE(zh, peak_genes$gene)
            expect_true(all(c('TSNE1','TSNE2') %in% names(rowData(zh))))
            expect_equal(sum(!is.na(rowData(zh)$TSNE1)), nrow(peak_genes))
          })

test_that("UMAP for sections",
          {
            zh <- UMAP(zh)
            expect_is(zh, 'SummarizedExperiment')
            expect_true(all(c('UMAP1','UMAP2') %in% names(colData(zh))))
          })

test_that("UMAP for genes",
          {
            zh <- UMAP(zh, peak_genes$gene)
            expect_true(all(c('UMAP1','UMAP2') %in% names(rowData(zh))))
            expect_equal(sum(!is.na(rowData(zh)$UMAP1)), nrow(peak_genes))
          })
