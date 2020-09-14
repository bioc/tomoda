context("Test finding peak genes.")

library(SummarizedExperiment)
data(zh.data)
zh <- CreateTomo(zh.data)
zh <- Normalize(zh)
zh <- Scale(zh)

test_that("find peak",
          {
            expect_equal(FindPeak(c(0:5, 5:0), threshold=1, length=4), c(3,10))
            expect_equal(FindPeak(c(0:2, 2:0), threshold=1, length=4), c(0,0))
          })

test_that("find peak genes",
          {
            peak_genes <- FindPeakGene(zh, nperm=0)
            expect_equal(ncol(peak_genes), 4)
            peak_genes <- FindPeakGene(zh, nperm=1e4)
            expect_equal(ncol(peak_genes), 6)
            peak_genes <- FindPeakGene(zh, length=ncol(zh))
            expect_true(is.null(peak_genes))
          })
