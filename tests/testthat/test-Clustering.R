context("Test clustering.")

library(SummarizedExperiment)
data(zh.data)
zh <- createTomo(zh.data)

test_that("hierarchical clustering",
          {
            expect_is(hierarchClust(zh), 'hclust')
          })

test_that("k-means",
          {
            zh_kmeans <- kmeansClust(zh, 3)
            expect_is(zh_kmeans, 'SummarizedExperiment')
            expect_true('kmeans_cluster' %in% names(colData(zh_kmeans)))
          })
