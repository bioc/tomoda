context("Test clustering.")

library(SummarizedExperiment)
data(zh.data)
zh <- CreateTomo(zh.data)
zh <- Normalize(zh)

test_that("hierarchical clustering",
          {
            expect_is(HClust(zh), 'hclust')
          })

test_that("k-means",
          {
            zh_kmeans <- KMeans(zh, 3)
            expect_is(zh_kmeans, 'SummarizedExperiment')
            expect_true('kmeans_cluster' %in% names(colData(zh_kmeans)))
          })
