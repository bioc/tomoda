context("Test creating a Tomo object.")

library(SummarizedExperiment)
data(zh.data)
library_size <- apply(zh.data, 2, sum)
normalized <- apply(zh.data, 2, function(x) x * mean(library_size) / sum(x))

test_that("create from count",
          {
            expect_is(CreateTomo(zh.data), 'SummarizedExperiment')
          })

test_that("create from normalized count",
          {
            expect_is(CreateTomo(matrix.normalized=normalized), 'SummarizedExperiment')
          })

test_that("create from count and normalized count",
          {
            tomo_object <- CreateTomo(matrix.count=zh.data, matrix.normalized=normalized)
            expect_is(tomo_object, 'SummarizedExperiment')
            expect_equal(names(assays(tomo_object)), c('count', 'normalized'))
          })
