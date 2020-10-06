context("Test creating a Tomo object.")

library(SummarizedExperiment)
data(zh.data)
library_size <- apply(zh.data, 2, sum)
normalized <- apply(zh.data, 2, function(x) x * mean(library_size) / sum(x))

test_that("create from count",
          {
            expect_is(createTomo(zh.data), 'SummarizedExperiment')
          })

test_that("create from normalized count",
          {
            expect_is(createTomo(matrix.normalized=normalized), 'SummarizedExperiment')
          })

test_that("create from count and normalized count",
          {
            tomo_object <- createTomo(object=zh.data, matrix.normalized=normalized)
            expect_is(tomo_object, 'SummarizedExperiment')
            expect_equal(assayNames(tomo_object), c('count', 'normalized', 'scaled'))
          })

test_that("create from SummarizeExperiment with count",
          {
              se <- SummarizedExperiment(assays=list(count=zh.data))
              tomo_object <- createTomo(se)
              expect_is(tomo_object, 'SummarizedExperiment')
              expect_equal(assayNames(tomo_object), c('count', 'normalized', 'scaled'))
          })

test_that("create from SummarizeExperiment with normalized count",
          {
              se <- SummarizedExperiment(assays=list(normalized=normalized))
              tomo_object <- createTomo(se)
              expect_is(tomo_object, 'SummarizedExperiment')
          })

test_that("create from SummarizeExperiment with count and normalized count",
          {
              se <- SummarizedExperiment(assays=list(count=zh.data, normalized=normalized))
              tomo_object <- createTomo(se)
              expect_is(tomo_object, 'SummarizedExperiment')
              expect_equal(assayNames(tomo_object), c('count', 'normalized', 'scaled'))
          })
