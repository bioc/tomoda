context("Test normalizing and scaling a Tomo object.")

library(SummarizedExperiment)
data(zh.data)
zh <- createTomo(zh.data, normalize=FALSE, scale=FALSE)

test_that("normalize",
          {
            zh_normalized <- normalizeTomo(zh)
            expect_is(zh_normalized, 'SummarizedExperiment')
            expect_equal(assayNames(zh_normalized), c('count', 'normalized'))
          })

test_that("scale",
          {
            zh_normalized <- normalizeTomo(zh)
            zh_scale <- scaleTomo(zh_normalized)
            expect_is(zh_scale, 'SummarizedExperiment')
            expect_equal(assayNames(zh_scale), c('count', 'normalized', 'scaled'))
          })

test_that("scale without normalization",
          {
            expect_identical(scaleTomo(zh), zh)
          })
