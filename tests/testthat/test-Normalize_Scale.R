context("Test normalizing and scaling a Tomo object.")

library(SummarizedExperiment)
data(zh.data)
zh <- CreateTomo(zh.data)

test_that("normalize",
          {
            zh_normalized <- Normalize(zh)
            expect_is(zh_normalized, 'SummarizedExperiment')
            expect_equal(names(assays(zh_normalized)), c('count', 'normalized'))
          })

test_that("scale",
          {
            zh_normalized <- Normalize(zh)
            zh_scale <- Scale(zh_normalized)
            expect_is(zh_scale, 'SummarizedExperiment')
            expect_equal(names(assays(zh_scale)), c('count', 'normalized', 'scaled'))
          })

test_that("scale without normalization",
          {
            expect_output(Scale(zh), "Normalized data does not exist! Please run function 'Normalize' before scaling data.")
          })
