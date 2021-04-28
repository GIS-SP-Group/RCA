library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

test_that("dataFilter returns RCA object even if no cells are filtered out", {
    expect_s4_class(dataFilter(rca.obj = pbmc_small_rca, nGene.thresholds = c(0,Inf), nUMI.thresholds = c(0,Inf), percent.mito.thresholds = c(0,Inf)), class = "RCA")
    expect_equal(ncol(dataFilter(rca.obj = pbmc_small_rca, nGene.thresholds = c(0,Inf), nUMI.thresholds = c(0,Inf), percent.mito.thresholds = c(0,Inf))$raw.data), 80)
})

test_that("dataFilter returns RCA object at default settings", {
    expect_s4_class(dataFilter(rca.obj = pbmc_small_rca), class = "RCA")
})


test_that("dataFilter returns RCA object even if all cells are filtered out", {
    expect_s4_class(dataFilter(rca.obj = pbmc_small_rca, nGene.thresholds = c(0,0)), class = "RCA")
    expect_equal(ncol(dataFilter(rca.obj = pbmc_small_rca, nGene.thresholds = c(0,Inf), nUMI.thresholds = c(0,Inf), percent.mito.thresholds = c(0,Inf))$raw.data), 0)
})
