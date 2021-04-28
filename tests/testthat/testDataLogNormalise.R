library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

test_that("dataLogNormalise returns RCA object at default settings", {
    expect_s4_class(dataLogNormalise(rca.obj = pbmc_small_rca), class = "RCA")
    expect_equal(ncol(dataLogNormalise(rca.obj = pbmc_small_rca)$data), ncol(pbmc_small_rca$raw.data))
    expect_equal(nrow(dataLogNormalise(rca.obj = pbmc_small_rca)$data), nrow(pbmc_small_rca$raw.data))
})
