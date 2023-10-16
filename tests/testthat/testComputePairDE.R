library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

pbmc_small_rca <- dataLogNormalise(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataProject(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataClust(rca.obj = pbmc_small_rca)

test_that("dataDE returns RCA object at default settings", {
    expect_s4_class(dataDE(rca.obj = pbmc_small_rca), class = "RCA")
    expect_equal(length(dataDE(rca.obj = pbmc_small_rca)$DE.genes), 2)
    expect_s3_class(dataDE(rca.obj = pbmc_small_rca)$DE.genes[[1]], class = "data.frame")
    expect_s3_class(dataDE(rca.obj = pbmc_small_rca)$DE.genes[[2]], class = "data.frame")
})
