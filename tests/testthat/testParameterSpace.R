library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

pbmc_small_rca <- dataLogNormalise(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataProject(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataSNN(rca.obj = pbmc_small_rca)
res <- parameterSpaceSNN(rca.obj = pbmc_small_rca)

test_that("parameterSpaceSNN returns data frame at default settings", {
    expect_s3_class(parameterSpaceSNN(rca.obj = pbmc_small_rca), class = "data.frame")
    expect_equal(nrow(parameterSpaceSNN(rca.obj = pbmc_small_rca)), 2016)
    expect_equal(ncol(parameterSpaceSNN(rca.obj = pbmc_small_rca)), 4)
    expect_equal(colnames(parameterSpaceSNN(rca.obj = pbmc_small_rca)), c("K", "eps", "minPts", "clusters"))
})

test_that("parameterSpaceSeurat returns RCA object at default settings", {
   expect_error(parameterSpaceSeurat(rca.obj = pbmc_small_rca))
})
