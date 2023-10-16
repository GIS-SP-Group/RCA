library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

pbmc_small_rca <- dataLogNormalise(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataProject(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataSNN(rca.obj = pbmc_small_rca)
res <- parameterSpaceSNN(rca.obj = pbmc_small_rca)

test_that("parameterSpaceSNN returns data frame at default settings", {
    expect_type(parameterSpaceSNN(rca.obj = pbmc_small_rca), type = "list")
    expect_length(parameterSpaceSNN(rca.obj = pbmc_small_rca), n = 2)

    expect_s3_class(parameterSpaceSNN(rca.obj = pbmc_small_rca)[[1]], class = "data.frame")
    expect_equal(nrow(parameterSpaceSNN(rca.obj = pbmc_small_rca)[[1]]), 2016)
    expect_equal(ncol(parameterSpaceSNN(rca.obj = pbmc_small_rca)[[1]]), 4)
    expect_equal(colnames(parameterSpaceSNN(rca.obj = pbmc_small_rca)[[1]]), c("K", "eps", "minPts", "clusters"))

    expect_s3_class(parameterSpaceSNN(rca.obj = pbmc_small_rca)[[2]], class = "plotly")
})

test_that("parameterSpaceSeurat returns RCA object at default settings", {
    expect_type(parameterSpaceSeurat(rca.obj = pbmc_small_rca), type = "list")
    expect_length(parameterSpaceSeurat(rca.obj = pbmc_small_rca), n = 2)

    expect_s3_class(parameterSpaceSeurat(rca.obj = pbmc_small_rca)[[1]], class = "data.frame")
    expect_equal(nrow(parameterSpaceSeurat(rca.obj = pbmc_small_rca)[[1]]), 11)
    expect_equal(ncol(parameterSpaceSeurat(rca.obj = pbmc_small_rca)[[1]]), 2)
    expect_equal(colnames(parameterSpaceSeurat(rca.obj = pbmc_small_rca)[[1]]), c("Resolution", "Clusters"))

    expect_s3_class(parameterSpaceSeurat(rca.obj = pbmc_small_rca)[[2]], class = "ggplot")
})
