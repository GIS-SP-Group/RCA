library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

pbmc_small_rca <- dataLogNormalise(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataProject(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataClust(rca.obj = pbmc_small_rca)

test_that("estimateCellTypeFromProjection returns RCA object at default settings", {
    expect_s4_class(estimateCellTypeFromProjection(rca.obj = pbmc_small_rca), class = "RCA")
    expect_length(unlist(estimateCellTypeFromProjection(rca.obj = pbmc_small_rca)$cell.Type.Estimate), n = ncol(pbmc_small_rca$projection.data))
})

test_that("estimateCellTypeFromProjectionPerCluster returns RCA object at default settings", {
    expect_s4_class(estimateCellTypeFromProjectionPerCluster(rca.obj = pbmc_small_rca), class = "RCA")
    expect_length(unlist(estimateCellTypeFromProjectionPerCluster(rca.obj = pbmc_small_rca)$cell.Type.Estimate), n = ncol(pbmc_small_rca$projection.data))
    expect_equal(names(unlist(estimateCellTypeFromProjectionPerCluster(rca.obj = pbmc_small_rca)$cell.Type.Estimate)), pbmc_small_rca$clustering.out$dynamicColorsList$`deepSplit 1`)
})
