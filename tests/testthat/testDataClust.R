library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

pbmc_small_rca <- dataLogNormalise(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataProject(rca.obj = pbmc_small_rca)

test_that("dataClust returns RCA object at default settings", {
    expect_s4_class(dataClust(rca.obj = pbmc_small_rca), class = "RCA")
    expect_equal(length(dataClust(rca.obj = pbmc_small_rca)$clustering.out$dynamicColorsList$`deepSplit 1`), ncol(pbmc_small_rca$projection.data))
})

test_that("dataSClust returns RCA object at default settings", {
    expect_s4_class(dataSClust(rca.obj = pbmc_small_rca, approx = FALSE), class = "RCA")
    expect_equal(length(dataSClust(rca.obj = pbmc_small_rca, approx = FALSE)$clustering.out$dynamicColorsList$Clusters), ncol(pbmc_small_rca$projection.data))
})

test_that("dataSNN returns RCA object at default settings", {
    expect_s4_class(dataSNN(rca.obj = pbmc_small_rca), class = "RCA")
    expect_equal(length(dataSNN(rca.obj = pbmc_small_rca)$clustering.out$dynamicColorsList$Clusters), ncol(pbmc_small_rca$projection.data))
})
