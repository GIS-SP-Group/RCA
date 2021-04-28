library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

pbmc_small_rca <- dataLogNormalise(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataProject(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataClust(rca.obj = pbmc_small_rca)

test_that("computeUMAP returns RCA object at default settings", {
    expect_s4_class(computeUMAP(rca.obj = pbmc_small_rca), class = "RCA")
    expect_equal(nrow(computeUMAP(rca.obj = pbmc_small_rca)$umap.coordinates), ncol(pbmc_small_rca$projection.data))
})

