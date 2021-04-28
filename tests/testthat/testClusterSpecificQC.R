library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

pbmc_small_rca <- dataLogNormalise(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataProject(rca.obj = pbmc_small_rca)
pbmc_small_rca <- dataClust(rca.obj = pbmc_small_rca)

test_that("performClusterSpecificQC returns RCA object when no cells are filtered out", {
    expect_s4_class(performClusterSpecificQC(rca.obj = pbmc_small_rca, cluster.labels = pbmc_small_rca$clustering.out$dynamicColorsList$`deepSplit 1`, nGene.low.thresholds = c(0,0,0,0,0,0), nGene.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf), nUMI.low.thresholds = c(0,0,0,0,0,0), nUMI.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf), pMito.low.thresholds = c(0,0,0,0,0,0), pMito.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf)), class = "RCA")
    expect_equal(ncol(performClusterSpecificQC(rca.obj = pbmc_small_rca, cluster.labels = pbmc_small_rca$clustering.out$dynamicColorsList$`deepSplit 1`, nGene.low.thresholds = c(0,0,0,0,0,0), nGene.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf), nUMI.low.thresholds = c(0,0,0,0,0,0), nUMI.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf), pMito.low.thresholds = c(0,0,0,0,0,0), pMito.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf))$raw.data), ncol(pbmc_small_counts))
    expect_equal(ncol(performClusterSpecificQC(rca.obj = pbmc_small_rca, cluster.labels = pbmc_small_rca$clustering.out$dynamicColorsList$`deepSplit 1`, nGene.low.thresholds = c(0,0,0,0,0,0), nGene.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf), nUMI.low.thresholds = c(0,0,0,0,0,0), nUMI.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf), pMito.low.thresholds = c(0,0,0,0,0,0), pMito.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf))$data), ncol(pbmc_small_counts))
    expect_equal(ncol(performClusterSpecificQC(rca.obj = pbmc_small_rca, cluster.labels = pbmc_small_rca$clustering.out$dynamicColorsList$`deepSplit 1`, nGene.low.thresholds = c(0,0,0,0,0,0), nGene.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf), nUMI.low.thresholds = c(0,0,0,0,0,0), nUMI.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf), pMito.low.thresholds = c(0,0,0,0,0,0), pMito.high.thresholds = c(Inf,Inf,Inf,Inf,Inf,Inf))$projection.data), ncol(pbmc_small_counts))
})
