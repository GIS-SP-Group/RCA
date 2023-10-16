library(RCAv2)

load("../testdata/pbmc_small_rca.rda")

pbmc_small_rca <- dataLogNormalise(rca.obj = pbmc_small_rca)

test_that("dataProject returns RCA object", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca)$projection.data), ncol(pbmc_small_rca$data))
})

test_that("dataProject GlobalPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "GlobalPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "GlobalPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "GlobalPanel")$projection.data), 179)
})

test_that("dataProject MonacoPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "MonacoPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "MonacoPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "MonacoPanel")$projection.data), 29)
})

test_that("dataProject MonacoBCellPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "MonacoBCellPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "MonacoBCellPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "MonacoBCellPanel")$projection.data), 5)
})

test_that("dataProject MonacoMonoPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "MonacoMonoPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "MonacoMonoPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "MonacoMonoPanel")$projection.data), 5)
})

test_that("dataProject MonacoTCellPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "MonacoTCellPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "MonacoTCellPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "MonacoTCellPanel")$projection.data), 15)
})

test_that("dataProject CITESeqPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "CITESeqPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "CITESeqPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "CITESeqPanel")$projection.data), 34)
})

test_that("dataProject ENCODEHumanPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "ENCODEHumanPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "ENCODEHumanPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "ENCODEHumanPanel")$projection.data), 93)
})

test_that("dataProject NovershternPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "NovershternPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "NovershternPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "NovershternPanel")$projection.data), 15)
})

test_that("dataProject NovershternTCellPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "NovershternTCellPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "NovershternTCellPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "NovershternTCellPanel")$projection.data), 6)
})

test_that("dataProject ENCODEMousePanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "ENCODEMousePanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "ENCODEMousePanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "ENCODEMousePanel")$projection.data), 15)
})

test_that("dataProject ZhangMouseBrainPanel", {
    expect_s4_class(dataProject(rca.obj = pbmc_small_rca, method = "ZhangMouseBrainPanel"), class = "RCA")
    expect_equal(ncol(dataProject(rca.obj = pbmc_small_rca, method = "ZhangMouseBrainPanel")$projection.data), ncol(pbmc_small_rca$data))
    expect_equal(nrow(dataProject(rca.obj = pbmc_small_rca, method = "ZhangMouseBrainPanel")$projection.data), 7)
})
