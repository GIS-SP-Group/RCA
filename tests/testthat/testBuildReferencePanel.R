library(RCAv2)
library(scuttle)

ct_num = 1
ct_colnames <- c()
for(i in 1:80) {
    ct_colnames <- append(ct_colnames, paste0("Celltype", ct_num, "_", (i %% 20)))
    if(i %% 20 == 0)
        ct_num = ct_num + 1
}

test_that("buildReferencePanel throws errors for random matrix", {
    mat <- matrix(sample(1000000, size = 230*80), nrow = 230, ncol = 80)

    colnames(mat) <- ct_colnames
    expect_error(buildReferencePanel(bulk.rna.data = mat, verbose = TRUE))
})

test_that("buildReferencePanel works for RNA matrix", {
    pbmc_tpm <- scuttle::calculateTPM(x = as.matrix(pbmc_small_counts))/1000
    colnames(pbmc_tpm) <- ct_colnames
    expect_equal(ncol(buildReferencePanel(bulk.rna.data = pbmc_tpm, verbose = FALSE)), 4)
    expect_equal(colnames(buildReferencePanel(bulk.rna.data = pbmc_tpm, verbose = FALSE)), c("Celltype1", "Celltype2", "Celltype3", "Celltype4"))
})
