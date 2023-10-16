library(RCAv2)

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
