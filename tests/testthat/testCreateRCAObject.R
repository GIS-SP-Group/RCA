library(RCAv2)

mat <- matrix(sample(1000000, size = 230*80), nrow = 230, ncol = 80)

test_that("createRCAObject throws errors for non dgCMatrix inputs", {
    expect_error(createRCAObject(rawData = mat))
    expect_error(createRCAObject(rawData = as.data.frame(mat)))
    expect_s4_class(createRCAObject(rawData = as(mat, "dgCMatrix")), class = "RCA")
})
