#' Log-transform and normalise data by sequencing depth
#'
#' @param rca.obj RCA object.
#' @param scale.factor scaling factor for log-normalization
#' @return RCA object.
#' @export
#'
dataLogNormalise <- function(rca.obj, scale.factor = 10000) {

    # Extract data from RCA object
    data <- rca.obj$raw.data

    # Compute sequencing depth vector
    seqDepthVec <- Matrix::colSums(data)

    # Normalise data by cell
    pb <- txtProgressBar(style = 3)

    norm.data <- sapply(seq_along(seqDepthVec), function(index) {
        norm <- scale.factor*data[,index]/seqDepthVec[index]
        setTxtProgressBar(pb = pb, value = index/length(seqDepthVec))
        return(norm)
    })
    colnames(norm.data) <- colnames(data)

    # Log-transform normalised data
    logNorm.data <- log(1+norm.data)

    # Store log-transformed normalised data in RCA object
    rca.obj$data <- as(logNorm.data, "dgCMatrix")

    # Return log-normalised data
    return(rca.obj)
}
