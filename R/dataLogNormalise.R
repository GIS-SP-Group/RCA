#' Log-transform and normalise data by sequencing depth
#'
#' @param rca.obj RCA object.
#' @param scale.factor scaling factor for log-normalization
#' @return RCA object.
#' @export
#'
dataLogNormalise <- function(rca.obj, scale.factor = 10000) {

    # Extract data from RCA object
    raw.data <- rca.obj$raw.data

    # Compute sequencing depth vector
    seqDepthVec <- Matrix::colSums(raw.data)

    # Normalise data by cell
    norm.data<-raw.data/seqDepthVec*scale.factor

    # Log-transform normalised data
    logNorm.data <- log(1+norm.data)

    colnames(logNorm.data) <- colnames(raw.data)

    # Store log-transformed normalised data in RCA object
    rca.obj$data <- as(logNorm.data, "dgCMatrix")

    # Return log-normalised data
    return(rca.obj)
}
