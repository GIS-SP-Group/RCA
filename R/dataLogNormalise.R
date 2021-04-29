#' Log-transform and normalise data by sequencing depth
#'
#' @param rca.obj RCA object.
#' @param scale.factor scaling factor for log-normalization (default 10,000).
#' @return RCA object.
#' @examples
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' print(RCA.pbmcs)
#'
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
    logNorm.data <- base::log(1+norm.data)

    base::colnames(logNorm.data) <- base::colnames(raw.data)

    # Store log-transformed normalised data in RCA object
    rca.obj$data <- methods::as(logNorm.data, "dgCMatrix")

    # Return log-normalised data
    return(rca.obj)
}
