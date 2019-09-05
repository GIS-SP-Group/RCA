#' Log-transform and normalise data by sequencing depth
#'
#' @param rca.obj RCA object.
#' @return RCA object.
#' @export
#' @examples
#'
#' logNormData = dataLogNormalize(data);
#'
dataLogNormalise <- function(rca.obj) {
    
    # Extract data from RCA object
    data <- rca.obj$data
    
    # Compute sequencing depth vector
    seqDepthVec <- colSums(data)
    
    # Normalise data by cell
    norm.data <- sapply(seq_along(seqDepthVec), function(index) {
        data[,index]/seqDepthVec[index]
    })
    
    # Log-transform normalised data
    logNorm.data <- log(1+norm.data)
    
    # Store log-transformed normalised data in RCA object
    rca.obj$data <- logNorm.data
    
    # Return log-normalised data
    return(rca.obj)
}
