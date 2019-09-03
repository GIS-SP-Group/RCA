#' Log-transform and normalise data by sequencing depth
#'
#' @param data data matrix (genes x cells)
#' @return log-normalised data matrix
#' @export
#' @examples
#'
#' logNormData = dataLogNormalize(data);
#'
dataLogNormalise <- function(data) {
    
    # Compute sequencing depth vector
    seqDepthVec <- colSums(data)
    
    # Normalise data by cell
    norm.data <- sapply(seq_along(seqDepthVec), function(index) {
        data[,index]/seqDepthVec[index]
    })
    
    # Log-transform normalised data
    logNorm.data <- log(1+norm.data)
    
    # Return log-normalised data
    return(logNorm.data)
}