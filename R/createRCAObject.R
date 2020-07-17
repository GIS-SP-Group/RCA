#' Create RCA object
#'
#' @param rawData A matrix of expression values
#' @param normData A matrix of normalized expression values
#' @return RCA object.
#' @export
#'
createRCAObject <- function(rawData, normData=NULL) {

    # Create RCA object using RCAConstruct and the raw data provided
    if (!(is.null(normData))){
    rca.obj <- RCAConstruct$new(raw.data = rawData, data=normData)
    }else
    {
    rca.obj <- RCAConstruct$new(raw.data = rawData)
    }
    # Return RCA object
    return(rca.obj)
}
