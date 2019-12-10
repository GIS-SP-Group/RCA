#' Create RCA object
#'
#' @param rawData A matrix of expression values
#' @return RCA object.
#' @export
#'
createRCAObject <- function(rawData,dataIsNormalized=FALSE) {

    # Create RCA object using RCAConstruct and the raw data provided
    if (dataIsNormalized){
    rca.obj <- RCAConstruct$new(raw.data = rawData, data=rawData)
    }else
    {
    rca.obj <- RCAConstruct$new(raw.data = rawData)
    }
    # Return RCA object
    return(rca.obj)
}
