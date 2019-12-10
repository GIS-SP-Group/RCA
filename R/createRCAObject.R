#' Create RCA object
#'
#' @param rawData A matrix of expression values
#' @return RCA object.
#' @export
#'
createRCAObject <- function(rawData) {

    # Create RCA object using RCAConstruct and the raw data provided
    rca.obj <- RCAConstruct$new(raw.data = rawData)

    # Return RCA object
    return(rca.obj)
}
