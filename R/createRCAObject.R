#' Create Reference Class
#'
#' @param raw_data A matrix of expression values
#' @return RCA object.
#' @export
#'
createRCAObject <- function(raw_data) {

    # Create RCA object using RCAConstruct and the raw data provided
    rca.obj <- RCAConstruct$new(raw.data = raw_data, data = raw_data)

    # Return RCA object
    return(rca.obj)
}
