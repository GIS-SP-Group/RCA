#' Create RCA object
#'
#' @param rawData A matrix of expression values
#' @param normData A matrix of normalized expression values
#' @return RCA object.
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#' @export
#'
#'
createRCAObject <- function(rawData, normData=NULL) {

    # Create RCA object using RCAConstruct and the raw data provided
    if (!(is.null(normData))){
    rca.obj <- methods::new(Class = "RCA", raw.data = rawData, data=normData)
    }else
    {
    rca.obj <- methods::new(Class = "RCA", raw.data = rawData)
    }
    # Return RCA object
    return(rca.obj)
}
