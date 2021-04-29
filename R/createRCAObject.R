#' Create RCA object
#'
#' @param rawData A matrix of expression values
#' @param normData A matrix of normalized expression values (default NULL)
#' @return RCA object.
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#' @examples
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' print(RCA.pbmcs)
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
