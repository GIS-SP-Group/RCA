#' Create Reference Class
#' 
#' @param data A matrix of expression values 
#' @return RCA object.
#' @export
#' @examples
#' 
#' rca.obj = dataConstruct(data);
#' 
createRCAObject <- function(data) {
    
    RCAConstruct <- setRefClass(Class = "RCA", fields = list(raw.data = "matrix", data = "matrix", projection.data = "matrix", clustering.out = "list"))
    
    rca.obj <- RCAConstruct(raw.data = data)
    
    return(rca.obj)
}