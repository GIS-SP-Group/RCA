#' RCA Class
#'
#' @export
#'
RCAConstruct <- setRefClass(Class = "RCA", fields = list(raw.data = "matrix", data = "matrix", projection.data = "matrix", clustering.out = "list"))
