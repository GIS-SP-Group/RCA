#' RCA Class
#'
#' @export
#'
RCAConstruct <- setRefClass(Class = "RCA", fields = list(raw.data = "Matrix", data = "Matrix", projection.data = "Matrix", clustering.out = "list"))
