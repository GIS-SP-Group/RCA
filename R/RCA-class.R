#' RCA Class
#'
#' @export
#'
if (!require(Matrix))
    install.packages("Matrix", repos = "http://cran.us.r-project.org")
require(Matrix)
RCAConstruct <- setRefClass(Class = "RCA", fields = list(raw.data = "Matrix", data = "Matrix", projection.data = "Matrix", clustering.out = "list"))
