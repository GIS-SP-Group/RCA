#' Plot bar plots showing the composition of RCA clusters
#'
#' @param rca.obj data matrix (genes x cells)
#' @param folderpath path to save heatmap to
#' @param filename file name of saved heatmap
#' @export
#'

plotRCAHeatmap <- function(rca.obj, folderpath = ".", filename = "RCA_Heatmap.pdf") {

    ### Extract projection data and clustering result from RCA object
    heatmapIn = as.matrix(rca.obj$projection.data)
    cellTree = rca.obj$clustering.out$cellTree
    clusterColorList = rca.obj$clustering.out$dynamicColorsList

    ### Check if package dependencies are available; if not, download from CRAN and require those packages
    # dplyr Package
    if (!require(ComplexHeatmap))
        install.packages("ComplexHeatmap", repos = "http://cran.us.r-project.org")
    require(ComplexHeatmap)

    # ggplot2 Package
    if (!require(ComplexHeatmap))
        install.packages("ComplexHeatmap", repos = "http://cran.us.r-project.org")
    require(ComplexHeatmap)



}
