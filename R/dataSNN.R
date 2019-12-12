#' Generate cell clusters using hierarchical clustering and dynamic tree cutting.
#'
#' @param rca.obj RCA object.
#' @param deepSplitValues integer vector indicating how deep the dendrogram should be cut. Values can range from 0 to 4.
#' @param minClustSize integer value indicating the minimum size of the resulting clusters. Default is 5.
#' @return RCA object.
#' @export
#'
dataSNN <- function(rca.obj,k=50,eps=20,minPts=10) {

    ### Extract projection data
    projection.data <- as.matrix(rca.obj$projection.data)
    ### Load packages

    # fastcluster
    if (!require(dbscan))
        install.packages("dbscan", repos = "http://cran.us.r-project.org")
    require(dbscan)

    # randomColorR
    if (!require(randomcoloR)) {
        install.packages("randomcoloR")
    }
    require(randomcoloR)

    # HiClimR
    if (!require(HiClimR)) {
        install.packages("HiClimR")
    }
    require(HiClimR)

    ### Cluster cells
    # If HiClimR is available, use fastCor to compute distance matrix
    if (require(HiClimR)) {
        d = as.dist(1 - fastCor(
            projection.data,
            upperTri = TRUE,
            nSplit = 5,
            optBLAS = T
        ))
    } else {
        # else, use cor
        d = as.dist(1 - cor(projection.data, method = "pearson"))
    }

    # Obtain cell tree using distance matrix
    clusteringResult<-sNNclust(d,k,eps,minPts,borderPoints = T)

    # Convert labels to colours for each tree cut

    clusterColors<-distinctColorPalette(length(unique(clusteringResult$cluster)))
    names(clusterColors)<-unique(clusteringResult$cluster)

    dynamicColorsList<-list(Colors=clusterColors[as.character(clusteringResult$cluster)])

    # Assign clustering result to RCA object
    rca.obj$clustering.out <- list(
        "d" = d,
        "cellTree" = clusteringResult$cluster,
        "dynamicColorsList" = dynamicColorsList
    )

    ### Return RCA object

    return(rca.obj)
}
