#' Generate cell clusters using hierarchical clustering and dynamic tree cutting.
#'
#' @param rca.obj RCA object.
#' @param k Number of cells to consider in the neirest neighbour graph. Default is 10.
#' @param eps Number of cells that have to be similar to connect 2 cells in a network. Default is 8
#' @param minPts Number of cells with at least eps neighbours to be considered a core point. Default is 5.
#' @return RCA object.
#' @export
#'
dataSNN <- function(rca.obj,k=10,eps=8,minPts=5) {

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

    # plotrix
    if (!require(plotrix)){
	install.packages("plotrix")
    }
    require(plotrix)


    pcaD = prcomp(projection.data)
    components=c(1:(max(which(summary(pcaD)$importance[3,]<0.99))+1))
    # Obtain cell tree using a reduced projection matrix.
    clusteringResult<-sNNclust(pcaD$rotation[,components],k,eps,minPts,borderPoints = T)

    # Convert labels to colours for each tree cut
    clusterColors<-distinctColorPalette(length(unique(clusteringResult$cluster)))
    clusterColors<-sapply(clusterColors,color.id)
    clusterColors<-sapply(clusterColors,function(x){return(x[1])})
    names(clusterColors)<-unique(clusteringResult$cluster)

    dynamicColorsList<-list(Colors=clusterColors[as.character(clusteringResult$cluster)])

    # Assign clustering result to RCA object
    rca.obj$clustering.out <- list(
        "cellTree" = clusteringResult$cluster,
        "dynamicColorsList" = dynamicColorsList
    )

    ### Return RCA object

    return(rca.obj)
}
