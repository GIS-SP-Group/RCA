#' Generate cell clusters using hierarchical clustering and dynamic tree cutting.
#'
#' @param rca.obj RCA object.
#' @param k Number of cells to consider in the neirest neighbour graph. Default is 10.
#' @param eps Number of cells that have to be similar to connect 2 cells in a network. Default is 8
#' @param minPts Number of cells with at least eps neighbours to be considered a core point. Default is 5.
#' @param dist.fun Either PCA or All. For PCA, a PCA reduction of the projection will be performed before the correlation matrix is computed.
#' @return RCA object.
#' @export
#'
dataSNN <- function(rca.obj,k=10,eps=8,minPts=5,dist.fun="All") {
    projection.data <- as.matrix(rca.obj$projection.data)
    if (dist.fun=="PCA"){
    # Extract projection data 
    print("Computing PCA on projection")
    pcaD = stats::prcomp(projection.data)
    components=c(1:(max(which(summary(pcaD)$importance[3,]<0.99))+1))

    if (require(HiClimR)) {
        d = as.dist(1 - HiClimR::fastCor(
            pcaD$rotation[,components],,
            upperTri = TRUE,
            nSplit = 5,
            optBLAS = T
        ))
    } else {
        # else, use cor
        d = as.dist(1 - cor(pcaD$rotation[,components], method = "pearson"))
    }
 
    # Obtain cell tree using a reduced projection matrix
    clusteringResult<-dbscan::sNNclust(d,k,eps,minPts,borderPoints = T)
    } else {
    print("Compute correlation distance matrix from projection")
    # If HiClimR is available, use fastCor to compute distance matrix
    if (require(HiClimR)) {
        d = as.dist(1 - HiClimR::fastCor(
            projection.data,
            upperTri = TRUE,
            nSplit = 5,
            optBLAS = T
        ))
    } else {
        # else, use cor
        d = as.dist(1 - cor(projection.data, method = "pearson"))
    }
    
    # Obtain cell tree using a reduced projection matrix
    clusteringResult<-dbscan::sNNclust(d,k,eps,minPts,borderPoints = T)
}
    
    # Convert labels to colours for each tree cut
    clusterColors<-randomcoloR::distinctColorPalette(length(unique(clusteringResult$cluster)))
    clusterColors<-sapply(clusterColors,plotrix::color.id)
    clusterColors<-sapply(clusterColors,function(x){return(x[1])})
    names(clusterColors)<-unique(clusteringResult$cluster)

    dynamicColorsList<-list(Colors=clusterColors[as.character(clusteringResult$cluster)])

    # Assign clustering result to RCA object
    rca.obj$clustering.out <- list(
        "cellTree" = clusteringResult$cluster,
        "dynamicColorsList" = dynamicColorsList
    )

    # Return RCA object
    return(rca.obj)
}
