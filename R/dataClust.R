#' Generate cell clusters using hierarchical clustering and dynamic tree cutting.
#'
#' @param rca.obj RCA object.
#' @param deepSplitValues integer vector indicating how deep the dendrogram should be cut. Values can range from 0 to 4.
#' @param minClustSize integer value indicating the minimum size of the resulting clusters. Default is 5.
#' @param corMeth Correlation method used to compute the distance matrix of the projection (pearson (default), spearman, kendal).
#' @return RCA object.
#' @export
#'
dataClust <- function(rca.obj, deepSplitValues = 1, minClustSize = 5, corMeth = "pearson") {

    ### Extract projection data
    projection.data <- base::as.matrix(rca.obj$projection.data)

    ### Cluster cells

    # If HiClimR is available, use fastCor to compute distance matrix
    if (("HiClimR" %in% .packages()) & (corMeth=="pearson")) {
        d = stats::as.dist(1 - HiClimR::fastCor(
            projection.data,
            upperTri = TRUE,
            nSplit = 5,
            optBLAS = T
        ))
    } else {
        # else, use cor
        d = stats::as.dist(1 - stats::cor(projection.data, method = corMeth))
    }

    # Obtain cell tree using distance matrix
    if ("fastcluster" %in% base::.packages()){
        cellTree = fastcluster::hclust(d, method = "average")
    } else {
        cellTree = stats::hclust(d, method = "average")
    }

    # For each deepsplit value given, compute dynamic groups
    dynamicGroupsList <-
        base::lapply(X = deepSplitValues, function(deepSplit) {
            dynamicTreeCut::cutreeDynamic(
                dendro = cellTree,
                distM = base::as.matrix(d),
                deepSplit = deepSplit,
                pamStage = FALSE,
                minClusterSize = minClustSize
            )
        })

    # Convert labels to colours for each tree cut
    dynamicColorsList <- base::lapply(dynamicGroupsList, WGCNA::labels2colors)
    base::names(dynamicColorsList) <- base::paste0("deepSplit ", deepSplitValues)

    # Assign clustering result to RCA object
    rca.obj$clustering.out <- base::list(
        "cellTree" = cellTree,
        "dynamicColorsList" = dynamicColorsList
    )

    ### Return RCA object
    return(rca.obj)
}
