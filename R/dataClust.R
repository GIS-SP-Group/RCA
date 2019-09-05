#' Generate cell clusters using hierarchical clustering and dynamic tree cutting.
#'
#' @param rca.obj RCA object.
#' @param deepSplitValues integer vector indicating how deep the dendrogram should be cut. Values can range from 0 to 4.
#' @param minClustSize integer value indicating the minimum size of the resulting clusters. Default is 5.
#' @return RCA object.
#' @export
#'
dataClust <- function(rca.obj, deepSplitValues = 1, minClustSize = 5) {

    ### Extract projection data
    projection.data <- rca.obj$projection.data
    ### Load packages

    # fastcluster
    if (!require(fastcluster))
        install.packages("fastcluster", repos = "http://cran.us.r-project.org")
    require(fastcluster)

    # WGCNA
    if (!require(WGCNA)) {
        if(!require(BiocManager))
            install.packages("BiocManager")
        BiocManager::install(c("impute", "GO.db", "preprocessCore"))
        install.packages("WGCNA")
    }
    require(WGCNA)

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
    cellTree = fastcluster::hclust(d, method = "average")

    # For each deepsplit value given, compute dynamic groups
    dynamicGroupsList <-
        lapply(X = deepSplitValues, function(deepSplit) {
            cutreeDynamic(
                dendro = cellTree,
                distM = as.matrix(d),
                deepSplit = deepSplit,
                pamStage = FALSE,
                minClusterSize = minClustSize
            )
        })

    # Convert labels to colours for each tree cut
    dynamicColorsList <- lapply(dynamicGroupsList, labels2colors)
    names(dynamicColorsList) <- paste0("deepSplit ", deepSplitValues)

    # Assign clustering result to RCA object
    rca.obj$clustering.out <- list(
        "d" = d,
        "cellTree" = cellTree,
        "dynamicColorsList" = dynamicColorsList
    )

    ### Return RCA object

    return(rca.obj)
}
