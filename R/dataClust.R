#' Generate cell clusters using hierarchical clustering and dynamic tree cutting.
#'
#' @param projection_data result of projection onto RCA panel
#' @param deepSplitValues integer vector indicating how deep the dendrogram should be cut. Values can range from 0 to 4.
#' @param minClustSize integer value indicating the minimum size of the resulting clusters. Default is 5.
#' @return list of clustering results.
#' @export
#' @examples
#'
#' clustering_output = cellClust(projection_data, 1:3, 10);
#'
dataClust <- function(projection_data, deepSplitValues = 1, minClustSize = 5) {
        
    ### Load packages
    
    # fastcluster
    if (!require(fastcluster))
        install.packages("fastcluster", repos = "http://cran.us.r-project.org")
    require(fastcluster)
    
    # WGCNA
    if (!require(WGCNA)) {
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
            projection_data,
            upperTri = TRUE,
            nSplit = 5,
            optBLAS = T
        ))
    } else {
        # else, use cor
        d = as.dist(1 - cor(projection_data, method = "pearson"))
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
    
    ### Return list of clustering results
    
    return(list(
        "d" = d,
        "cellTree" = cellTree,
        "dynamicColorsList" = dynamicColorsList
    ))
}
