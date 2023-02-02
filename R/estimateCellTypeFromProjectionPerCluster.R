#' Estimate the most likely cell type from the projection to the reference panel
#'
#' @param rca.obj RCA object.
#' @param homogeneity a parameter indicating the homogeneity of the cluster. If the difference is below this threshold, the cell type will be set to unknown (default NULL).
#' @return RCA object.
#'
#' @examples
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
#' RCA.pbmcs <- dataClust(RCA.pbmcs)
#' RCA.pbmcs <- estimateCellTypeFromProjectionPerCluster(RCA.pbmcs)
#' @export
#'

estimateCellTypeFromProjectionPerCluster <- function(rca.obj, homogeneity=NULL) {

    # Assign data to tmp variables
    projection <- rca.obj$projection.data
    clusterColors <- rca.obj$clustering.out$dynamicColorsList[[1]]
    
    cellTypes <- as.list(rownames(projection)[max.col(Matrix::t(projection))])
  
    # Compute the cell type compositions of each cluster using per cell cell type predictions
    enrichmentAll <- base::c()
    Count <- Cluster <- NULL
    for (type in base::unique(clusterColors)) {
	index = base::which(clusterColors == type)
        enrichmentAll <- rbind(enrichmentAll,(base::cbind(type,base::table(base::unlist(cellTypes)[index]))))
    }
    enrichmentAll <- base::data.frame(base::cbind(base::row.names(enrichmentAll), enrichmentAll))
    base::colnames(enrichmentAll) <- base::c("CT", "Cluster", "Count")
    base::rownames(enrichmentAll) <- base::c(1:base::dim(enrichmentAll)[1])
    enrichmentAll$Count <- base::as.numeric(base::as.character(enrichmentAll$Count))
    totalCounts <- base::data.frame(dplyr::count(enrichmentAll, wt = Count, Cluster))
    enrichmentAll <- dplyr::left_join(enrichmentAll, totalCounts, by = "Cluster")
    enrichmentAll <- cbind(enrichmentAll, enrichmentAll$Count / enrichmentAll$n * 100)
    base::colnames(enrichmentAll)[5] <- "Ratio"
    enrichmentAll$Ratio <- base::as.numeric(base::as.character(enrichmentAll$Ratio))

    # Determine the identity of a cluster with or without the option to set the cluster identity to unkown
    maxCellTypeCluster <- base::list()
    homologyScores <- base::list()
    for (type in base::unique(clusterColors)) {
        subset <- enrichmentAll[base::which(enrichmentAll$Cluster == type), ]
        maxIndex <- base::which(subset$Ratio == base::max(subset$Ratio))
        homologyScores <-
            base::c(homologyScores, subset$Ratio[maxIndex])
        if (base::is.null(homogeneity)) {
            maxCellTypeCluster <-
                base::c(maxCellTypeCluster, base::as.character(subset$CT[maxIndex]))
        } else{
            if (subset$Ratio[maxIndex] > 0.5) {
                maxCellTypeCluster <-
                    base::c(maxCellTypeCluster,
                            base::as.character(subset$CT[maxIndex]))
            } else{
                maxCellTypeCluster <- base::c(maxCellTypeCluster, "Unkown")
            }
        }
    }
    base::names(homologyScores) <- base::unique(clusterColors)
    base::names(maxCellTypeCluster) <- base::unique(clusterColors)

    clusterConfidence <- homologyScores[clusterColors]
    estimatedCellTypes <- maxCellTypeCluster[clusterColors]

    rca.obj$cell.Type.Estimate.per.cluster <- estimatedCellTypes
    rca.obj$cScore <- clusterConfidence


    # Return RCA object
    return(rca.obj)
}


