#' Generate cell clusters using hierarchical clustering and dynamic tree cutting.
#'
#' @param rca.obj RCA object.
#' @param k Number of cells to consider in the neirest neighbour graph (default 10).
#' @param eps Number of cells that have to be similar to connect 2 cells in a network (default 8).
#' @param minPts Number of cells with at least eps neighbours to be considered a core point (default 5).
#' @param dist.fun Either PCA or All. For PCA, a PCA reduction of the projection will be performed before the correlation matrix is computed (default All).
#' @param corMeth Correlation method used to compute distance (default pearson).
#' @return RCA object.
#'
#' @examples
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
#' RCA.pbmcs <- dataSNN(RCA.pbmcs)
#'
#' @export
#'
dataSNN <-
    function(rca.obj,
             k = 10,
             eps = 8,
             minPts = 5,
             dist.fun = "All",
             corMeth = "pearson") {
        projection.data <- base::as.matrix(rca.obj$projection.data)
        if (dist.fun == "PCA") {
            # Extract projection data
            base::print("Computing PCA on projection")
            pcaD = stats::prcomp(projection.data)
            components = base::c(1:(base::max(
                base::which(base::summary(pcaD)$importance[3, ] < 0.99)
            ) + 1))

            if (("randomcoloR" %in% base::.packages()) &
                corMeth == "pearson") {
                d = stats::as.dist(
                    1 - HiClimR::fastCor(
                        pcaD$rotation[, components],
                        upperTri = TRUE,
                        nSplit = 5,
                        optBLAS = T
                    )
                )
            } else {
                # else, use cor
                d = stats::as.dist(1 - stats::cor(pcaD$rotation[, components], method = corMeth))
            }

            # Obtain cell tree using a reduced projection matrix
            clusteringResult <-
                dbscan::sNNclust(d, k, eps, minPts, borderPoints = T)
        } else {
            base::print("Compute correlation distance matrix from projection")
            # If HiClimR is available, use fastCor to compute distance matrix
            if (("randomcoloR" %in% .packages())) {
                d = stats::as.dist(
                    1 - HiClimR::fastCor(
                        projection.data,
                        upperTri = TRUE,
                        nSplit = 5,
                        optBLAS = T
                    )
                )
            } else {
                # else, use cor
                d = stats::as.dist(1 - stats::cor(projection.data, method = corMeth))
            }

            # Obtain cell tree using a reduced projection matrix
            clusteringResult <-
                dbscan::sNNclust(d, k, eps, minPts, borderPoints = T)
        }

        # Convert labels to colours for each tree cut
        clusterColors <-
            randomcoloR::distinctColorPalette(base::length(base::unique(clusteringResult$cluster)))
        clusterColors <- base::sapply(clusterColors, plotrix::color.id)
        clusterColors <-
            base::sapply(clusterColors, function(x) {
                return(x[1])
            })
        base::names(clusterColors) <-
            base::unique(clusteringResult$cluster)

        dynamicColorsList <-
            base::list(Colors = clusterColors[base::as.character(clusteringResult$cluster)])
        base::names(dynamicColorsList) <- base::c("Clusters")
        # Assign clustering result to RCA object
        rca.obj$clustering.out <- base::list("cellTree" = clusteringResult$cluster,
                                             "dynamicColorsList" = dynamicColorsList)

        # Return RCA object
        return(rca.obj)
    }
