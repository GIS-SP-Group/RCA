#' Filter the dataset and remove genes and cells that could be of poor quality using user-defined thresholds
#'
#' @param rca.obj RCA object.
#' @param cluster.labels vector of cluster labels
#' @param nGene.low.thresholds numeric vector with lower nGene thresholds. Length of vector should be length(unique(cluster.labels)).
#' @param nGene.high.thresholds numeric vector with upper nGene thresholds. Length of vector should be length(unique(cluster.labels)).
#' @param nUMI.low.thresholds numeric vector with lower nUMI thresholds. Length of vector should be length(unique(cluster.labels)).
#' @param nUMI.high.thresholds numeric vector withupper nUMI thresholds. Length of vector should be length(unique(cluster.labels)).
#' @param pMito.low.thresholds numeric vector with lower pMito thresholds. Length of vector should be length(unique(cluster.labels)).
#' @param pMito.high.thresholds numeric vector and upper pMito thresholds. Length of vector should be length(unique(cluster.labels)).
#'
#' @return RCA object.
#' @export
#'
performClusterSpecificQC <- function(rca.obj, cluster.labels, nGene.low.thresholds = NULL, nGene.high.thresholds = NULL, nUMI.low.thresholds = NULL, nUMI.high.thresholds = NULL, pMito.low.thresholds = NULL, pMito.high.thresholds = NULL) {

    # Create data frame for cell-cluster mapping
    cluster.df <- data.frame(Cell = colnames(rca.obj$data), Cluster = cluster.labels)

    # Create dataframe for QC parameters
    qc.df <- data.frame(row.names = unique(cluster.labels), nGene.low = if(is.null(nGene.low.thresholds)) rep(0, length(unique(cluster.labels))) else nGene.low.thresholds, nGene.high = if(is.null(nGene.high.thresholds)) rep(Inf, length(unique(cluster.labels))) else nGene.high.thresholds, nUMI.low = if(is.null(nUMI.low.thresholds)) rep(0, length(unique(cluster.labels))) else nUMI.low.thresholds, nUMI.high = if(is.null(nUMI.high.thresholds)) rep(Inf, length(unique(cluster.labels))) else nUMI.high.thresholds, pMito.low = if(is.null(pMito.low.thresholds)) rep(0, length(unique(cluster.labels))) else pMito.low.thresholds, pMito.high = if(is.null(pMito.high.thresholds)) rep(1, length(unique(cluster.labels))) else pMito.high.thresholds)

    # Create empty filtered cell-cluster mapping data frame
    filt.cluster.df <- data.frame(Cell = character(), Cluster = character())

    # For each cluster
    for(cluster in unique(cluster.labels)) {

        # Subset cluster data
        data <- rca.obj$raw.data[, subset(cluster.df$Cell, cluster.df$Cluster == cluster)]

        # Assign filtered data
        filt.cells <- colnames(data)

        # Get nGene thresholds
        nGene.thresholds = c(qc.df[cluster, ]$nGene.low, qc.df[cluster, ]$nGene.high)

        # Filter nGene
        if(!is.null(nGene.thresholds)) {

            # Compute nGene vector
            nGeneVec <- Matrix::colSums(data>0)

            if(is.numeric(nGene.thresholds) & (length(nGene.thresholds) == 2)) {

                nGene.filt.cells <- filt.cells[which((nGeneVec >= nGene.thresholds[1]) & (nGeneVec <= nGene.thresholds[2]))]
            } else if(is.numeric(nGene.thresholds) & (length(nGene.thresholds) == 1)) {
                nGene.filt.cells <- filt.cells[which(nGeneVec >= nGene.thresholds)]
            } else {
                warning("nGene.thresholds was not of the appropriate format. Please enter a numeric vector with lower and upper thresholds.")
                nGene.filt.cells <- filt.cells
            }
        } else {
            warning("No nGene.thresholds provided.")
            nGene.filt.cells <- filt.cells
        }


        # Get nUMI thresholds
        nUMI.thresholds = c(qc.df[cluster, ]$nUMI.low, qc.df[cluster, ]$nUMI.high)

        # Filter by nUMI
        if(!is.null(nUMI.thresholds)) {

            # Compute nUMI vector
            nUMIVec <- Matrix::colSums(data)

            if(is.numeric(nUMI.thresholds) & (length(nUMI.thresholds) == 2)) {

                nUMI.filt.cells <- filt.cells[which((nUMIVec >= nUMI.thresholds[1]) & (nUMIVec <= nUMI.thresholds[2]))]
            } else if(is.numeric(nGene.thresholds) & (length(nUMI.thresholds) == 1)) {
                nUMI.filt.cells <- filt.cells[which(nUMIVec >= nUMI.thresholds)]
            } else {
                warning("nUMI.thresholds was not of the appropriate format. Please enter a numeric vector with lower and upper thresholds.")
                nUMI.filt.cells <- filt.cells
            }
        } else {
            warning("No nUMI.thresholds provided.")
            nUMI.filt.cells <- filt.cells
        }

        # Get pMito thresholds
        pMito.thresholds = c(qc.df[cluster, ]$pMito.low, qc.df[cluster, ]$pMito.high)

        # Filter by pMito
        if(!is.null(pMito.thresholds)) {

            # Select mito genes
            mito.genes = grep(pattern = "^MT-", x = rownames(data), value = T)

            # Compute percent.mito vector
            pMitoVec <- Matrix::colSums(data[mito.genes, ])/Matrix::colSums(data)


            if(is.numeric(pMito.thresholds) & (length(pMito.thresholds) == 2)) {

                pMito.filt.cells <- filt.cells[which((pMitoVec >= pMito.thresholds[1]) & (pMitoVec <= pMito.thresholds[2]))]
            } else if(is.numeric(pMito.thresholds) & (length(pMito.thresholds) == 1)) {
                pMito.filt.cells <- filt.cells[, which(pMitoVec >= pMito.thresholds)]
            } else {
                warning("percent.mito.thresholds was not of the appropriate format. Please enter a numeric vector with lower and upper thresholds.")
                pMito.filt.cells <- filt.cells

            }
        } else {
            warning("No percent.mito.thresholds provided.")
            pMito.filt.cells <- filt.cells
        }

        # Select only the cells that satisfy all 3 cell filtering criteria
        filt.cells <- intersect(intersect(nGene.filt.cells, nUMI.filt.cells), pMito.filt.cells)
        filt.cluster.df <- rbind(filt.cluster.df, subset(cluster.df, Cell %in% filt.cells))

    }

    # Subset data in RCA object
    rca.obj$raw.data <- rca.obj$raw.data[, filt.cluster.df$Cell]
    rca.obj$data <- rca.obj$data[, filt.cluster.df$Cell]

    # Subset cluster labels
    for(i in 1:length(rca.obj$clustering.out$dynamicColorsList)) {
        rca.obj$clustering.out$dynamicColorsList[[i]] <- rca.obj$clustering.out$dynamicColorsList[[i]][which(cluster.df$Cell %in% filt.cluster.df$Cell)]
    }

    # Return
    return(rca.obj)

}
