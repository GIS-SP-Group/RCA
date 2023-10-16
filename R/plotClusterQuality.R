#' Perform cluster specific QC
#'
#' @param rca.obj RCA object
#' @param cluster.labels vector of cluster labels for each cell
#' @param width width of plot in inches (default 20).
#' @param height height of plot in inches (default 20).
#' @param filename file name of saved scatter plots
#'
#' @examples
#' \dontrun{
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
#' RCA.pbmcs <- dataClust(RCA.pbmcs)
#' plotClusterQuality(RCA.pbmcs,
#' RCA.pbmcs$clustering.out$dynamicColorsList[[1]])
#' }
#' @export
#'

plotClusterQuality <- function(rca.obj, cluster.labels,
                               width = 20, height = 20,
                               filename = "RCA_Cluster_Quality.pdf") {
    nGene <- pMito <- nUMI <- NULL
    # Create data frame for cell-cluster mapping
    cluster.df <- base::data.frame(Cell = base::colnames(rca.obj$data), Cluster = cluster.labels)

    # Create empty data frame for QC parameters
    quality.df <- base::data.frame(Cluster = base::character(),
                                   nGene = base::numeric(),
                                   nUMI = base::numeric(),
                                   pMito = base::numeric(),
                                   stringsAsFactors = FALSE)

    # For each cluster
    for(cluster in base::unique(cluster.labels)){

        # Subset cluster data
        data <- rca.obj$raw.data[, base::subset(cluster.df$Cell, cluster.df$Cluster == cluster), drop = FALSE]

        # Compute nGene vector
        nGeneVec <- Matrix::colSums(data > 0)

        # Compute nUMI vector
        nUMIVec <- Matrix::colSums(data)

        # Select mito genes
        mito.genes = base::grep(pattern = "^MT-", x = base::rownames(data), value = T)

        # Compute percent.mito vector
        pMitoVec <- Matrix::colSums(data[mito.genes, , drop = FALSE])/Matrix::colSums(data)*100

        # Append QC vectors to data frame
        cluster.quality.df <- base::data.frame(Cluster = base::rep(cluster, base::ncol(data)), nGene = nGeneVec, nUMI = nUMIVec, pMito = pMitoVec)
        quality.df <- base::rbind(quality.df, cluster.quality.df)
    }

    # nGene vs pMito plot
    grDevices::pdf(file = base::paste0("nGene_pMito_", filename), width = width, height = height)
    nGene_pMito_plot <- ggplot2::ggplot(data = quality.df, ggplot2::aes(x = nGene, y = pMito)) +
                        ggplot2::geom_point(size = 1) +
                        ggplot2::geom_jitter() +
                        ggplot2::facet_wrap(.~Cluster, scales = "free",ncol = 5) +
                        ggplot2::theme_bw() +
                        ggplot2::ggtitle("nGene vs pMito") +
                        ggplot2::geom_density2d()
    base::print(nGene_pMito_plot)
    grDevices::dev.off()


    # nUMI vs pMito plot
    grDevices::pdf(file = base::paste0("nUMI_pMito_", filename), width = width, height = height)
    nUMI_pMito_plot <- ggplot2::ggplot(data = quality.df, ggplot2::aes(x = nUMI, y = pMito)) +
                       ggplot2::geom_point(size = 1) +
                       ggplot2::geom_jitter() +
                       ggplot2::facet_wrap(.~Cluster, scales = "free",ncol = 5) +
                       ggplot2::theme_bw() +
                       ggplot2::ggtitle("nUMI vs pMito") +
                       ggplot2::geom_density2d()
    base::print(nUMI_pMito_plot)
    grDevices::dev.off()

    # nGene vs nUMI plot
    grDevices::pdf(file = base::paste0("nGene_nUMI_", filename), width = width, height = height)
    nGene_nUMI_plot <- ggplot2::ggplot(data = quality.df,ggplot2::aes(x = nGene, y = nUMI)) +
                       ggplot2::geom_point(size = 1) +
                       ggplot2::geom_jitter() +
                       ggplot2::facet_wrap(.~Cluster, scales = "free",ncol = 5) +
                       ggplot2::theme_bw() +
                       ggplot2::ggtitle("nGene vs nUMI") +
                       ggplot2::geom_density2d()
    base::print(nGene_nUMI_plot)
    grDevices::dev.off()

}
