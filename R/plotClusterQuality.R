#' Perform cluster specific QC
#'
#' @param rca.obj RCA object
#' @param cluster.labels vector of cluster labels
#' @param width width of plot in inches. Default is 20.
#' @param height height of plot in inches. Default is 20.
#' @param folderpath path to save heatmap to
#' @param filename file name of saved heatmap
#'
#' @export
#'

plotClusterQuality <- function(rca.obj, cluster.labels, width = 20, height = 20, folderpath = ".", filename = "RCA_Cluster_Quality.pdf") {

    # Create data frame for cell-cluster mapping
    cluster.df <- data.frame(Cell = colnames(rca.obj$data), Cluster = cluster.labels)

    # Create empty data frame for QC parameters
    quality.df <- data.frame(Cluster = character(), nGene = numeric(), nUMI = numeric(), pMito = numeric(), stringsAsFactors = FALSE)

    # For each cluster
    for(cluster in unique(cluster.labels)) {

        # Subset cluster data
        data <- rca.obj$raw.data[, subset(cluster.df$Cell, cluster.df$Cluster == cluster), drop = FALSE]

        # Compute nGene vector
        nGeneVec <- Matrix::colSums(data>0)

        # Compute nUMI vector
        nUMIVec <- Matrix::colSums(data)

        # Select mito genes
        mito.genes = grep(pattern = "^MT-", x = rownames(data), value = T)

        # Compute percent.mito vector
        pMitoVec <- Matrix::colSums(data[mito.genes, , drop = FALSE])/Matrix::colSums(data)

        # Append QC vectors to data frame
        cluster.quality.df <- data.frame(Cluster = rep(cluster, ncol(data)), nGene = nGeneVec, nUMI = nUMIVec, pMito = pMitoVec)
        quality.df <- rbind(quality.df, cluster.quality.df)
    }

    # nGene vs pMito plot
    pdf(file = paste0(folderpath, "/", "nGene_pMito_", filename), width = width, height = height)
    nGene_pMito_plot <- ggplot2::ggplot(data = quality.df, ggplot2::aes(x = nGene, y = pMito)) + ggplot2::geom_point(size = 1) + ggplot2::geom_jitter() + ggplot2::facet_wrap(.~Cluster, scales = "free",nrow=5) + ggplot2::theme_bw() + ggplot2::ggtitle("nGene vs pMito")+ggplot2::geom_density2d()
    print(nGene_pMito_plot)
    dev.off()


    # nUMI vs pMito plot
    pdf(file = paste0(folderpath, "/", "nUMI_pMito_", filename), width = width, height = height)
    nUMI_pMito_plot <- ggplot2::ggplot(data = quality.df, ggplot2::aes(x = nUMI, y = pMito)) + ggplot2::geom_point(size = 1) +ggplot2::geom_jitter() + ggplot2::facet_wrap(.~Cluster, scales = "free",nrow=5) + ggplot2::theme_bw() + ggplot2::ggtitle("nUMI vs pMito")+ ggplot2::geom_density2d()
    print(nUMI_pMito_plot)
    dev.off()

    # nGene vs nUMI plot
    pdf(file = paste0(folderpath, "/", "nGene_nUMI_", filename), width = width, height = height)
    nGene_nUMI_plot <- ggplot2::ggplot(data = quality.df, ggplot2::aes(x = nGene, y = nUMI)) + ggplot2::geom_point(size = 1) + ggplot2::geom_jitter() + ggplot2::facet_wrap(.~Cluster, scales = "free",nrow=5) + ggplot2::theme_bw() + ggplot2::ggtitle("nGene vs nUMI")+ ggplot2::geom_density2d()
    print(nGene_nUMI_plot)
    dev.off()

}
