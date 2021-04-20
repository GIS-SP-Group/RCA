#' Filter the dataset and remove genes and cells that could be of poor quality using user-defined thresholds
#'
#' @param rca.obj RCA object.
#' @param nGene.thresholds numeric vector with lower and upper nGene thresholds
#' @param nUMI.thresholds numeric vector with lower and upper nUMI thresholds
#' @param percent.mito.thresholds numeric vector with lower and upper pMito thresholds
#' @param min.cell.exp minimum number of cells a gene must be expressed in
#' @param plot boolean - if True, plot data filter metrics
#' @param filename file name of saved plots
#' @return RCA object.
#' @export
#'
dataFilter <- function(rca.obj, nGene.thresholds = c(100, NULL),
                       nUMI.thresholds = c(1000, NULL), percent.mito.thresholds = c(0.0, 0.2),
                       min.cell.exp = 10, plot = T,  filename = "RCA_Filter.pdf") {

    if(!(("gridExtra" %in% base::.packages()) | !("ggplot2" %in% base::.packages()))){
	print("Ensure that grid.arrange and ggplot2 are installed")
	}
    # Extract data from RCA object
    data <- rca.obj$raw.data

    ### Gene filter ###

    if(!base::is.null(min.cell.exp)) {

        # Select genes with expression in a minimum number of cells
        geneExpVec <- Matrix::rowSums(data>0)
        filt.genes <- base::rownames(data)[base::which(geneExpVec > min.cell.exp)]

    } else {
        filt.genes <- base::rownames(data)
    }

    ### Cell filter - nGene ###
    # Assign filtered data
    filt.cells <- base::colnames(data)

    if(!base::is.null(nGene.thresholds)) {

        # Compute nGene vector
        nGeneVec <- Matrix::colSums(data>0)

        if(base::is.numeric(nGene.thresholds) & (base::length(nGene.thresholds) == 2)) {

            nGene.filt.cells <- filt.cells[base::which((nGeneVec >= nGene.thresholds[1]) & (nGeneVec <= nGene.thresholds[2]))]
        } else if(base::is.numeric(nGene.thresholds) & (base::length(nGene.thresholds) == 1)) {
            nGene.filt.cells <- filt.cells[base::which(nGeneVec >= nGene.thresholds)]
        } else {
            warning("nGene.thresholds was not of the appropriate format. Please enter a numeric vector with lower and upper thresholds.")
            nGene.filt.cells <- filt.cells
        }
    } else {
        warning("No nGene.thresholds provided.")
        nGene.filt.cells <- filt.cells
    }

    ### Cell filter - nUMI ###

    if(!is.null(nUMI.thresholds)) {

        # Compute nUMI vector
        nUMIVec <- Matrix::colSums(data)

        if(base::is.numeric(nUMI.thresholds) & (base::length(nUMI.thresholds) == 2)) {

            nUMI.filt.cells <- filt.cells[base::which((nUMIVec >= nUMI.thresholds[1]) & (nUMIVec <= nUMI.thresholds[2]))]
        } else if(base::is.numeric(nGene.thresholds) & (base::length(nUMI.thresholds) == 1)) {
            nUMI.filt.cells <- filt.cells[base::which(nUMIVec >= nUMI.thresholds)]
        } else {
            warning("nUMI.thresholds was not of the appropriate format. Please enter a numeric vector with lower and upper thresholds.")
            nUMI.filt.cells <- filt.cells
        }
    } else {
        warning("No nUMI.thresholds provided.")
        nUMI.filt.cells <- filt.cells
    }

    ### Cell filter - Mitochondrial Rate ###

    if(!base::is.null(percent.mito.thresholds)) {

        # Select mito genes
        mito.genes = base::grep(pattern = "^MT-|^Mt-", x = base::rownames(data), value = T)

        # Compute percent.mito vector
        pMitoVec <- Matrix::colSums(data[mito.genes, ])/Matrix::colSums(data)


        if(base::is.numeric(percent.mito.thresholds) & (base::length(percent.mito.thresholds) == 2)) {

            pMito.filt.cells <- filt.cells[base::which((pMitoVec >= percent.mito.thresholds[1]) & (pMitoVec <= percent.mito.thresholds[2]))]
        } else if(base::is.numeric(percent.mito.thresholds) & (base::length(percent.mito.thresholds) == 1)) {
            pMito.filt.cells <- filt.cells[, base::which(pMitoVec >= percent.mito.thresholds)]
        } else {
            warning("percent.mito.thresholds was not of the appropriate format. Please enter a numeric vector with lower and upper thresholds.")
            pMito.filt.cells <- filt.cells

        }
    } else {
        warning("No percent.mito.thresholds provided.")
        pMito.filt.cells <- filt.cells
    }

    # Select only the cells that satisfy all 3 cell filtering criteria
    filt.cells <- base::intersect(base::intersect(nGene.filt.cells, nUMI.filt.cells), pMito.filt.cells)

    # If plot is True
    if(plot) {

        # Create dataframe for 3D plot
        cellFiltDf <- base::data.frame(nGene = nGeneVec, nUMI = nUMIVec, pMito = pMitoVec)
        cellFiltDf$isFilt <- base::ifelse(test = base::rownames(cellFiltDf) %in% filt.cells, yes = "Retained", no = "Discarded")

        # Colors for scatter plot
        colors <- c("Discarded" = "grey", "Retained" = "black")

        # Create 2D plots for all pairwise combinations of nGene, nUMI and pMito

        # nGene x pMito

        # Create threshold lines
        nGene.thres.lines <- base::lapply(nGene.thresholds, function(thres){
            ggplot2::geom_vline(ggplot2::aes(xintercept = thres), color = "red", linetype="dashed", size=1)
        })

        pMito.thres.lines <- base::lapply(percent.mito.thresholds, function(thres){
            ggplot2::geom_hline(ggplot2::aes(yintercept = thres), color = "red", linetype="dashed", size=1)
        })

        # Create plot
        nGene <- pMito <- isFilt <- nUMI <- NULL
        nGene_pMito_plot <- ggplot2::ggplot(data = cellFiltDf, ggplot2::aes(x = nGene, y = pMito, colour = isFilt)) +
            ggplot2::geom_point(size = 1) + ggplot2::geom_jitter() + ggplot2::scale_color_manual(values = colors) +
            nGene.thres.lines + pMito.thres.lines + ggplot2::theme_bw() + ggplot2::ggtitle("a) nGene vs pMito")+
            ggplot2::theme(legend.position="none")

        # nUMI x pMito

        # Create threshold lines
        nUMI.thres.lines <- base::lapply(nUMI.thresholds, function(thres){
            ggplot2::geom_vline(ggplot2::aes(xintercept = thres), color = "red", linetype="dashed", size=1)
        })

        # Create plot
        nUMI_pMito_plot <- ggplot2::ggplot(data = cellFiltDf, ggplot2::aes(x = nUMI, y = pMito, colour = isFilt)) +
            ggplot2::geom_point(size = 1) + ggplot2::geom_jitter() +
            ggplot2::scale_color_manual(values = colors) + nUMI.thres.lines + pMito.thres.lines + ggplot2::theme_bw() +
            ggplot2::ggtitle("b) nUMI vs pMito")+ggplot2::theme(legend.position="none")

        # nGene x nUMI

        # Change nUMI threshold lines to horizontal lines
        nUMI.thres.lines <- base::lapply(nUMI.thresholds, function(thres){
            ggplot2::geom_hline(ggplot2::aes(yintercept = thres), color = "red", linetype="dashed", size=1)
        })

        # Create plot
        nGene_nUMI_plot <- ggplot2::ggplot(data = cellFiltDf, ggplot2::aes(x = nGene, y = nUMI, colour = isFilt)) +
            ggplot2::geom_point(size = 1) + ggplot2::geom_jitter() +
            ggplot2::scale_color_manual(values = colors) + nGene.thres.lines + nUMI.thres.lines +
            ggplot2::theme_bw() + ggplot2::ggtitle("c) nGene vs nUMI")+ ggplot2::theme(legend.position="right")+
            ggplot2::theme(legend.key.height=ggplot2::unit(1.5,"cm"))+ggplot2::labs(colour="")

        grDevices::pdf(file = filename, width = 15, height = 5)

        # Arrange scatter plots
        gridExtra::grid.arrange(nGene_pMito_plot, nUMI_pMito_plot, nGene_nUMI_plot, nrow = 1,widths=base::c(1,1,1.3))

        # Save plot
        grDevices::dev.off()
    }

    # Combine results from cell and gene filtering
    rca.obj$raw.data <- methods::as(data[filt.genes, filt.cells], "dgCMatrix")

    # Return
    return(rca.obj)
}
