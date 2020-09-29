#' Filter the dataset and remove genes and cells that could be of poor quality using user-defined thresholds
#'
#' @param rca.obj RCA object.
#' @param nGene.thresholds numeric vector with lower and upper nGene thresholds
#' @param nUMI.thresholds numeric vector with lower and upper nUMI thresholds
#' @param percent.mito.thresholds numeric vector with lower and upper pMito thresholds
#' @param min.cell.exp minimum number of cells a gene must be expressed in
#' @param plot boolean - if True, plot data filter metrics
#' @param folderpath path to save data filter metric plots to
#' @param filename file name of saved plots
#' @return RCA object.
#' @export
#'
dataFilter <- function(rca.obj, nGene.thresholds = c(100, NULL), nUMI.thresholds = c(1000, NULL), percent.mito.thresholds = c(0.0, 0.2), min.cell.exp = 10, plot = T, folderpath = ".", filename = "RCA_Filter.pdf") {

    if(!(require(gridExtra) && require(ggplot2))){
	print("Ensure that grid.arrange and ggplot2 are installed")
	}
    # Extract data from RCA object
    data <- rca.obj$raw.data

    ### Gene filter ###

    if(!is.null(min.cell.exp)) {

        # Select genes with expression in a minimum number of cells
        geneExpVec <- Matrix::rowSums(data>0)
        filt.genes <- rownames(data)[which(geneExpVec > min.cell.exp)]

    } else {
        filt.genes <- rownames(data)
    }



    ### Cell filter - nGene ###

    # Assign filtered data
    filt.cells <- colnames(data)

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

    ### Cell filter - nUMI ###

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

    ### Cell filter - Mitochondrial Rate ###

    if(!is.null(percent.mito.thresholds)) {

        # Select mito genes
        mito.genes = grep(pattern = "^MT-", x = rownames(data), value = T)

        # Compute percent.mito vector
        pMitoVec <- Matrix::colSums(data[mito.genes, ])/Matrix::colSums(data)


        if(is.numeric(percent.mito.thresholds) & (length(percent.mito.thresholds) == 2)) {

            pMito.filt.cells <- filt.cells[which((pMitoVec >= percent.mito.thresholds[1]) & (pMitoVec <= percent.mito.thresholds[2]))]
        } else if(is.numeric(percent.mito.thresholds) & (length(percent.mito.thresholds) == 1)) {
            pMito.filt.cells <- filt.cells[, which(pMitoVec >= percent.mito.thresholds)]
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

    # If plot is True
    if(plot) {

        # Create dataframe for 3D plot
        cellFiltDf <- data.frame(nGene = nGeneVec, nUMI = nUMIVec, pMito = pMitoVec)
        cellFiltDf$isFilt <- ifelse(test = rownames(cellFiltDf) %in% filt.cells, yes = "Retained", no = "Discarded")

        # Colors for scatter plot
        colors <- c("Discarded" = "grey", "Retained" = "black")

        # Create 2D plots for all pairwise combinations of nGene, nUMI and pMito

        # nGene x pMito

        # Create threshold lines
        nGene.thres.lines <- lapply(nGene.thresholds, function(thres){
            geom_vline(aes(xintercept = thres), color = "red", linetype="dashed", size=1)
        })

        pMito.thres.lines <- lapply(percent.mito.thresholds, function(thres){
            geom_hline(aes(yintercept = thres), color = "red", linetype="dashed", size=1)
        })

        # Create plot
        nGene_pMito_plot <- ggplot(data = cellFiltDf, aes(x = nGene, y = pMito, colour = isFilt)) + geom_point(size = 1) + geom_jitter() + scale_color_manual(values = colors) + nGene.thres.lines + pMito.thres.lines + theme_bw() + ggtitle("a) nGene vs pMito")+theme(legend.position="none")

        # nUMI x pMito

        # Create threshold lines
        nUMI.thres.lines <- lapply(nUMI.thresholds, function(thres){
            geom_vline(aes(xintercept = thres), color = "red", linetype="dashed", size=1)
        })

        # Create plot
        nUMI_pMito_plot <- ggplot(data = cellFiltDf, aes(x = nUMI, y = pMito, colour = isFilt)) + geom_point(size = 1) + geom_jitter() + scale_color_manual(values = colors) + nUMI.thres.lines + pMito.thres.lines + theme_bw() + ggtitle("b) nUMI vs pMito")+theme(legend.position="none")

        # nGene x nUMI

        # Change nUMI threshold lines to horizontal lines
        nUMI.thres.lines <- lapply(nUMI.thresholds, function(thres){
            geom_hline(aes(yintercept = thres), color = "red", linetype="dashed", size=1)
        })

        # Create plot
        nGene_nUMI_plot <- ggplot(data = cellFiltDf, aes(x = nGene, y = nUMI, colour = isFilt)) + geom_point(size = 1) + geom_jitter() + scale_color_manual(values = colors) + nGene.thres.lines + nUMI.thres.lines + theme_bw() + ggtitle("c) nGene vs nUMI")+theme(legend.position="right")+theme(legend.key.height=unit(1.5,"cm"))+labs(colour="")

        pdf(file = paste0(folderpath,"/",filename), width = 15, height = 5)

        # Arrange scatter plots
        grid.arrange(nGene_pMito_plot, nUMI_pMito_plot, nGene_nUMI_plot, nrow = 1,widths=c(1,1,1.3))

        # Save plot
        dev.off()
    }

    # Combine results from cell and gene filtering
    rca.obj$raw.data <- as(data[filt.genes, filt.cells], "dgCMatrix")

    # Return
    return(rca.obj)
}
