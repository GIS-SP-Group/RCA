#' Plot umap of projection to the RCA panel
#'
#' @param rca.obj RCA object
#' @param cellPropertyList list of cell properties to plot
#' @param folderpath path to save umap to
#' @param filename file name of saved umap
#' @export
#'

plotRCAUMAP <- function(rca.obj, cellPropertyList = NULL, folderpath = ".", filename = "RCA_UMAP.pdf") {

    ### Extract projection data from RCA object
    projection = as.matrix(rca.obj$projection.data)
    clusterColorList = rca.obj$clustering.out$dynamicColorsList

    ### Check if package dependencies are available; if not, download from CRAN and require those packages
    # umap
    if (!require(umap))
        install.packages("umap", repos = "http://cran.us.r-project.org")
    require(umap)
    if (!require(ggplot2))
        install.packages("ggplot2", repos = "http://cran.us.r-project.org")
    require(ggplot2)
    if (!require(ggpubr))
        install.packages("ggpubr", repos = "http://cran.us.r-project.org")
    require(ggpubr)

    # Compute UMAP projection from cell type projection
    umap.projection <- umap(t(projection))

    # Store UMAP layout in data frame for plotting
    umap.df <- as.data.frame(umap.projection$layout)
    colnames(umap.df) <- c("UMAP1", "UMAP2")

    umapPlots<-list()
    # If no cluster colors or cell properties are to be plotted
    if(is.null(clusterColorList) & is.null(cellPropertyList)) {

        # Plot UMAP of cells without annotations
        umap.plot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2)) + geom_point(size = .5) + theme_classic()
	umapPlots<-c(umapPlots,list(umap.plot))
        # Save UMAP
        ggsave(filename = paste0(folderpath, "/", filename), plot = umap.plot)

    } else {

        # If cluster colors are to be plotted
        if(!is.null(clusterColorList)) {

            # Create a list of UMAP plots for each cluster coloring
            for(index in seq_along(clusterColorList)) {

                # Get the name of this cluster annotation
                clusterColorName = names(clusterColorList[index])

                # Set the data frame column to the color vector
                umap.df[[clusterColorName]] <- clusterColorList[[index]]

                # Create the plot
                umapClusterColorsPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, colour = umap.df[[clusterColorName]])) + geom_point(size = .5) + scale_color_manual(values = sort(unique(umap.df[[clusterColorName]]))) + labs(colour = clusterColorName) + theme_bw()

                # Save plot
                ggsave(filename = paste0(folderpath, "/", "ClusterColors_", clusterColorName,"_", filename), plot = umapClusterColorsPlot,width=10,height=8,units="in")
		umapPlots<-c(umapPlots,list(umapClusterColorsPlot))

            }

        }

        # If cell properties are to be plotted
        if(!is.null(cellPropertyList)) {

            # Create a list of UMAP plots for each cluster coloring
            for(index in seq_along(cellPropertyList)) {

                # Get the name of this cluster annotation
                CellPropertyName = names(cellPropertyList)[index]

                # Set the data frame column to the color vector
                umap.df[[CellPropertyName]] <- cellPropertyList[[index]]

                # Create the plot
		if (class(umap.df[[CellPropertyName]])=="numeric"){
	                umapCellPropertyPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, alpha= umap.df[[CellPropertyName]], colour = umap.df[[CellPropertyName]])) + geom_point(size = .5) + labs(colour = CellPropertyName,alpha="") + theme_bw()+scale_color_gradient(low="grey",high="blue")

		}
		else{
                umapCellPropertyPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, colour = umap.df[[CellPropertyName]])) + geom_point(size = .5) + labs(colour = CellPropertyName) + theme_bw()
		}

                # Save plot
                ggsave(filename = paste0(folderpath, "/", "CellProperty_", CellPropertyName,"_", filename), plot = umapCellPropertyPlot,width=10,height=8,units="in")
		umapPlots<-c(umapPlots,list(umapCellPropertyPlot))

            }
        }

    }
	return(umapPlots)
}
