#' Plot umap of projection to the RCA panel
#'
#' The presence of cell type estimates, relative ranks and confindence scores are detected automatically and are plotted accordingly.
#' @param rca.obj RCA object
#' @param cellPropertyList list of cell properties to plot
#' @param folderpath path to save umap to
#' @param filename file name of saved umap
#' @param fontsize Size of the font used for plotting
#' @export
#'

plotRCAUMAP <- function(rca.obj, cellPropertyList = NULL, folderpath = ".", filename = "RCA_UMAP.pdf",fontsize=10) {

    # Extract projection data from RCA object
    clusterColorList = rca.obj$clustering.out$dynamicColorsList
    rRank=rca.obj$rRank
    rBaseColors<-rca.obj$baseColors
    confScore=unlist(rca.obj$cScore)

    # Compute UMAP projection from cell type projection
    if (is.null(rca.obj$umap.coordinates)){
	    print("UMAP coordinates have not been computed yet")
	    return(NA)
    }
    else{
        # Store UMAP layout in data frame for plotting
         umap.df <- as.data.frame(rca.obj$umap.coordinates[,1:2])
         colnames(umap.df) <- c("UMAP1", "UMAP2")
    }
    umapPlots<-list()
    # If no cluster colors or cell properties are to be plotted
    if(is.null(clusterColorList) & is.null(cellPropertyList)) {

        # Plot UMAP of cells without annotations
        umap.plot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2)) + geom_point(size = .5) + theme_classic(fontsize)
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
                umapClusterColorsPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, colour = .data[[clusterColorName]])) + geom_point(size = .5) + scale_color_manual(values = sort(unique(umap.df[[clusterColorName]]))) + labs(colour = clusterColorName) + theme_bw(fontsize) + guides(colour = guide_legend(override.aes = list(size=5)))

                # Save plot
                ggsave(filename = paste0(folderpath, "/", "ClusterColors_", clusterColorName,"_", filename), plot = umapClusterColorsPlot,width=9,height=7,units="in")
		umapPlots<-c(umapPlots,list(umapClusterColorsPlot))

            }

        }

      # If cluster rank is to be plotted
        if(!is.null(rRank) & length(rRank) != 0) {

            #Get the name of this cluster annotation
            clusterColorName = names(clusterColorList[index])

            # Set the data frame column to the color vector
            umap.df[[clusterColorName]] <- clusterColorList[[index]]

            # Create the plot
            umapClusterColorsPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, colour = unlist(rRank))) + geom_point(size = .5) + scale_color_identity() +  theme_bw(fontsize) +ggtitle("a)")


	    names(rBaseColors)<-NULL
	    colorOrder<-order(unique(unlist(rBaseColors)))
	    colorVec<-unique(unlist(rBaseColors))[colorOrder]
	    names(colorVec)<-unique(names(unlist(rBaseColors)))[colorOrder]
            umapClusterColorsPlot2 <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, colour = unlist(rBaseColors))) + geom_point(size = .5) + scale_color_identity(labels=names(colorVec),guide="legend") +  theme_bw(fontsize) + ggtitle("b)")+   theme(legend.position="right")+labs(color="Cell type")+guides(colour = guide_legend(override.aes = list(size=4)))

            # Save plot
	    pdf(paste0(folderpath, "/", "RelativeRank_", filename),width=14,height=7)
	    grid.arrange(umapClusterColorsPlot,umapClusterColorsPlot2,widths=c(1,1.2))
	    dev.off()
	    umapPlots<-c(umapPlots,list(umapClusterColorsPlot, umapClusterColorsPlot2))

        }

      # If cluster confidence is to be plotted
        if(!is.null(confScore) & length(confScore) != 0 & length(rca.obj$cell.Type.Estimate != 0)) {

            #Get the name of this cluster annotation
            clusterColorName = names(clusterColorList[index])

            # Set the data frame column to the color vector
            umap.df[[clusterColorName]] <- clusterColorList[[index]]

            # Create the plot
            umapClusterColorsPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, alpha=confScore, colour = umap.df[[clusterColorName]])) + geom_point(size = .5) + theme_bw(fontsize) + ggtitle("a)")+theme(legend.position="none")+ scale_color_manual(values = sort(unique(umap.df[[clusterColorName]]))) + labs(colour = clusterColorName)


	    umapClusterColorsPlot2<-ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, colour = names(rca.obj$cell.Type.Estimate))) +
	    geom_point(size = .5) +  theme(legend.position="right")+labs(color="Cell type")+theme_bw(fontsize)+ggtitle("b)")+
	    guides(colour = guide_legend(override.aes = list(size=4)))+ scale_color_identity(labels=unlist(rca.obj$cell.Type.Estimate),guide="legend")

            # Save plot
	    pdf(paste0(folderpath, "/", "ConfidenceScore_", filename),width=14,height=7)
	    grid.arrange(umapClusterColorsPlot,umapClusterColorsPlot2,widths=c(1,1.2))
	    dev.off()
	    umapPlots<-c(umapPlots,list(umapClusterColorsPlot, umapClusterColorsPlot2))

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
	                umapCellPropertyPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, alpha= umap.df[[CellPropertyName]], colour = umap.df[[CellPropertyName]])) + geom_point(size = .5) + labs(colour = CellPropertyName,alpha="") + theme_bw(fontsize)+scale_color_gradient(low="grey",high="blue")+guides(colour = guide_legend(override.aes = list(size=5)))

		}
		else{
                umapCellPropertyPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, colour = umap.df[[CellPropertyName]])) + geom_point(size = .5) + labs(colour = CellPropertyName) + theme_bw(fontsize)+guides(colour = guide_legend(override.aes = list(size=5)))
		}

                # Save plot
                ggsave(filename = paste0(folderpath, "/", "CellProperty_", CellPropertyName,"_", filename), plot = umapCellPropertyPlot,width=9,height=7,units="in")
		umapPlots<-c(umapPlots,list(umapCellPropertyPlot))

            }
        }

    }
	return(umapPlots)
}
