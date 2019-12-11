#' Plot umap of projection to the RCA panel
#'
#' @param rca.obj RCA object
#' @param cellPropertyList list of cell properties to plot
#' @param folderpath path to save umap to
#' @param filename file name of saved umap
#' @export
#'

plotRCAUMAP3D <- function(rca.obj, cellPropertyList = NULL, folderpath = ".", filename = "RCA_UMAP.pdf") {

    ### Extract projection data from RCA object
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
    if (!require(plotly))
        install.packages("plotly", repos = "http://cran.us.r-project.org")
    require(plotly)

    if (is.null(rca.obj$umap.coordinates)){
	print("UMAP has not been computed yet")
    	return NULL
    }
    if (dim(rca.obj$umap.coordinates)[2]<3){
	print("UMAP projection does not have 3 dimensions")
    	return NULL
    }
    else{
        # Store UMAP layout in data frame for plotting
        umap.df <- as.data.frame(rca.obj$umap.coordinates)
        colnames(umap.df) <- c("UMAP1", "UMAP2","UMAP3")
    }
    # If no cluster colors or cell properties are to be plotted
    if(is.null(clusterColorList) & is.null(cellPropertyList)) {
        # Plot UMAP of cells without annotations
        	umap3dPlot<-plot_ly(data = umap.df,
	            x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
		    type = "scatter3d",
		    mode = "markers",
		    marker = list(size = 5, width=2))

	# Save UMAP
	htmlwidgets::saveWidget(as_widget(umap3dPlot),  paste0(folderpath, "/", filename))


    } else {

        # If cluster colors are to be plotted
        if(!is.null(clusterColorList)) {
	    print("Color by cluster id")
            # Create a list of UMAP plots for each cluster coloring
            for(index in seq_along(clusterColorList)) {

                # Get the name of this cluster annotation
                clusterColorName = names(clusterColorList[index])

                # Set the data frame column to the color vector
                umap.df[[clusterColorName]] <- clusterColorList[[index]]

                # Create the plot
		if (require("randomcoloR")){
		dColors=distinctColorPalette(length(unique(umap.df[[clusterColorName]])))
		umap3dPlot<-plot_ly(data = umap.df,
		            x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
			            color = ~umap.df[[clusterColorName]],
			            colors = dColors,
				            type = "scatter3d",
				            mode = "markers",
					            marker = list(size = 5, width=2))
		}
		else{
		umap3dPlot<-plot_ly(data = umap.df,
		            x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
			            color = ~umap.df[[clusterColorName]],
				            type = "scatter3d",
				            mode = "markers",
					            marker = list(size = 5, width=2))

		}
    
                # Save plot
    	        htmlwidgets::saveWidget(as_widget(umap3dPlot),  paste0(folderpath, "/ClusterColors_",clusterColorName,"_",filename))
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
	                umapCellPropertyPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, alpha= umap.df[[CellPropertyName]], colour = umap.df[[CellPropertyName]])) + geom_point(size = .5) + labs(colour = CellPropertyName,alpha="") + theme_bw()+scale_color_gradient(low="grey",high="blue")+guides(colour = guide_legend(override.aes = list(size=2.5)))

		}
		else{
		if (require("randomColorR")){
		dColors=distinctColorPalette(length(unique(umap.df[[CellPropertyName]])))
		umap3dPlot<-plot_ly(data = umap.df,
		            x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
			            color = ~umap.df[[CellPropertyName]],
			            colors = dColors,
				            type = "scatter3d",
				            mode = "markers",
					            marker = list(size = 5, width=2))
		}
		else{
		umap3dPlot<-plot_ly(data = umap.df,
		            x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
			            color = ~umap.df[[CellPropertyName]],
				            type = "scatter3d",
				            mode = "markers",
					            marker = list(size = 5, width=2))
		}
 
		}

                # Save plot
                htmlwidgets::saveWidget(as_widget(umap3dPlot),  paste0(folderpath, "/CellProperty_",CellPropertyName,"_",filename))
            }
        }

    }
}
