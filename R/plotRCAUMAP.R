#' Plot heatmap of projection to the RCA panel
#'
#' @param projection data matrix (genes x cells)
#' @param clusterColorList list of cluster colour annotations
#' @param cellPropertyList list of cell properties to plot
#' @param folderpath path to save heatmap to
#' @param filename file name of saved heatmap
#' @export
#' @examples
#'
#' plotRCAUMAP(rca_projection, cellTree, dynamicColorsList, nGeneList, folderpath = ".", filename = "RCA_Heatmap.pdf")
#'

plotRCAUMAP <- function(projection, clusterColorList = NULL, cellPropertyList = NULL, folderpath = ".", filename = "RCA_UMAP.pdf") {
    
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
    
    
    # If no cluster colors or cell properties are to be plotted
    if(is.null(clusterColorList) & is.null(cellPropertyList)) {
        
        # Plot UMAP of cells without annotations
        umap.plot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2)) + geom_point(size = 1.5) + theme_classic()
        
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
                umapClusterColorsPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, colour = umap.df[[clusterColorName]])) + geom_point(size = 1.5) + scale_color_manual(values = sort(unique(umap.df[[clusterColorName]]))) + labs(colour = clusterColorName) + theme_classic() + ggtitle(label = paste("Cluster Colors", clusterColorName))
                
                # Save plot
                ggsave(filename = paste0(folderpath, "/", "ClusterColors_", clusterColorName,"_", filename), plot = umapClusterColorsPlot)
                
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
                umapCellPropertyPlot <- ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, colour = umap.df[[CellPropertyName]])) + geom_point(size = 1.5) + scale_colour_gradient(low = "lightgrey", high = "blue") + labs(colour = CellPropertyName) + theme_classic() + ggtitle(label = paste("Cell Property", CellPropertyName))
                
                # Save plot
                ggsave(filename = paste0(folderpath, "/", "CellProperty_", CellPropertyName,"_", filename), plot = umapCellPropertyPlot)
               
            }
        }

    }
    
}