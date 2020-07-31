#' Plot heatmap of projection to the RCA panel
#'
#' @param rca.obj data matrix (genes x cells)
#' @param var.thrs Minimum threshold of variance for cell type correlations.
#' @param width width of plot in inches. Default is 20.
#' @param height height of plot in inches. Default is 20.
#' @param folderpath path to save heatmap to
#' @param filename file name of saved heatmap
#' @param extraCellProperty vector indicating cell property to be plotted in heatmap
#' @param SeuratColorScheme Use the color scheme known from Seurat (default). Otherwise, use a blue to red color scheme
#'
#' @export
#'

plotRCAHeatmap <- function(rca.obj, var.thrs = 0.1, width = 20, height = 20, folderpath = ".", filename = "RCA_Heatmap.pdf", extraCellProperty = NULL,SeuratColorScheme = TRUE) {
    require(ComplexHeatmap)
    # Extract projection data and clustering result from RCA object
    heatmapIn = as.matrix(rca.obj$projection.data)
    cellTree = rca.obj$clustering.out$cellTree
    clusterColorList = rca.obj$clustering.out$dynamicColorsList

    # Subset projection data to remove unnecessary cell types
    varVec <- apply(heatmapIn, 1, var)


    # Set color scheme of heatmap
    if (SeuratColorScheme){
    colorScheme <-  colorRamp2(c(min(heatmapIn), 0, max(heatmapIn)),c("purple", "black", "yellow"))
    }else{
    colorScheme <-
        circlize::colorRamp2(
            c(min(
                heatmapIn
            ), 0, max(
                heatmapIn
            )),c("#7777FF",
              "white",
              "red")
        )
    }

    heatmapIn <- heatmapIn[varVec >= var.thrs, ]

    if (class(cellTree) == "hclust"){
    # If no cluster colors or cell properties are to be plotted
    if(is.null(clusterColorList)) {

        # Initialize heatmap object
        ht <- ComplexHeatmap::Heatmap(
            matrix = heatmapIn,
            col = colorScheme,

            cluster_columns = as.dendrogram(cellTree),
            cluster_rows = TRUE,

            column_dend_side = "top",
            row_dend_side = "left",

            show_row_dend = TRUE,
            row_dend_width = grid::unit(30, "mm"),

            show_column_dend = TRUE,
            column_dend_height = grid::unit(100, "mm"),
            column_dend_reorder = FALSE,

            show_column_names = FALSE,

            show_row_names = TRUE,
            row_names_gp = grid::gpar(fontsize = 15),

            heatmap_legend_param = list(title = "legend", color_bar = "continuous"),
            use_raster = TRUE,
            raster_device = "png",
            raster_quality = 1
        )

    } else {

        # Ensure each cluster color list is a list of named vectors
        for(index in 1:length(clusterColorList)) {
            names(clusterColorList[[index]]) <- clusterColorList[[index]]
        }

        clusterColorDf <- data.frame(clusterColorList)
        names(clusterColorDf) <- names(clusterColorList)

        # Create cell property list - list of NODG, nUMI and percent.mito
        nUMI <- Matrix::colSums(rca.obj$raw.data)
        nodg <- Matrix::colSums(rca.obj$raw.data > 0)

        mito.genes = grep(pattern = "^MT-", x = rownames(rca.obj$raw.data), value = T)
        if(length(mito.genes) == 0) {
	    if (is.null(extraCellProperty)){
	            cellPropertyList <- list("nUMI" = nUMI, "NODG" = nodg)
	    }else{
	            cellPropertyList <- list("nUMI" = nUMI, "NODG" = nodg,"extra" = extraCellProperty)
	    }
	} else {
	        pMito <- Matrix::colSums(rca.obj$raw.data[mito.genes, ])/Matrix::colSums(rca.obj$raw.data)
		if (is.null(extraCellProperty)){
	            cellPropertyList <- list("nUMI" = nUMI, "NODG" = nodg, "pMito" = pMito)
		}else{
	            cellPropertyList <- list("nUMI" = nUMI, "NODG" = nodg, "pMito" = pMito, "extra" = extraCellProperty)
		}
	}


        # Create list of annotation bar plots from cell property list
        annoBarPlotList <- lapply(cellPropertyList, function(cellPropertyVec){
            ComplexHeatmap::anno_barplot(
                cellPropertyVec,
                gp = grid::gpar(fill = "#777777", col = "#777777"),
                axis = TRUE,
                axis_param = list(side = "right"),
                which = "column"
            )
        })

        # Set names of annotation bar plots
        names(annoBarPlotList) <- names(cellPropertyList)

        # Set parameter list for dynamic number of cell property plots
        paramList <- list(df = clusterColorDf,
                          col = clusterColorList,
                          show_annotation_name = TRUE,
                          annotation_name_side = "left",
                          gap = grid::unit(5, "mm"),
                          which = "column")

        # Add cell property plots
        paramList <- append(paramList, annoBarPlotList)

        # Create HeatmapAnnotation object
        columnColorBar <- do.call(what = ComplexHeatmap::HeatmapAnnotation, args = paramList)

        # Initialize heatmap object
        ht <- ComplexHeatmap::Heatmap(
            matrix = heatmapIn,
            col = colorScheme,

            cluster_columns = as.dendrogram(cellTree),
            cluster_rows = TRUE,

            column_dend_side = "top",
            row_dend_side = "left",

            show_row_dend = TRUE,
            row_dend_width = grid::unit(30, "mm"),

            show_column_dend = TRUE,
            column_dend_height = grid::unit(100, "mm"),
            column_dend_reorder = FALSE,

            show_column_names = FALSE,

            show_row_names = TRUE,
            row_names_gp = grid::gpar(fontsize = 10),

            top_annotation = columnColorBar,

            heatmap_legend_param = list(title = "legend", color_bar = "continuous"),
            use_raster = TRUE,
            raster_device = "png",
            raster_quality = 1
        )
    }
    } else{
    #graph based clustering
    if(is.null(clusterColorList)) {
        # Initialize heatmap object
        ht <- ComplexHeatmap::Heatmap(
            matrix = heatmapIn,
            col = colorScheme,

            cluster_columns = FALSE,
	    column_order = order(cellTree),
            cluster_rows = TRUE,

            row_dend_side = "left",
            show_row_dend = TRUE,
            row_dend_width = grid::unit(30, "mm"),

            show_column_dend = FALSE,
            show_column_names = FALSE,

            show_row_names = TRUE,
            row_names_gp = grid::gpar(fontsize = 15),

            heatmap_legend_param = list(title = "legend", color_bar = "continuous"),
            use_raster = TRUE,
            raster_device = "png",
            raster_quality = 1
        )

    } else {

        # Ensure each cluster color list is a list of named vectors
        for(index in 1:length(clusterColorList)) {
            names(clusterColorList[[index]]) <- clusterColorList[[index]]
        }
        clusterColorDf <- data.frame(clusterColorList)
        names(clusterColorDf) <- names(clusterColorList)

        # Create cell property list - list of NODG, nUMI and percent.mito
        nUMI <- Matrix::colSums(rca.obj$raw.data)
        nodg <- Matrix::colSums(rca.obj$raw.data > 0)

        mito.genes = grep(pattern = "^MT-", x = rownames(rca.obj$raw.data), value = T)
        pMito <- Matrix::colSums(rca.obj$raw.data[mito.genes, ])/Matrix::colSums(rca.obj$raw.data)

	if (is.null(extraCellProperty)){
	        cellPropertyList <- list("nUMI" = nUMI, "NODG" = nodg, "pMito" = pMito)
	}else{
		cellPropertyList <- list("nUMI" = nUMI, "NODG" = nodg, "pMito" = pMito,"extra" = extraCellProperty)
	}

        # Create list of annotation bar plots from cell property list
        annoBarPlotList <- lapply(cellPropertyList, function(cellPropertyVec){
            ComplexHeatmap::anno_barplot(
                cellPropertyVec,
                gp = grid::gpar(fill = "#777777", col = "#777777"),
                axis = TRUE,
                axis_param = list(side = "right"),
                which = "column"
            )
        })

        # Set names of annotation bar plots
        names(annoBarPlotList) <- names(cellPropertyList)

        # Set parameter list for dynamic number of cell property plots
        paramList <- list(df = clusterColorDf,
                          col = clusterColorList,
                          show_annotation_name = TRUE,
                          annotation_name_side = "left",
                          gap = grid::unit(5, "mm"),
                          which = "column")

        # Add cell property plots
        paramList <- append(paramList, annoBarPlotList)

        # Create HeatmapAnnotation object
        columnColorBar <- do.call(what = ComplexHeatmap::HeatmapAnnotation, args = paramList)

        # Initialize heatmap object
        ht <- ComplexHeatmap::Heatmap(
            matrix = heatmapIn,
            col = colorScheme,

	    cluster_columns = FALSE,
	    column_order = order(cellTree),

            cluster_rows = TRUE,

            row_dend_side = "left",

            show_row_dend = TRUE,
            row_dend_width = grid::unit(30, "mm"),

            show_column_dend = FALSE,
            show_column_names = FALSE,

            show_row_names = TRUE,
            row_names_gp = grid::gpar(fontsize = 10),

            top_annotation = columnColorBar,

            heatmap_legend_param = list(title = "legend", color_bar = "continuous"),
            use_raster = TRUE,
            raster_device = "png",
            raster_quality = 1
        )

    }

    }
    # Create pdf object to hold heatmap
    pdf(paste0(folderpath, "/", filename),
        width = width,
        height = height)

    # Draw heatmap inside pdf device
    ComplexHeatmap::draw(ht,
         heatmap_legend_side = "left",
         annotation_legend_side = "left")

    # Shut down device
    dev.off()

}
