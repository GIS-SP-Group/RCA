#' Plot heatmap of projection to the RCA panel
#'
#' @param rca.obj data matrix (genes x cells)
#' @param cellPropertyList list of cell properties to plot
#' @param folderpath path to save heatmap to
#' @param filename file name of saved heatmap
#' @export
#'

plotRCAHeatmap <- function(rca.obj, cellPropertyList = NULL, folderpath = ".", filename = "RCA_Heatmap.pdf") {

    ### Extract projection data and clustering result from RCA object
    heatmapIn = as.matrix(rca.obj$projection.data)
    cellTree = rca.obj$clustering.out$cellTree
    clusterColorList = rca.obj$clustering.out$dynamicColorsList

    ### Check if package dependencies are available; if not, download from CRAN and require those packages
    # ComplexHeatmap
    if (!require(ComplexHeatmap))
        install.packages("ComplexHeatmap", repos = "http://cran.us.r-project.org")
    require(ComplexHeatmap)

    # ComplexHeatmap
    if (!require(circlize))
        install.packages("circlize", repos = "http://cran.us.r-project.org")
    require(circlize)

    # Set color scheme of heatmap
    colorScheme <-
        colorRamp2(
            seq(-max(abs(
                heatmapIn
            )), max(abs(
                heatmapIn
            )), length.out = 9),
            c(
                "#00007F",
                "blue",
                "#007FFF",
                "cyan",
                "#7FFF7F",
                "yellow",
                "#FF7F00",
                "red",
                "#7F0000"
            )
        )

    # If no cluster colors or cell properties are to be plotted
    if(is.null(clusterColorList) & is.null(cellPropertyList)) {

        # Initialize heatmap object
        ht <- Heatmap(
            matrix = heatmapIn,
            col = colorScheme,

            cluster_columns = as.dendrogram(cellTree),
            cluster_rows = TRUE,

            column_dend_side = "top",
            row_dend_side = "left",

            show_row_dend = TRUE,
            row_dend_width = unit(30, "mm"),

            show_column_dend = TRUE,
            column_dend_height = unit(100, "mm"),
            column_dend_reorder = FALSE,

            show_column_names = TRUE,

            show_row_names = TRUE,
            row_names_gp = gpar(fontsize = 15),

            heatmap_legend_param = list(title = "legend", color_bar = "continuous"),
            use_raster = TRUE,
            raster_device = "png",
            raster_quality = 1
        )

    } else {

        # If cluster colors are to be plotted
        if(!is.null(clusterColorList)) {

            # Ensure each cluster color list is a list of named vectors
            for(index in 1:length(clusterColorList)) {
                names(clusterColorList[[index]]) <- clusterColorList[[index]]
            }

            clusterColorDf <- data.frame(clusterColorList)
            names(clusterColorDf) <- names(clusterColorList)
        }

        # If cell properties are to be plotted
        if(!is.null(cellPropertyList)) {

            # Create list of annotation bar plots from cell property list
            annoBarPlotList <- lapply(cellPropertyList, function(cellPropertyVec){
                anno_barplot(
                    cellPropertyVec,
                    gp = gpar(fill = "#777777", col = "#777777"),
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
                              gap = unit(5, "mm"),
                              which = "column")

            # Add cell property plots
            paramList <- append(paramList, annoBarPlotList)

            # Create HeatmapAnnotation object
            columnColorBar <- do.call(what = HeatmapAnnotation, args = paramList)

        } else {

            # Create HeatmapAnnotation object without cell property plots
            columnColorBar <-
                HeatmapAnnotation(
                    df = clusterColorDf,
                    col = clusterColorList,
                    show_annotation_name = TRUE,
                    annotation_name_side = "left",
                    gap = unit(5, "mm"),
                    which = "column"
                )
        }

        # Initialize heatmap object
        ht <- Heatmap(
            matrix = heatmapIn,
            col = colorScheme,

            cluster_columns = as.dendrogram(cellTree),
            cluster_rows = TRUE,

            column_dend_side = "top",
            row_dend_side = "left",

            show_row_dend = TRUE,
            row_dend_width = unit(30, "mm"),

            show_column_dend = TRUE,
            column_dend_height = unit(100, "mm"),
            column_dend_reorder = FALSE,

            show_column_names = TRUE,

            show_row_names = TRUE,
            row_names_gp = gpar(fontsize = 5),

            top_annotation = columnColorBar,

            heatmap_legend_param = list(title = "legend", color_bar = "continuous"),
            use_raster = TRUE,
            raster_device = "png",
            raster_quality = 1
        )

    }

    # Create pdf object to hold heatmap
    pdf(paste0(folderpath, "/", filename),
        width = 20,
        height = 20)

    # Draw heatmap inside pdf device
    draw(ht,
         heatmap_legend_side = "left",
         annotation_legend_side = "left")

    # Shut down device
    dev.off()

}
