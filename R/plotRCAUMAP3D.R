#' Plot umap of projection to the RCA panel
#'
#' @param rca.obj RCA object
#' @param cellPropertyList list of cell properties to plot
#' @param folderpath path to save umap to
#' @param filename file name of saved umap
#' @export
#'

plotRCAUMAP3D <-
    function(rca.obj,
             cellPropertyList = NULL,
             folderpath = ".",
             filename = "RCA_UMAP.html") {
        UMAP1 <- UMAP2 <- NULL
        ### Extract projection data from RCA object
        clusterColorList = rca.obj$clustering.out$dynamicColorsList
        rRank = rca.obj$rRank
        rBaseColors <- rca.obj$baseColors

        if (base::is.null(rca.obj$umap.coordinates)) {
            base::print("UMAP has not been computed yet")
            return(NA)
        }
        if (base::dim(rca.obj$umap.coordinates)[2] < 3) {
            base::print("UMAP projection does not have 3 dimensions")
            return(NA)
        }
        else{
            # Store UMAP layout in data frame for plotting
            umap.df <- base::as.data.frame(rca.obj$umap.coordinates)
            base::colnames(umap.df) <- base::c("UMAP1", "UMAP2", "UMAP3")
        }
        # If no cluster colors or cell properties are to be plotted
        if (base::is.null(clusterColorList) & base::is.null(cellPropertyList)) {
            # Plot UMAP of cells without annotations
            umap3dPlot <- plotly::plot_ly(
                data = umap.df,
                x = ~ UMAP1,
                y = ~ UMAP2,
                z = ~ UMAP3,
                type = "scatter3d",
                mode = "markers",
                marker = base::list(size = 5, width = 2)
            )

            # Save UMAP
            htmlwidgets::saveWidget(plotly::as_widget(umap3dPlot),
                                    base::paste0(folderpath, "/", filename))


        } else {
            # If cluster colors are to be plotted
            if (!base::is.null(clusterColorList)) {
                base::print("Color by cluster id")
                # Create a list of UMAP plots for each cluster coloring
                for (index in base::seq_along(clusterColorList)) {
                    # Get the name of this cluster annotation
                    clusterColorName = base::names(clusterColorList[index])

                    # Set the data frame column to the color vector
                    umap.df[[clusterColorName]] <-
                        clusterColorList[[index]]

                    dColors = base::sort(base::unique(umap.df[[clusterColorName]]))
                    # Create the plot
                    umap3dPlot <- plotly::plot_ly(
                        data = umap.df,
                        x = ~ UMAP1,
                        y = ~ UMAP2,
                        z = ~ UMAP3,
                        color = ~ umap.df[[clusterColorName]],
                        colors = dColors,
                        type = "scatter3d",
                        mode = "markers",
                        marker = base::list(size = 5, width = 2),
                        text = ~ umap.df[[clusterColorName]],
                        hoverinfo = "text"
                    )


                    # Save plot
                    htmlwidgets::saveWidget(
                        plotly::as_widget(umap3dPlot),
                        base::paste0(
                            folderpath,
                            "/ClusterColors_",
                            clusterColorName,
                            "_",
                            filename
                        )
                    )
                }

            }


            if (!base::is.null(rRank) & base::length(rRank) != 0) {
                # Get the name of this cluster annotation
                clusterColorName = base::names(clusterColorList[index])


                umap.df[[clusterColorName]] <- clusterColorList[[index]]

                # Create the plot
                umapClusterColorsPlot <-
                    ggplot2::ggplot(data = umap.df,
                           mapping = ggplot2::aes(
                               x = UMAP1,
                               y = UMAP2,
                               colour = base::unlist(rRank)
                           )) + ggplot2::geom_point(size = .5) + ggplot2::scale_color_identity() +
                    ggplot2::theme_bw() + ggplot2::ggtitle("a)")


                base::names(rBaseColors) <- NULL
                colorOrder <- base::order(base::unique(base::unlist(rBaseColors)))
                colorVec <- base::unique(base::unlist(rBaseColors))[colorOrder]
                base::names(colorVec) <- base::unique(base::names(base::unlist(rBaseColors)))[colorOrder]
                umap3dPlot <- plotly::plot_ly(
                    data = umap.df,
                    x = ~ UMAP1,
                    y = ~ UMAP2,
                    z = ~ UMAP3,
                    color = ~ base::unlist(rRank),
                    colors = colorVec,
                    type = "scatter3d",
                    mode = "markers",
                    marker =  base::list(size = 5, width = 2),
                    text = ~  base::names(base::unlist(rRank)),
                    hoverinfo = "text"
                )


                # Save plot
                htmlwidgets::saveWidget(
                    plotly::as_widget(umap3dPlot),
                    base::paste0(
                        folderpath,
                        "/RelativeColoring_CellTypes_",
                        filename
                    )
                )


            }




            # If cell properties are to be plotted
            if (! base::is.null(cellPropertyList)) {
                # Create a list of UMAP plots for each cluster coloring
                for (index in  base::seq_along(cellPropertyList)) {
                    # Get the name of this cluster annotation
                    CellPropertyName =  base::names(cellPropertyList)[index]

                    # Set the data frame column to the color vector
                    umap.df[[CellPropertyName]] <-
                        cellPropertyList[[index]]

                    # Create the plot
                    if ( base::class(umap.df[[CellPropertyName]]) == "numeric" ) {
                        umapCellPropertyPlot <-
                            ggplot2::ggplot(
                                data = umap.df,
                                mapping = ggplot2::aes(
                                    x = UMAP1,
                                    y = UMAP2,
                                    alpha = umap.df[[CellPropertyName]],
                                    colour = umap.df[[CellPropertyName]]
                                )
                            ) + ggplot2::geom_point(size = .5) +
                            ggplot2::labs(colour = CellPropertyName, alpha = "") + ggplot2::theme_bw() +
                            ggplot2::scale_color_gradient(low = "grey", high = "blue") +
                            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = base::list(size = 2.5)))

                    }
                    else{
                        if (("randomcoloR" %in% base::.packages())) {
                            dColors = randomcoloR::distinctColorPalette(base::length(base::unique(umap.df[[CellPropertyName]])))
                            umap3dPlot <- plotly::plot_ly(
                                data = umap.df,
                                x = ~ UMAP1,
                                y = ~ UMAP2,
                                z = ~ UMAP3,
                                color = ~ umap.df[[CellPropertyName]],
                                colors = dColors,
                                type = "scatter3d",
                                mode = "markers",
                                marker = base::list(size = 5, width = 2),
                                text = ~ umap.df[[CellPropertyName]],
                                hoverinfo = "text"
                            )
                        }
                        else{
                            umap3dPlot <- plotly::plot_ly(
                                data = umap.df,
                                x = ~ UMAP1,
                                y = ~ UMAP2,
                                z = ~ UMAP3,
                                color = ~ umap.df[[CellPropertyName]],
                                type = "scatter3d",
                                mode = "markers",
                                marker = base::list(size = 5, width = 2),
                                text = ~ umap.df[[CellPropertyName]],
                                hoverinfo = "text"
                            )
                        }

                    }

                    # Save plot
                    htmlwidgets::saveWidget(
                        plotly::as_widget(umap3dPlot),
                        base::paste0(
                            folderpath,
                            "/CellProperty_",
                            CellPropertyName,
                            "_",
                            filename
                        )
                    )
                }
            }

        }
    }
