#' Plot umap of projection to the RCA panel
#'
#' The presence of cell type estimates, relative ranks and confindence scores are detected automatically and are plotted accordingly.
#' @param rca.obj RCA object
#' @param cellPropertyList list of cell properties to plot
#' @param filename file name of saved umap
#' @param fontsize Size of the font used for plotting
#'
#' @return list of UMAP plots
#'
#' @examples
#' \dontrun{
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
#' RCA.pbmcs <- computeUMAP(RCA.pbmcs)
#' plotRCAUMAP(RCA.pbmcs)
#' }
#'
#' @export
#'

plotRCAUMAP <-
    function(rca.obj,
             cellPropertyList = NULL,
             filename = "RCA_UMAP.pdf",
             fontsize = 10) {
        UMAP1 <- UMAP2 <- .data <- NULL
        # Extract projection data from RCA object
        clusterColorList = rca.obj$clustering.out$dynamicColorsList
        rRank = rca.obj$rRank
        rBaseColors <- rca.obj$baseColors
        confScore = base::unlist(rca.obj$cScore)

        # Compute UMAP projection from cell type projection
        if (base::is.null(rca.obj$umap.coordinates)) {
            base::print("UMAP coordinates have not been computed yet")
            return(NA)
        }
        else{
            # Store UMAP layout in data frame for plotting
            umap.df <- base::as.data.frame(rca.obj$umap.coordinates[, 1:2])
            base::colnames(umap.df) <- base::c("UMAP1", "UMAP2")
        }
        umapPlots <- base::list()
        # If no cluster colors or cell properties are to be plotted
        if (base::is.null(clusterColorList) & base::is.null(cellPropertyList)) {
            # Plot UMAP of cells without annotations
            umap.plot <-
                ggplot2::ggplot(data = umap.df,
                       mapping = ggplot2::aes(x = UMAP1, y = UMAP2)) +
                ggplot2::geom_point(size = .5) +
                ggplot2::theme_classic(fontsize)
            umapPlots <- c(umapPlots, list(umap.plot))
            # Save UMAP
            ggplot2::ggsave(filename = base::paste0(filename),
                   plot = umap.plot)

        } else {
            # If cluster colors are to be plotted
            if (!base::is.null(clusterColorList)) {
                # Create a list of UMAP plots for each cluster coloring
                for (index in base::seq_along(clusterColorList)) {
                    # Get the name of this cluster annotation
                    clusterColorName = base::names(clusterColorList[index])

                    # Set the data frame column to the color vector
                    umap.df[[clusterColorName]] <-
                        clusterColorList[[index]]

                    # Create the plot
                    umapClusterColorsPlot <-
                        ggplot2::ggplot(
                            data = umap.df,
                            mapping = ggplot2::aes(
                                x = UMAP1,
                                y = UMAP2,
                                colour = .data[[clusterColorName]]
                            )
                        ) + ggplot2::geom_point(size = .5) +
                        ggplot2::scale_color_manual(values = base::sort(base::unique(umap.df[[clusterColorName]]))) +
                        ggplot2::labs(colour = clusterColorName) + ggplot2::theme_bw(fontsize) +
                        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = base::list(size = 5)))

                    # Save plot
                    ggplot2::ggsave(
                        filename = base::paste0(
                            "ClusterColors_",
                            clusterColorName,
                            "_",
                            filename
                        ),
                        plot = umapClusterColorsPlot,
                        width = 9,
                        height = 7,
                        units = "in"
                    )
                    umapPlots <- base::c(umapPlots, base::list(umapClusterColorsPlot))

                }

            }

            # If cluster rank is to be plotted
            if ( base::length(rRank) != 0) {
                #Get the name of this cluster annotation
                clusterColorName = base::names(clusterColorList[index])

                # Set the data frame column to the color vector
                umap.df[[clusterColorName]] <- clusterColorList[[index]]

                # Create the plot
                umapClusterColorsPlot <-
                    ggplot2::ggplot(data = umap.df,
                           mapping = ggplot2::aes(
                               x = UMAP1,
                               y = UMAP2,
                               colour = base::unlist(rRank)
                           )) + ggplot2::geom_point(size = .5) + ggplot2::scale_color_identity() +  ggplot2::theme_bw(fontsize) +
                    ggplot2::ggtitle("a)")


                base::names(rBaseColors) <- NULL
                colorOrder <-  base::order( base::unique( base::unlist(rBaseColors)))
                colorVec <-  base::unique( base::unlist(rBaseColors))[colorOrder]
                base::names(colorVec) <-  base::unique( base::names( base::unlist(rBaseColors)))[colorOrder]
                umapClusterColorsPlot2 <-
                    ggplot2::ggplot(data = umap.df,
                           mapping = ggplot2::aes(
                               x = UMAP1,
                               y = UMAP2,
                               colour = base::unlist(rBaseColors)
                           )) + ggplot2::geom_point(size = .5) +
                    ggplot2::scale_color_identity(labels = names(colorVec), guide = "legend") +
                    ggplot2::theme_bw(fontsize) + ggplot2::ggtitle("b)") +
                    ggplot2::theme(legend.position = "right") +
                    ggplot2::labs(color = "Cell type") +
                    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = base::list(size = 4)))

                # Save plot
                grDevices::pdf(
                    base::paste0("RelativeRank_", filename),
                    width = 14,
                    height = 7
                )
                gridExtra::grid.arrange(umapClusterColorsPlot,
                             umapClusterColorsPlot2,
                             widths = base::c(1, 1.2))
                grDevices::dev.off()
                umapPlots <-
                    base::c(umapPlots,
                            base::list(umapClusterColorsPlot, umapClusterColorsPlot2))

            }

            # If cluster confidence is to be plotted
            if (base::length(confScore) != 0 & (base::length(rca.obj$cell.Type.Estimate.per.cluster) != 0)) {
                #Get the name of this cluster annotation
                clusterColorName = base::names(clusterColorList[index])

                # Set the data frame column to the color vector
                umap.df[[clusterColorName]] <- clusterColorList[[index]]

                # Create the plot
                umapClusterColorsPlot <-
                    ggplot2::ggplot(
                        data = umap.df,
                        mapping = ggplot2::aes(
                            x = UMAP1,
                            y = UMAP2,
                            alpha = confScore,
                            colour = umap.df[[clusterColorName]]
                        )
                    ) + ggplot2::geom_point(size = .5) + ggplot2::theme_bw(fontsize) +
                    ggplot2::ggtitle("a)") + ggplot2::theme(legend.position = "none") +
                    ggplot2::scale_color_manual(values = base::sort(base::unique(umap.df[[clusterColorName]]))) +
                    ggplot2::labs(colour = clusterColorName)


                umapClusterColorsPlot2 <-
                    ggplot2::ggplot(data = umap.df,
                           mapping = ggplot2::aes(
                               x = UMAP1,
                               y = UMAP2,
                               colour = names(rca.obj$cell.Type.Estimate.per.cluster)
                           )) +
                    ggplot2::geom_point(size = .5) +  ggplot2::theme(legend.position = "right") +
                    ggplot2::labs(color = "Cell type") +
                    ggplot2::theme_bw(fontsize) + ggplot2::ggtitle("b)") +
                    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = base::list(size = 4))) +
                    ggplot2::scale_color_identity(labels = base::unlist(rca.obj$cell.Type.Estimate.per.cluster), guide = "legend")

                # Save plot
                grDevices::pdf(
                    base::paste0("ConfidenceScore_", filename),
                    width = 14,
                    height = 7
                )
                gridExtra::grid.arrange(umapClusterColorsPlot,
                             umapClusterColorsPlot2,
                             widths = base::c(1, 1.2))
                grDevices::dev.off()
                umapPlots <-
                    base::c(umapPlots,
                            base::list(umapClusterColorsPlot, umapClusterColorsPlot2))

            }





            # If cell properties are to be plotted
            if (!base::is.null(cellPropertyList)) {
                # Create a list of UMAP plots for each cluster coloring
                for (index in base::seq_along(cellPropertyList)) {
                    # Get the name of this cluster annotation
                    CellPropertyName = base::names(cellPropertyList)[index]

                    # Set the data frame column to the color vector
                    umap.df[[CellPropertyName]] <-
                        cellPropertyList[[index]]

                    # Create the plot
                    if (base::class(umap.df[[CellPropertyName]]) == "numeric") {
                        umapCellPropertyPlot <-
                            ggplot2::ggplot(
                                data = umap.df,
                                mapping = ggplot2::aes(
                                    x = UMAP1,
                                    y = UMAP2,
                                    alpha = umap.df[[CellPropertyName]],
                                    colour = umap.df[[CellPropertyName]]
                                )
                            ) + ggplot2::geom_point(size = .5) + ggplot2::labs(colour = CellPropertyName, alpha = "") +
                            ggplot2::theme_bw(fontsize) +
                            ggplot2::scale_color_gradient(low = "grey", high = "blue") +
                            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = base::list(size =  5)))

                    }
                    else{
                        umapCellPropertyPlot <-
                            ggplot2::ggplot(
                                data = umap.df,
                                mapping = ggplot2::aes(
                                    x = UMAP1,
                                    y = UMAP2,
                                    colour = umap.df[[CellPropertyName]]
                                )
                            ) + ggplot2::geom_point(size = .5) + ggplot2::labs(colour = CellPropertyName) + ggplot2::theme_bw(fontsize) +
                            ggplot2::guides(colour = ggplot2::guide_legend(override.aes = base::list(size = 5)))
                    }

                    # Save plot
                    ggplot2::ggsave(
                        filename = base::paste0(
                            "CellProperty_",
                            CellPropertyName,
                            "_",
                            filename
                        ),
                        plot = umapCellPropertyPlot,
                        width = 9,
                        height = 7,
                        units = "in"
                    )
                    umapPlots <- base::c(umapPlots, base::list(umapCellPropertyPlot))

                }
            }

        }
        return(umapPlots)
    }
