#' Plot bar plots showing the composition of RCA clusters
#'
#' @param rca.obj data matrix (genes x cells)
#' @param deepSplit provides the index of the clustering to use for different cuts in the hierarchical clustering (default 1)
#' @param filename file name of saved heatmap
#'
#' @return figure showing absolute and relative cluster compositions.
#'
#' @examples
#' \dontrun{
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
#' RCA.pbmcs <- dataClust(RCA.pbmcs)
#' plotRCAClusterComposition(RCA.pbmcs,
#' RCA.pbmcs$clustering.out$dynamicColorsList[[1]])
#' }
#'
#' @export
#'


plotRCAClusterComposition <-
    function(rca.obj,
             deepSplit = 1,
             filename = "Cluster_Composition.pdf") {
        Count <- Cluster <- Ratio <- CT <- NULL
        if (!(base::is.null(rca.obj$cell.Type.Estimate))) {
            # Extract projection data and clustering result from RCA object
            heatmapIn = base::as.matrix(rca.obj$projection.data)
            cellTree = rca.obj$clustering.out$cellTree
            clusterColors = rca.obj$clustering.out$dynamicColorsList[[deepSplit]]

            # Compute the composition of each clust with respect to per cell cell type predictions
            enrichmentAll <- c()
            for (type in base::unique(clusterColors)) {
                index = base::which(clusterColors == type)
                enrichmentAll <- base::rbind(enrichmentAll,
                                             (base::cbind(
                                                 type, base::table(base::unlist(rca.obj$cell.Type.Estimate)[index])
                                             )))
            }
            enrichmentAll <-
                base::data.frame(base::cbind(base::row.names(enrichmentAll), enrichmentAll))
            base::colnames(enrichmentAll) <- c("CT", "Cluster", "Count")
            base::rownames(enrichmentAll) <- c(1:dim(enrichmentAll)[1])
            enrichmentAll$Count <-
                base::as.numeric(base::as.character(enrichmentAll$Count))
            totalCounts <-
                base::data.frame(dplyr::count(enrichmentAll, wt = Count, Cluster))
            enrichmentAll <-
                dplyr::left_join(enrichmentAll, totalCounts, by = "Cluster")
            enrichmentAll <-
                base::cbind(enrichmentAll,
                            enrichmentAll$Count / enrichmentAll$n * 100)
            base::colnames(enrichmentAll)[5] <- "Ratio"
            enrichmentAll$Ratio <-
                base::as.numeric(base::as.character(enrichmentAll$Ratio))
            nCols <- ceiling(length(unique(enrichmentAll$CT)) / 26)
            #Generate the cluster composition plots using the randomcolorR package if available
            if (("randomcoloR" %in% base::.packages())) {
                dColors <-
                    randomcoloR::distinctColorPalette(base::length(base::unique(enrichmentAll$CT)))

                ratioPlot <-
                    ggplot2::ggplot(enrichmentAll,
                                    ggplot2::aes(
                                        x = Cluster,
                                        y = Ratio,
                                        fill = CT
                                    )) +
                    ggplot2::geom_bar(stat = "identity") +
                    ggplot2::theme_bw(15) + ggplot2::ylab("Percentage") +
                    ggplot2::coord_flip() + ggplot2::ggtitle("a)") +
                    ggplot2::theme(legend.position = "none") + ggplot2::scale_fill_manual(values =
                                                                                              dColors)
                countPlot <-
                    ggplot2::ggplot(enrichmentAll,
                                    ggplot2::aes(
                                        x = Cluster,
                                        y = Count,
                                        fill = CT
                                    )) +
                    ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw(15) +
                    ggplot2::ylab("Count") +
                    ggplot2::coord_flip() + ggplot2::ggtitle("b)") + ggplot2::scale_fill_manual(values =
                                                                                                    dColors) +
                    ggplot2::labs(fill = "Cell type") + ggplot2::guides(fill = ggplot2::guide_legend(ncol =
                                                                                                         nCols))
            }
            else{
                ratioPlot <-
                    ggplot2::ggplot(enrichmentAll,
                                    ggplot2::aes(
                                        x = Cluster,
                                        y = Ratio,
                                        fill = CT
                                    )) +
                    ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw(15) +
                    ggplot2::ylab("Percentage") +
                    ggplot2::coord_flip() + ggplot2::ggtitle("a)") + ggplot2::theme(legend.position = "none")
                countPlot <-
                    ggplot2::ggplot(enrichmentAll,
                                    ggplot2::aes(
                                        x = Cluster,
                                        y = Count,
                                        fill = CT
                                    )) +
                    ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw(15) +
                    ggplot2::ylab("Count") +
                    ggplot2::coord_flip() + ggplot2::ggtitle("b)") + ggplot2::labs(fill =
                                                                                       "Cell type") +
                    ggplot2::guides(fill = ggplot2::guide_legend(ncol = nCols))
            }
            fig <-
                gridExtra::grid.arrange(ratioPlot, countPlot, widths = c(1, (1.4 + (nCols -
                                                                                        1) / 10), nrow = 1))
            if (!(is.null(filename))) {
                grDevices::pdf(filename,
                               width = 15 + (nCols - 1),
                               height = 7)
                fig
                grDevices::dev.off()
            }
            return(fig)
        } else{
            base::print("Cell specific estimates are not computed yet")
        }
    }
