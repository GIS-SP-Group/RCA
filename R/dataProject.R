#' Compute Reference Component features for clustering analysis
#'
#' @param sc_data data object from a RCA object.
#' @param method Either "GlobalPanel"(default) or "GlobalPanel_CellTypes", "ColonEpiPanel","MonacoPanel","ENCODEMousePanel","ENCODEHumanPanel","ZhangMouseBrainPanel","NovershternPanel","NovershternTCellPanel","MonacoTCellPanel","MonacoBCellPanel","MonacoMonoPanel", "CITESeqPanel" or "Custom"
#' @param customPath directory path (including filename) to any custom panel stored in RDS format. Only used if method == "Custom".
#' @param corMeth Any of the correlation measures supported by R, defaults to pearson
#' @param power power to raise up to for the RCA features before clustering, default is 4
#' @param scale True if the data should be scaled, False otherwise
#' @param min.cell.number.expressing Minimum percent of cells (0 percent -100 percent) expressing a gene such that it is considered in the projection step, default is 1 percent
#' @return a projection matrix.
#'
dataProjectWorker <-
    function(sc_data,
             method = "GlobalPanel",
             customPath = NULL,
             corMeth = "pearson",
             power = 4,
             scale = T,
             min.cell.number.expressing = 1) {
        # If panel for correlation is GlobalPanel

        if (method == "GlobalPanel_CellTypes") {
            # Initialise variable to store projection data from the two fragments of the Global Panel
            projection_list = base::list()

            # For each fragment of the Global Panel

            # Initialise panel
            panel = ReferencePanel[[1]][[1]]

            # Select genes with expression in a minimum number of cells
            geneExpVec <-
                Matrix::rowSums(sc_data > 0) / base::dim(sc_data)[2] * 100
            if (min.cell.number.expressing == 0) {
                shared_genes <-
                    base::intersect(base::rownames(sc_data), base::rownames(panel))
            } else{
                filt.genes <- base::which(geneExpVec < min.cell.number.expressing)
                # Select genes that are shared by the input data and the panel
                if (base::length(filt.genes) > 0) {
                    shared_genes <-
                        base::intersect(base::rownames(sc_data)[-filt.genes],
                                        base::rownames(panel))
                } else {
                    shared_genes <-
                        base::intersect(base::rownames(sc_data),
                                        base::rownames(panel))
                }
            }

            # Reduce the panel and input data to the shared genes
            subset_panel = panel[shared_genes,]
            subset_data = sc_data[shared_genes, , drop = FALSE]

            # For values in the panel below the minimum threshold, set those values to threshold
            subset_panel[subset_panel <= (ReferencePanel$at)[1]] = (ReferencePanel$at)[1]

            # Compute projection of input data with the panel fragment
            if (corMeth == "pearson") {
                subset_panel = base::as.matrix(subset_panel)
                projection_fragment <-
                    qlcMatrix::corSparse(X = subset_panel, Y = subset_data)
            } else {
                projection_fragment <-
                    stats::cor(subset_panel,
                               base::as.matrix(subset_data),
                               method = corMeth)
            }


            # Reattach dimnames
            base::colnames(projection_fragment) <-
                base::colnames(subset_data)
            base::rownames(projection_fragment) <-
                base::colnames(subset_panel)

            # Raise the projection fragment to power
            projection_fragment = base::abs(projection_fragment) ^ (power) * base::sign(projection_fragment)

            # If scaling is required
            if (scale) {
                # Scale
                projection_fragment = base::scale(projection_fragment,
                                                  center = TRUE,
                                                  scale = TRUE)
            }

            # Store projection data of fragment of Global Panel

            # Combine the projection result of multiple Global Panel fragments
            projection = projection_fragment

        }
        # If panel for correlation is ColonEpitheliumPanel

        else if (method == "GlobalPanel") {
            # Initialise variable to store projection data from the two fragments of the Global Panel
            projection_list = base::list()

            # For each fragment of the Global Panel
            for (i in 1:base::length(ReferencePanel[[1]])) {
                # Initialise panel
                panel = ReferencePanel[[1]][[i]]

                # Select genes with expression in a minimum number of cells
                geneExpVec <-
                    Matrix::rowSums(sc_data > 0) / base::dim(sc_data)[2] * 100

                # Select genes that are shared by the input data and the panel
                if (min.cell.number.expressing == 0) {
                    shared_genes <-
                        base::intersect(base::rownames(sc_data),
                                        base::rownames(panel))
                } else{
                    filt.genes <- base::which(geneExpVec < min.cell.number.expressing)
                    # Select genes that are shared by the input data and the panel
                    if (base::length(filt.genes) > 0) {
                        shared_genes <-
                            base::intersect(base::rownames(sc_data)[-filt.genes],
                                            base::rownames(panel))
                    } else {
                        shared_genes <-
                            base::intersect(base::rownames(sc_data),
                                            base::rownames(panel))
                    }

                }

                # Reduce the panel and input data to the shared genes
                subset_panel = panel[shared_genes,]
                subset_data = sc_data[shared_genes, , drop = FALSE]

                # For values in the panel below the minimum threshold, set those values to threshold
                subset_panel[subset_panel <= (ReferencePanel$at)[i]] = (ReferencePanel$at)[i]

                # Compute projection of input data with the panel fragment
                if (corMeth == "pearson") {
                    subset_panel = base::as.matrix(subset_panel)
                    projection_fragment <-
                        qlcMatrix::corSparse(X = subset_panel, Y = subset_data)
                } else {
                    projection_fragment <-
                        stats::cor(subset_panel,
                                   base::as.matrix(subset_data),
                                   method = corMeth)
                }


                # Reattach dimnames
                base::colnames(projection_fragment) <-
                    base::colnames(subset_data)
                base::rownames(projection_fragment) <-
                    base::colnames(subset_panel)

                # Raise the projection fragment to power
                projection_fragment =  base::abs(projection_fragment) ^ (power) *  base::sign(projection_fragment)

                # If scaling is required
                if (scale) {
                    # Scale
                    projection_fragment =  base::scale(projection_fragment,
                                                       center = TRUE,
                                                       scale = TRUE)
                }

                # Store projection data of fragment of Global Panel
                projection_list[[i]] = projection_fragment
            }

            # Combine the projection result of multiple Global Panel fragments
            projection =  base::do.call("rbind", projection_list)

        }
        # If panel for correlation is ColonEpitheliumPanel
        else if (method == "ColonEpitheliumPanel") {
            # Select genes with expression in a minimum number of cells
            geneExpVec <-
                Matrix::rowSums(sc_data > 0) / base::dim(sc_data)[2] * 100
            if (min.cell.number.expressing == 0) {
                shared_genes <-
                    base::intersect(
                        base::rownames(sc_data),
                        base::rownames(ReferencePanel$ColonEpiPanel)
                    )
            } else{
                filt.genes <- base::which(geneExpVec < min.cell.number.expressing)
                # Select genes that are shared by the input data and the panel
                if (base::length(filt.genes) > 0) {
                    shared_genes <-
                        base::intersect(
                            base::rownames(sc_data)[-filt.genes],
                            base::rownames(ReferencePanel$ColonEpiPanel)
                        )
                } else {
                    shared_genes <-
                        base::intersect(
                            base::rownames(sc_data),
                            base::rownames(ReferencePanel$ColonEpiPanel)
                        )
                }
            }

            # Reduce the panel and input data to the shared genes
            subset_panel = ReferencePanel$ColonEpiPanel[shared_genes,]
            subset_data = sc_data[shared_genes, , drop = FALSE]

            # Compute projection of input data with the panel
            if (corMeth == "pearson") {
                subset_panel = base::as.matrix(subset_panel)
                projection <-
                    qlcMatrix::corSparse(X = subset_panel, Y = subset_data)
            } else {
                projection <-
                    stats::cor(subset_panel,
                               base::as.matrix(subset_data),
                               method = corMeth)
            }
            base::rownames(projection) <- base::colnames(subset_panel)
            base::colnames(projection) <- base::colnames(subset_data)

            # Raise the projection to power
            projection = base::abs(projection) ^ (power) * base::sign(projection)

            # If scaling is required
            if (scale) {
                # Scale
                projection = base::scale(projection,
                                         center = TRUE,
                                         scale = TRUE)
            }
        }
        # If any other panel is chosen
        else if (method %in% base::names(ReferencePanel)) {
            panel <- ReferencePanel[[method]]

            # Select genes with expression in a minimum number of cells
            geneExpVec <-
                Matrix::rowSums(sc_data > 0) / base::dim(sc_data)[2] * 100
            if (min.cell.number.expressing == 0) {
                shared_genes <-
                    base::intersect(base::rownames(sc_data), base::rownames(panel))
            } else{
                filt.genes <- base::which(geneExpVec < min.cell.number.expressing)
                # Select genes that are shared by the input data and the panel
                if (base::length(filt.genes) > 0) {
                    shared_genes <-
                        base::intersect(base::rownames(sc_data)[-filt.genes],
                                        base::rownames(panel))
                } else {
                    shared_genes <-
                        base::intersect(base::rownames(sc_data),
                                        base::rownames(panel))
                }
            }

            # Reduce the panel and input data to the shared genes
            subset_panel = panel[shared_genes,]
            subset_data = sc_data[shared_genes, , drop = FALSE]

            # Compute projection of input data with the panel
            if (corMeth == "pearson") {
                subset_panel = base::as.matrix(subset_panel)
                projection <-
                    qlcMatrix::corSparse(X = subset_panel, Y = subset_data)
            } else {
                projection <-
                    stats::cor(subset_panel,
                               base::as.matrix(subset_data),
                               method = corMeth)
            }
            base::rownames(projection) <- base::colnames(subset_panel)
            base::colnames(projection) <- base::colnames(subset_data)

            # Raise the projection to power
            projection = base::abs(projection) ^ (power) * base::sign(projection)

            # If scaling is required
            if (scale) {
                # Scale
                projection = base::scale(projection,
                                         center = TRUE,
                                         scale = TRUE)
            }

        }

        # If no provided method is chosen, it is assumed that the user wishes to use a custom panel
        else {
            # Load panel from path provided
            panel <- base::readRDS(customPath)


            # Select genes with expression in a minimum number of cells
            geneExpVec <-
                Matrix::rowSums(sc_data > 0) / base::dim(sc_data)[2] * 100
            if (min.cell.number.expressing == 0) {
                shared_genes <-
                    base::intersect(base::rownames(sc_data), base::rownames(panel))
            } else{
                filt.genes <- which(geneExpVec < min.cell.number.expressing)
                # Select genes that are shared by the input data and the panel
                if (base::length(filt.genes) > 0) {
                    shared_genes <-
                        base::intersect(base::rownames(sc_data)[-filt.genes],
                                        base::rownames(panel))
                } else {
                    shared_genes <-
                        base::intersect(base::rownames(sc_data),
                                        base::rownames(panel))
                }
            }

            # Reduce the panel and input data to the shared genes
            subset_panel = panel[shared_genes,]
            subset_data = sc_data[shared_genes, , drop = FALSE]

            # Compute projection of input data with the panel
            if (corMeth == "pearson") {
                subset_panel = base::as.matrix(subset_panel)
                projection <-
                    qlcMatrix::corSparse(X = subset_panel, Y = subset_data)
            } else {
                projection <-
                    stats::cor(subset_panel,
                               base::as.matrix(subset_data),
                               method = corMeth)
            }
            base::rownames(projection) <- base::colnames(subset_panel)
            base::colnames(projection) <- base::colnames(subset_data)

            # Raise the projection to power
            projection = base::abs(projection) ^ (power) * base::sign(projection)

            # If scaling is required
            if (scale) {
                # Scale
                projection = base::scale(projection,
                                         center = TRUE,
                                         scale = TRUE)
            }
        }
        # Store projection result as Matrix
        return(methods::as(base::as.matrix(projection), "dgCMatrix"))
    }


#' Compute Reference Component features for clustering analysis
#'
#' @param rca.obj RCA object.
#' @param method Either "GlobalPanel" (default) or "GlobalPanel_CellTypes", "ColonEpiPanel","MonacoPanel","ENCODEMousePanel","ENCODEHumanPanel","ZhangMouseBrainPanel","NovershternPanel","NovershternTCellPanel","MonacoTCellPanel","MonacoBCellPanel","MonacoMonoPanel", "CITESeqPanel" or "Custom"
#' @param customPath directory path (including filename) to any custom panel stored in RDS format. Only used if method == "Custom".
#' @param corMeth Any of the correlation measures supported by R (default pearson)
#' @param power power to raise up to for the RCA features before clustering (default 4)
#' @param scale True if the data should be scaled, False otherwise (default TRUE)
#' @param min.cell.number.expressing Minimum number of cells (0-100) expressing a gene such that it is considered in the projection step (default 1).
#' @return RCA object.
#'
#' @examples
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
#' print(RCA.pbmcs)
#'
#' @export


dataProject <-
    function(rca.obj,
             method = "GlobalPanel",
             customPath = NULL,
             corMeth = "pearson",
             power = 4,
             scale = T,
             min.cell.number.expressing = 1) {
        # Run the worker function
        rca.obj$projection.data <-
            dataProjectWorker(
                rca.obj$data,
                method,
                customPath,
                corMeth,
                power,
                scale,
                min.cell.number.expressing
            )

        ### Return RCA object
        return(rca.obj)
    }


#' Compute Reference Component features for clustering analysis against a list of panels
#'
#' @param rca.obj RCA object.
#' @param method List of panel identifiers containing "GlobalPanel", "GlobalPanel_CellTypes", "ColonEpiPanel","MonacoPanel","ENCODEMousePanel","ENCODEHumanPanel","ZhangMouseBrainPanel","NovershternPanel","NovershternTCellPanel","MonacoTCellPanel","MonacoBCellPanel","MonacoMonoPanel", "CITESeqPanel" or "Custom". Default is: list("NovershternPanel","MonacoPanel","GlobalPanel_CellTypes").
#' @param customPath directory path (including filename) to any custom panel stored in RDS format. Only used if method == "Custom".
#' @param corMeth Any of the correlation measures supported by R (default pearson)
#' @param power power to raise up to for the RCA features before clustering (default 4)
#' @param scale True if the data should be scaled, False otherwise (default TRUE)
#' @param min.cell.number.expressing Minimum number of cells (0-100) expressing a gene such that it is considered in the projection step (default 1).
#' @return RCA object.
#'
#' @examples
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' RCA.pbmcs <- dataProjectMultiPanel(RCA.pbmcs)
#' print(RCA.pbmcs)
#'
#' @export


dataProjectMultiPanel <-
    function(rca.obj,
             method = list("NovershternPanel", "MonacoPanel", "GlobalPanel_CellTypes"),
             customPath = NULL,
             corMeth = "pearson",
             power = 4,
             scale = T,
             min.cell.number.expressing = 1) {
        # Extract data
        tmp <- base::c()
        if (!(base::is.null(method))) {
            for (element in method) {
                tmp <-
                    base::rbind(
                        tmp,
                        dataProjectWorker(
                            rca.obj$data,
                            element,
                            customPath,
                            corMeth,
                            power,
                            T,
                            min.cell.number.expressing
                        )
                    )
            }
        }

        if (!(is.null(customPath))) {
            for (element in customPath) {
                tmp <-
                    base::rbind(
                        tmp,
                        dataProjectWorker(
                            rca.obj$data,
                            "Custom",
                            element,
                            corMeth,
                            power,
                            T,
                            min.cell.number.expressing
                        )
                    )
            }
        }
        # Assign projection result to RCA object
        rca.obj$projection.data <- tmp

        ### Return RCA object
        return(rca.obj)

    }
