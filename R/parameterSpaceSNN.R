#' Generates an overview on the dependence of SNN parameters to the number of clusters
#'
#' @param rca.obj RCA object.
#' @param kL vector of k values to be used. Default is 10-50.
#' @param epsL vector of values indicating the number of neighbours between two nodes such that cells are shared. Default is 5 to 20.
#' @param minPtsL Vector of minimum points with at least eps neighbours to be considered a core point. Default is 5 to 10.
#' @param folderpath Path to store the html file. Defaults to working directory.
#' @param filename name of the html file produced that visuzalizes the parameter space.
#' @return a data frame holding parameter values and resulting cluster numbers.
#' @export
#'
parameterSpaceSNN <- function(rca.obj, kL=base::c(30:50), epsL=base::c(5:20),
                              minPtsL=base::c(5:10), folderpath=".",
                              filename="Graph_based_Clustering_Parameter_Space.html") {

    # Extract projection data
    projection.data <- base::as.matrix(rca.obj$projection.data)


    pcaD = stats::prcomp(projection.data)
    components = base::c(1:(base::max(base::which(base::summary(pcaD)$importance[3,] < 0.99)) + 1))
    d = pcaD$rotation[,components]

    # Obtain parameter space
    kList <- base::c()
    epsList <- base::c()
    minPtsList <- base::c()
    cNumbers <- base::c()
    for (k in kL) {
      for (eps in epsL) {
          for (minPts in minPtsL) {
	        kList <- base::c(kList, k)
		      epsList <- base::c(epsList, eps)
		      minPtsList <- base::c(minPtsList, minPts)
		      clusteringResult <- dbscan::sNNclust(d, k, eps, minPts, borderPoints = T)
		      cNumbers <- base::c(cNumbers,base::length(base::unique(clusteringResult$cluster)))
		      }
	    }
    }

    # Generate a complete data frame
    paramcolors <- grDevices::colorRampPalette(base::c("blue","red"))(base::length(base::unique(cNumbers)))
    cNumbersf <- base::factor(cNumbers)
    hoverInfo <- base::paste0("k: ",kList,"\neps: ",epsList, "\nminPts: ",minPtsList,"\n#clusters: ", cNumbers)
    snnDataO <- base::data.frame(base::cbind(kList,epsList,minPtsList))
    parameterSpace3D <- plotly::plot_ly(data = snnDataO,
	        x = ~kList, y = ~epsList, z = ~minPtsList,
	        color = ~cNumbersf,
	        colors = paramcolors,
	        type = "scatter3d",
	        mode = "markers",
	        marker = base::list(size = 5, width = 2),
	        text = ~hoverInfo,
            hoverinfo = "text")
    htmlwidgets::saveWidget(plotly::as_widget(parameterSpace3D), base::paste0(folderpath, "/",filename))
    pspace <- base::data.frame(base::cbind(snnDataO,cNumbers))
    base::colnames(pspace) <- base::c("K","eps","minPts","clusters")
    return(pspace)

}
