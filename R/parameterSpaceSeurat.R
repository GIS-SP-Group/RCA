#' Generates a data frame showcasing the dependence between resolution and cluster quantity
#'
#' @param rca.obj RCA object.
#' @param stepsize Stepsize used to generate an overview on clustering results (default 0.1)
#' @param filename name of the pdf file produced that visuzalizes the parameter space (default NULL)
#' @return list of a a data frame holding parameter values and resulting cluster numbers and the parameter figure
#'
#' @examples
#' \dontrun{
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
#' parameterSpaceSeurat(RCA.pbmcs)
#' }
#'
#' @export
#'
parameterSpaceSeurat <- function(rca.obj, stepsize = 0.1, filename = NULL) {

	# Extract projection data
	nClusters <- base::c()
	stepsize = 0.1
	for (RES in base::seq(0,1,stepsize)) {
		nClusters <- base::rbind(nClusters,
		                       base::cbind(Resolution = RES,
		                                   Clusters = base::length(base::unique(RCAv2::dataSClust(rca.obj,res = RES)$clustering.out$dynamicColorsList$Clusters))))
	}
	nClusters <- base::data.frame(nClusters)
	Resolution <- Clusters <- NULL
	parameterFigure <- ggplot2::ggplot(nClusters, ggplot2::aes(x = Resolution,y = Clusters)) +
		ggplot2::geom_point() +
		ggplot2::geom_line() +
		ggplot2::theme_bw(15) +
		ggplot2::ylab("#Clusters") +
		ggplot2::xlab("Seurat resolution")


	if (!(is.null(filename))) {
	    grDevices::png(filename, width = 800, height = 800)
	    base::print(parameterFigure)
	    grDevices::dev.off()
	}
	return(list(nClusters, parameterFigure))
}
