#' Generates a data frame showcasing the dependence between resolution and cluster quantity
#'
#' @param rca.obj RCA object.
#' @param stepsize Stepsize used to generate an overview on clustering results.
#' @param folderpath Path to store the html file. Defaults to working directory.
#' @param filename name of the pdf file produced that visuzalizes the parameter space.
#' @return a data frame holding parameter values and resulting cluster numbers.
#' @export
#'
parameterSpaceSeurat <- function(rca.obj,stepsize=0.1,folderpath="./",filename="Seurat_Parameter_Space.pdf") {

	# Extract projection data
	nClusters<-c()
	stepsize=0.1
	for (RES in seq(0,1,stepsize)){
		nClusters<-rbind(nClusters,cbind(Resolution=RES,Clusters=length(unique(RCAv2::dataSClust(rca.obj,res = RES)$clustering.out$dynamicColorsList$Clusters))))
	}
	nClusters<-data.frame(nClusters)
	parameterFigure<-ggplot2::ggplot(nClusters,ggplot2::aes(x=Resolution,y=Clusters))+
		ggplot2::geom_point()+
		ggplot2::geom_line()+
		ggplot2::theme_bw(15)+
		ggplot2::ylab("#Clusters")+
		ggplot2::xlab("Seurat resolution")

	pdf(paste0(folderpath,"/",filename),width=12,height=12)
	parameterFigure
	dev.off()
	return (parameterFigure)
}
