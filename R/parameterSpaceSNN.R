#' Generate cell clusters using hierarchical clustering and dynamic tree cutting.
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
parameterSpaceSNN <- function(rca.obj,kL=c(10:50),epsL=c(5:20),minPtsL=c(5:10),folderpath=".",filename="Graph_based_Clustering_Parameter_Space.html") {

    ### Extract projection data
    projection.data <- as.matrix(rca.obj$projection.data)
    ### Load packages

    # fastcluster
    if (!require(dbscan))
        install.packages("dbscan", repos = "http://cran.us.r-project.org")
    require(dbscan)
    if (!require(ggplot2))
        install.packages("ggplot2", repos = "http://cran.us.r-project.org")
    require(ggplot2)
    if (!require(ggpubr))
        install.packages("ggpubr", repos = "http://cran.us.r-project.org")
    require(ggpubr)
    if (!require(plotly))
        install.packages("plotly", repos = "http://cran.us.r-project.org")
    require(plotly)


    pcaD = prcomp(projection.data)
    components=c(1:(max(which(summary(pcaD)$importance[3,]<0.99))+1))
    d=pcaD$rotation[,components]

    # Obtain parameter space
    kList<-c()
    epsList<-c()
    minPtsList<-c()
    cNumbers<-c()
    for (k in kL){
      for (eps in epsL){
          for (minPts in minPtsL){
	        kList<-c(kList,k)
		      epsList<-c(epsList,eps)
		      minPtsList<-c(minPtsList,minPts)
		      clusteringResult<-sNNclust(d,k,eps,minPts,borderPoints = T)
		      cNumbers<-c(cNumbers,length(unique(clusteringResult$cluster)))
		      }
	    }
    }

    # Generate a complet data frame
    paramcolors<-colorRampPalette(c("blue","red"))(length(unique(cNumbers)))
    cNumbersf<-factor(cNumbers)
    hoverInfo<-paste0("k: ",kList,"\neps: ",epsList, "\nminPts: ",minPtsList,"\n#clusters: ", cNumbers)
    snnDataO<-data.frame(cbind(kList,epsList,minPtsList))
    parameterSpace3D<-plot_ly(data = snnDataO,
	        x = ~kList, y = ~epsList, z = ~minPtsList,
	        color = ~cNumbersf,
	        colors = paramcolors,
	        type = "scatter3d",
	        mode = "markers",
	        marker = list(size = 5, width=2),
	        text=~hoverInfo, #This is that extra column we made earlier for which we will use for cell ID
                hoverinfo="text")
    htmlwidgets::saveWidget(as_widget(parameterSpace3D),  paste0(folderpath, "/",filename))
    pspace<-data.frame(cbind(snnDataO,cNumbers))
    colnames(pspace)<-c("K","eps","minPts","clusters")
    return (pspace)

}
