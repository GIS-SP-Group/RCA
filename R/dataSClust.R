#' Generate cell clusters using Seurat graph based clustering.
#'
#' @param rca.obj RCA object.
#' @param res Resolution parameter (between 0.0 and 1.0) to be used for clustering. Default: 0.5. Check the Seurat documentation for details.
#' @param corMeth Correlation method used to compute the distance matrix of the projection (pearson (default), spearman, kendal).
#' @return RCA object.
#' @export
#'
dataSClust <- function(rca.obj,res=0.5,corMeth="pearson") {
	projection.data <- as.matrix(rca.obj$projection.data)
	tempS<-Seurat::CreateSeuratObject(as.matrix(rca.obj$raw.data))
	if (require(HiClimR) & (corMeth=="pearson")){
		projection<-as.dist(1-HiClimR::fastCor(projection.data))
	}else{
		projection<-as.dist(1-cor(projection.data,method="corMeth"))
	}
	str(projection)
	tempS@reductions[["pca"]]<-new(Class = "DimReduc", cell.embeddings = matrix(0,0,0), assay.used = "RNA")
	tempS@reductions$pca@cell.embeddings<-as.matrix(projection)
	tempS<-Seurat::FindNeighbors(object = tempS)
	tempS<-Seurat::FindClusters(tempS,resolution = res)

	# Convert labels to colours for each tree cut
	if (length(unique(tempS$seurat_clusters))<41){
		dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
	} else {
		if (require(randomcolorR) & require(plotrix)){
		     clusterColors<-randomcoloR::distinctColorPalette(length(unique(tempS$seurat_clusters)))
		     clusterColors<-sapply(clusterColors,plotrix::color.id)
		     clusterColors<-sapply(clusterColors,function(x){return(x[1])})
		     names(clusterColors)<-unique(clusteringResult$cluster)
		     dynamicColorsList<-list(Colors=clusterColors[as.character(clusteringResult$cluster)])
		    } else{
			dynamicColorsList<-list(WGCNA::labels2colors(tempS$seurat_clusters))
			}
	}
	names(dynamicColorsList)<-c("Clusters")
	# Assign clustering result to RCA object
	rca.obj$clustering.out <- list(
	"cellTree" = tempS$seurat_clusters,
	"dynamicColorsList" = dynamicColorsList
	)
	# Return RCA object
	return(rca.obj)
}
