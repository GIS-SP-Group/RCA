#' Generate elbow plot for a PCA of the projection matrix.
#'
#' @param rca.obj RCA object.
#' @param nPCs Number of PCs to be used if distance should be computed in PC embedding of the projection (default 0, computation of distance in projection space)
#' @param filename Name of figure to be stored. Must end with '.png'.
#' @param approx Flag indicating wheter an approximation of the SVD should be used
#' @return NULL.
#' @export
#'
elbowPlot <- function(rca.obj,nPCs=50,filename="Projection_Elbow.png",approx=F) {

  projection.data <- base::as.matrix(rca.obj$projection.data)

  if (approx){
    tmp<-irlba::irlba(projection.data,nv=nPCs)
    tmp$sdev <- tmp$d/base::sqrt(base::max(1, base::ncol(projection.data) - 1))
  }
  else{
	  tmp<-irlba::prcomp_irlba(projection.data,n=nPCs,center=F,scale.=F)
  }
	grDevices::png(filename)
	graphics::plot(tmp$sdev,ylab="Standard deviation",xlab="Component")
	grDevices::dev.off()
	return()
}


#' Generate cell clusters using Seurat graph based clustering.
#'
#' @param rca.obj RCA object.
#' @param res Resolution parameter (between 0.0 and 1.0) to be used for clustering. Default: 0.5. Check the Seurat documentation for details.
#' @param corMeth Correlation method used to compute the distance matrix of the projection (none(default), pearson , spearman, kendal).
#' @param nPCs Number of PCs to be used if distance should be computed in PC embedding of the projection (default 10)
#' @param approx Using approximation for PCA computation (default T)
#' @return RCA object.
#' @export
#'
dataSClust <- function(rca.obj,res=0.5, corMeth="none", nPCs=10, approx=T) {
 tempS<-Seurat::CreateSeuratObject(rca.obj$raw.data)
 projection.data <- base::as.matrix(rca.obj$projection.data)
 if(!(approx)){
  if ((nPCs==0) & (corMeth !="none")){
	  if (("HiClimR" %in% base::.packages()) & (corMeth=="pearson")){
		  projection<-1-HiClimR::fastCor(projection.data)
		  }else{
			projection<-1-stats::cor(projection.data,method=corMeth)
		  }
  } else{
    if (corMeth != "none"){
		  pca_result<-stats::cor(base::t(irlba::prcomp_irlba(projection.data,n=nPCs,center=F,scale.=F)$rotation),method=corMeth)
		  base::colnames(pca_result)<-base::colnames(projection.data)
		  base::row.names(pca_result)<-base::colnames(projection.data)
		  projection<-stats::as.dist(1-pca_result)
		} else {
		  if ((nPCs != 0)&(corMeth=="none")){
        pca_result<-irlba::prcomp_irlba(projection.data,n=nPCs,center=F,scale.=F)$rotation
        base::row.names(pca_result)<-base::colnames(projection.data)
        projection<-pca_result
		  }
      }
  }
   tempS@reductions[["pca"]]<-methods::new(Class = "DimReduc", cell.embeddings = base::matrix(0,0,0), assay.used = "RNA")
   tempS@reductions$pca@cell.embeddings<-projection
 }else{
   if (corMeth!="none"){
     base::print("Ignoring distance metric in approximative computation")
   }
   #From Seurat RunPCA
   npcs <- base::min(nPCs, base::nrow(projection.data) - 1)
   pca.results <- irlba::irlba(projection.data, nv = npcs)
   feature.loadings <- pca.results$v
   sdev <- pca.results$d/base::sqrt(base::max(1, base::ncol(projection.data) - 1))
   projection <- pca.results$u %*% diag(pca.results$d)

   base::rownames(x = feature.loadings) <- base::rownames(x = projection.data)
   base::colnames(x = feature.loadings) <- base::paste0("PC_", 1:npcs)
   base::rownames(x = projection) <- base::colnames(x = projection.data)
   base::colnames(x = projection) <- base::colnames(x = feature.loadings)
   total.variance <- base::sum(base::apply(X=projection.data,MARGIN=2,FUN=stats::var))
   tempS@reductions[["pca"]] <- Seurat::CreateDimReducObject(
     embeddings = projection,
     loadings = feature.loadings,
     assay = "RNA",
     stdev = sdev,
     key = "PC_",
     misc = base::list(total.variance = total.variance))
	}

	tempS<-Seurat::FindNeighbors(object = tempS)
	tempS<-Seurat::FindClusters(tempS,resolution = res)

	# Convert labels to colours for each tree cut
	if (base::length(base::unique(tempS$seurat_clusters))<41){
		dynamicColorsList<-base::list(WGCNA::labels2colors(tempS$seurat_clusters))
	} else {
		if (("randomcoloR" %in% base::.packages()) & ("plotrix" %in% base::.packages())){
		     clusterColors<-randomcoloR::distinctColorPalette(base::length(base::unique(tempS$seurat_clusters)))
		     clusterColors<-base::sapply(clusterColors,plotrix::color.id)
		     clusterColors<-base::sapply(clusterColors,function(x){return(x[1])})
		     base::names(clusterColors)<-base::unique(tempS$seurat_clusters)
		     dynamicColorsList<-base::list(Colors=clusterColors[base::as.character(tempS$seurat_clusters)])
		    } else{
			dynamicColorsList<-base::list(WGCNA::labels2colors(tempS$seurat_clusters))
			}
	}
	base::names(dynamicColorsList)<-c("Clusters")
	# Assign clustering result to RCA object
	rca.obj$clustering.out <- base::list(
	"cellTree" = tempS$seurat_clusters,
	"dynamicColorsList" = dynamicColorsList
	)
	# Return RCA object
	return(rca.obj)
}
