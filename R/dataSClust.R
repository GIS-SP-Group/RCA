#' Generate elbow plot for a PCA of the projection matrix.
#'
#' @param rca.obj RCA object.
#' @param nPCs Number of PCs to be used if distance should be computed in PC embedding of the projection (default 0, computation of distance in projection space)
#' @param approx Flag indicating wheter an approximation of the SVD should be used
#' @return NULL.
#' @export
#'
elbowPlot <- function(rca.obj,nPCs=50,filename="Projection_Elbow.png",approx=F) {
  if (approx){
    tmp<-irlba::irlba(rca.obj$projection.data,nv=nPCs)
    tmp$sdev <- tmp$d/sqrt(max(1, ncol(rca.obj$projection.data) - 1))
  }
  else{
	  tmp<-irlba::prcomp_irlba(rca.obj$projection.data,n=nPCs,center=F,scale.=F)
  }
	png(filename)
	plot(tmp$sdev,ylab="Standard deviation",xlab="Component")
	dev.off()
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
 if(!(approx)){
  if ((nPCs==0) & (corMeth !="none")){
	  if (require(HiClimR) & (corMeth=="pearson")){
		  projection<-as.dist(1-HiClimR::fastCor(as.matrix(rca.obj$projection.data)))
		  }else{
			projection<-as.dist(1-cor(rca.obj$projection.data,method=corMeth))
		  }
  } else{
    if (corMeth != "none"){
		  pca_result<-cor(t(irlba::prcomp_irlba(rca.obj$projection.data,n=nPCs,center=F,scale.=F)$rotation),method=corMeth)
		  colnames(pca_result)<-colnames(rca.obj$projection.data)
		  row.names(pca_result)<-colnames(rca.obj$projection.data)
		  projection<-as.dist(1-pca_result)
		} else {
		  if ((nPCs != 0)&(corMeth=="none")){
        pca_result<-irlba::prcomp_irlba(rca.obj$projection.data,n=nPCs,center=F,scale.=F)$rotation
        row.names(pca_result)<-colnames(rca.obj$projection.data)
        projection<-pca_result}
      }
  }
   tempS@reductions[["pca"]]<-new(Class = "DimReduc", cell.embeddings = matrix(0,0,0), assay.used = "RNA")
   tempS@reductions$pca@cell.embeddings<-as.matrix(projection)
 }else{
   if (corMeth!="none"){
     print("Ignoring distance metric in approximative computation")
   }
   #From Seurat RunPCA
   npcs <- min(nPCs, nrow(rca.obj$projection.data) - 1)
   pca.results <- irlba::irlba(A =t(rca.obj$projection.data), nv = npcs)
   feature.loadings <- pca.results$v
   sdev <- pca.results$d/sqrt(max(1, ncol(rca.obj$projection.data) - 1))
   projection <- pca.results$u %*% diag(pca.results$d)
   
   rownames(x = feature.loadings) <- rownames(x = rca.obj$projection.data)
   colnames(x = feature.loadings) <- paste0("PC_", 1:npcs)
   rownames(x = projection) <- colnames(x = rca.obj$projection.data)
   colnames(x = projection) <- colnames(x = feature.loadings)
   total.variance <- sum(apply(X=rca.obj$projection.data,MARGIN=2,FUN=var))
   tempS@reductions[["pca"]] <- Seurat::CreateDimReducObject(
     embeddings = projection,
     loadings = feature.loadings,
     assay = "RNA",
     stdev = sdev,
     key = "PC_",
     misc = list(total.variance = total.variance))
	}

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
