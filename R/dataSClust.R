#Modified from https://www.r-bloggers.com/2013/02/bigcor-large-correlation-matrices-in-r/
bigcor <- function(x, nblocks = 10, verbose = TRUE, corMeth="pearson",...)
{
  require(ff)
  NCOL <- ncol(x)
    ## test if ncol(x) %% nblocks gives remainder 0
    if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
    ## preallocate square matrix of dimension
    ## ncol(x) in 'ff' single format
    corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
      ## split column numbers into 'nblocks' groups
      SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
      ## create all unique combinations of blocks
      COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
        COMBS <- t(apply(COMBS, 1, sort))
        COMBS <- unique(COMBS)
	  ## iterate through each block combination, calculate correlation matrix
	  ## between blocks and store them in the preallocated matrix on both
	  ## symmetric sides of the diagonal
	  for (i in 1:nrow(COMBS)) {
		COMB <- COMBS[i, ]
		G1 <- SPLIT[[COMB[1]]]
	        G2 <- SPLIT[[COMB[2]]]
		COR <- cor(x[, G1], x[, G2], method=corMeth)
		corMAT[G1, G2] <- COR
		corMAT[G2, G1] <- t(COR)
		COR <- NULL
		}
	  gc()
	  return(corMAT)
}

#' Generate cell clusters using Seurat graph based clustering.
#'
#' @param rca.obj RCA object.
#' @param res Resolution parameter (between 0.0 and 1.0) to be used for clustering. Default: 0.5. Check the Seurat documentation for details.
#' @param corMeth Correlation method used to compute the distance matrix of the projection (pearson (default), spearman, kendal).
#' @param bigCor Number of blocks to be used in the divide and conquer bigcor function to compute the correlation distance (advisable for > 100.000 cells if PC embedding is not used, default 0)
#' @param nPCs Number of PCs to be used if distance should be computed in PC embedding of the projection (default 0, computation of distance in projection space)
#' @return RCA object.
#' @export
#'
dataSClust <- function(rca.obj,res=0.5,corMeth="pearson",bigCor=0,nPCs=0) {
	projection.data <- as.matrix(rca.obj$projection.data)
	tempS<-Seurat::CreateSeuratObject(as.matrix(rca.obj$raw.data))
	if (nPCs==0){
		if (bigCor!=0){
			bigM<-as.matrix(bigcor(projection.data,nblocks=bigCor,method=corMeth))
			cat(head(bigM))
			cat(dim(bigM))
			projection<-as.dist(1-bigM)
		}else{
			if (require(HiClimR) & (corMeth=="pearson")){
				projection<-as.dist(1-HiClimR::fastCor(projection.data))
			}else{
				projection<-as.dist(1-cor(projection.data,method=corMeth))
				}
			}
	}
	else{
		pca_result<-cor(t(irlba::prcomp_irlba(projection.data,n=nPCs,center=F,scale.=F)$rotation),method=corMeth)
		colnames(pca_result)<-colnames(projection.data)
		row.names(pca_result)<-colnames(projection.data)
		projection<-as.dist(1-pca_result)
	}
	cat(dim(as.matrix(projection)))
	cat(head(colnames(as.matrix(projection))))
	cat(head(rownames(as.matrix(projection))))
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
