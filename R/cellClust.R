#' Generate cell clusters
#' 
#' Hierarchical clustering and dynamic tree cutting. 
#' 
#' @param obj_in data object.
#' @param method can only be "hclust" (default). 
#' @param deepSplit_wgcna integer value indicating how deep the dendrogram should be cut. Values can be 0, 1 (default), 2, 3 and 4. 
#' @param min_group_Size_wgcna integer value indicating the minimum size of the resulting clusters. Default is 5.
#' @return data object.
#' @export
#' @examples
#' 
#' data_obj = cellClust(data_obj);
#' 
cellClust <- function(obj_in,method="hclust",deepSplit_wgcna=1,min_group_Size_wgcna=5)
{
  fpkm_temp = obj_in$fpkm_for_clust
  
  if (!require(fastcluster)) install.packages("fastcluster",repos = "http://cran.us.r-project.org") 
  require(fastcluster)
  if (!require(WGCNA)){
    source("http://bioconductor.org/biocLite.R")
    biocLite(c("impute", "GO.db", "preprocessCore"))
    install.packages("WGCNA")
  }
  require(WGCNA)  
  if (require(HiClimR)){
  d = as.dist(1-fastCor(fpkm_temp,upperTri=TRUE, nSplit=5, optBLAS = T))
  }else{
  d = as.dist(1-cor(fpkm_temp,method="pearson"))
  }
  cellTree = fastcluster::hClust(d,method = "average")     
  dynamicGroups = cutreeDynamic(dendro = cellTree,distM = as.matrix(d),deepSplit = deepSplit_wgcna,pamStage = FALSE,minClusterSize= min_group_Size_wgcna)
  dynamicColors = labels2colors(dynamicGroups)
    
  return(list("d" = d,"cellTree" = cellTree,"dynamicColors" = dynamicColors,"group_labels_color" = group_labels_color))
}
