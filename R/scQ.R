#' scQ (single cell quantile normalization) for scRNAseq data
#' @param data Expression data frame with genes in rows. 
#' @param thr cells with detected genes less than \code{thr+1} would be  thrown away. Also after normalization
#' only \code{thr} number of genes would have non zero epression in each cell.
#' @return normalized data as a data frame retaining orginal row and column names
#' @examples 
#' 
#' norm_Data <- scQ(fpkm_data,1000) # scQ normalization 

scQ <- function(data,thr)
{
  # loading required libraries
  if (!require("preprocessCore")) 
  {
    source("http://bioconductor.org/biocLite.R")
    biocLite("preprocessCore")
  }
  require(preprocessCore)
  FPKMbackup <- as.matrix(data)
  # keeping genes that are expressed in at least one cell
  FPKMbackup1 <- apply(FPKMbackup>0,1,function(x) sum(x))>=1
  FPKMbackup1<-subset(FPKMbackup,FPKMbackup1)
  # finding cells which have thr+1 genes
  bu = apply(FPKMbackup1>0,2,function(x) sum(x>0)) >= thr+1
  FPKMreduced<-subset(FPKMbackup1,select=bu)
  #class(FPKMreduced)
  kk<-apply(FPKMreduced,2, function(x) sort(as.numeric(x),decreasing = TRUE)[thr+1])
  for(i in 1:dim(FPKMreduced)[2])
  {
    #i<-3
    ind<-order(FPKMreduced[,i],decreasing=TRUE)[(thr+1):length(FPKMreduced[,i])]
    FPKMreduced[,i][ind]<-kk[i]
  }
  storage.mode(FPKMreduced)<-"double"
  X <- normalize.quantiles.robust(as.matrix(FPKMreduced),use.median=TRUE,use.log2=FALSE)
  #X<-normalize.quantiles(as.matrix(FPKMreduced))
  X<-data.frame(X)
  rownames(X)<-rownames(FPKMreduced)
  colnames(X)<-colnames(FPKMreduced)
  # dividing the whole matrix by  the minimum value of the matrix 
  # This will make the pseudo count 1
  X <- data.frame(X)/min(X)
  B<-apply(X,1,function(x) sd(x))>0
  X<-X[B,]
  #View(head(X))
  return(X)
  
}
