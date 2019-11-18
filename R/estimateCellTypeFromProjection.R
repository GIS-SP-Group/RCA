#' Estimate the most likely cell type from the projection to the reference panel
#'
#' @param rca.obj RCA object.
#' @param confidence: a parameter indicating the difference between z-scores. If the difference is below this threshold, the cell type will be set to unknown. Default is NULL.
#' @return RCA object.
#' @export
#'
estimateCellTypeFromProjection <- function(rca.obj, confidence=NULL) {

    projection <- rca.obj$projection.data

    cTIdf<-function(x,confidence){
	  temp<-x
          tempMax<-max(temp)
          index<-which(temp==tempMax)
          temp<-temp[-index]
          deltaMax<-max(temp)/tempMax
          if (deltaMax < confidence)
            return(names(x)[index])
          else return("Unkown")
    }

    cTIdfWU<-function(x){
         return(names(x)[which(x==max(x))])
    }



    cellTypes<-c()
    for (i in c(1:dim(rca.obj$projection.data)[2])){
    	  if (confidence == NULL){
          cellTypesWU<-c(cellTypesWU,cTIdfWU(rca.obj$projection.data[,i]))
	  }
	  else{
	  cellTypes<-c(cellTypes,cTIdf(rca.obj$projection.data[,i],confidence))
	  }
    }

    # Assign projection result to RCA object
    rca.obj$cell.Type.Estimate <- cellTypes

    ### Return RCA object

    return(rca.obj)
}
