#' Estimate the most likely cell type from the projection to the reference panel
#'
#' @param rca.obj RCA object.
#' @param confidence: a parameter indicating the difference between z-scores. If the difference is below this threshold, the cell type will be set to unknown. Default is NULL.
#' @param ctRank: parameter indicating whether a relative rank coloring for each cell type shall be computed. Default is FALSE.
#' @param cSCompute: parameter indicating wheter the confidence score should be computed for each cell. Default is FALSE.
#' @return RCA object.
#' @export
#'
estimateCellTypeFromProjection <- function(rca.obj, confidence=NULL, ctRank=F, cSCompute=F) {

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

    cTIdfAlpha<-function(x){
        temp<-x
        tempMax<-max(temp)
        index<-which(temp==tempMax)
        temp<-temp[-index]
        deltaMax<-max(temp)/tempMax
        return(1.0-abs(deltaMax))
    }

    cTIdfConfCol<-function(x,index,bC){
        colorVec<-colorRampPalette(c("grey",bC))(100)
        maxVal<-max(x[,index])
        maxIndex<-which(x[,index]==maxVal)
        cellTypeOrder<-order(x[maxIndex,])
	ratio<-max(1,(which(cellTypeOrder==index)/length(cellTypeOrder))*100)
	result<-colorVec[ratio]
	return(result)
    }

    cellTypes<-list()
    confidenceScore<-list()
    relativeColorRank<-list()
    for (i in c(1:dim(rca.obj$projection.data)[2])){
    	  if (is.null(confidence)){
          cellTypes<-c(cellTypes,cTIdfWU(rca.obj$projection.data[,i]))
	  }
	  else{
	  cellTypes<-c(cellTypes,cTIdf(rca.obj$projection.data[,i],confidence))
	  }
    }
    if (cSCompute){
        for (i in c(1:dim(rca.obj$projection.data)[2])){
            confidenceScore<-c(confidenceScore,cTIdfAlpha(rca.obj$projection.data[,i]))
        }
	rca.obj$cScore <-confidenceScore
    }else{
    rca.obj$cScore <- list()
    }
    if (ctRank){
	require("randomcoloR")
        myColors<-distinctColorPalette(length(unique(cellTypes)))
        names(myColors)<-unique(cellTypes)
        baseColors<-myColors[unlist(cellTypes)]
	rca.obj$baseColors<-list(Colors=baseColors)
        for (i in c(1:dim(rca.obj$projection.data)[2])){
            relativeColorRank<-c(relativeColorRank,cTIdfConfCol(rca.obj$projection.data,i,baseColors[i]))
        }
	rca.obj$rRank <- relativeColorRank
    }else{
    rca.obj$rRank <- list()
    rca.obj$baseColors<-list()
    }


    # Assign projection result to RCA object
    rca.obj$cell.Type.Estimate <- cellTypes
    ### Return RCA object
    return(rca.obj)
}
