#' RCA Class
#'
#' @export
#'
RCAConstruct <- setRefClass(Class = "RCA", 
			    fields = list(raw.data = "Matrix", 
					  data = "Matrix", 
					  projection.data = "Matrix", 
					  clustering.out = "list",
					  umap.coordinates = "data.frame",
					  cell.Type.Estimate = "list"))

RCAConstruct$methods(show=function(){
		     print("RCA reference class object")
		     if (!(is.null(raw.data))){
			     dataSize=dim(raw.data);
			     print(paste0("Raw data: ",dataSize[2]," cells and ",dataSize[1]," features."))
			    }
                     if (!(is.null(data))){
			     dataSize=dim(data);
			     print(paste0("Filtered data: ",dataSize[2]," cells and ",dataSize[1]," features."))
		     }
		     if (!(is.null(projection.data))){
			     dataSize=dim(projection.data);
			     print(paste0("Projection data: ",dataSize[2]," cells to ",dataSize[1]," cell-types."))
		     }
		     if (!(is.null(clustering.out))){
			     dataSize=length(unique(clustering.out$dynamicColorsList[[1]]));
			     print(paste0("The data set contains ",dataSize," RCA clusters."))}
		     if (!(is.null(cell.Type.Estimate))){
			     dataSize=length(unique(cell.Type.Estimate));
			     print(paste0("The data set contains ",dataSize," unique cell types."))}
		})

		      
