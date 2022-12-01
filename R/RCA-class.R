#' RCA Class
#'
#' @importFrom methods new
#' @export
#'
RCAConstruct <- setRefClass(Class = "RCA",
			    fields = list(raw.data = "Matrix",
					  data = "Matrix",
					  projection.data = "Matrix",
					  clustering.out = "list",
					  umap.coordinates = "data.frame",
					  cell.Type.Estimate.per.cell = "list",
					  cell.Type.Estimate.per.cluster = "list",
					  baseColors = "list",
					  rRank = "list",
					  cScore = "list",
					  DE.genes = "list"))



RCAConstruct$methods(show = function(){
		     print("RCA reference class object")
			     dataSize = dim(raw.data);
			     print(paste0("Raw data: ",dataSize[2]," cells and ",dataSize[1]," features."))

			     dataSize = dim(data);
			     print(paste0("Normalized data: ",dataSize[2]," cells and ",dataSize[1]," features."))

			     dataSize = dim(projection.data);
			     print(paste0("Projection data: ",dataSize[2]," cells to ",dataSize[1]," cell-types."))

			     nDim = dim(umap.coordinates)[2]
			     print(paste0("UMAP coordinates are available for ",nDim," dimensions."))

			     dataSize = length(unique(clustering.out$dynamicColorsList[[1]]));
			     print(paste0("The data set contains ",dataSize," RCA clusters."))

			     dataSize = length(unique(cell.Type.Estimate.per.cell));
			     print(paste0("The data set contains ",dataSize," unique cell types. (per cell annotation)"))
			     
			     
			     dataSize = length(unique(cell.Type.Estimate.per.cluster));
			     print(paste0("The data set contains ",dataSize," unique cell types. (per cell annotation)"))

		})


