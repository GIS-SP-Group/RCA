#' RCA Class
#'
#' @export
#'
RCAConstruct <- setRefClass(Class = "RCA", 
			    fields = list(raw.data = "Matrix", 
					  data = "Matrix", 
					  projection.data = "Matrix", 
					  clustering.out = "list"))

RCAConstruct$methods(show=function(){dataSize=dim(raw.data);
		     print(paste0("RCA reference class object (version 2.0) holding ",dataSize[2]," cells and ",dataSize[1],"features."))})
