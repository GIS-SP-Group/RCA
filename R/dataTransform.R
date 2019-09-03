#' Data Transformation
#' 
#' Transform raw values into log10 scale
#' 
#' @param obj_in data object.
#' @param method can only be "log10" (default), indicating a log10 transformation of the data.
#' @return data object.
#' @export
#' @examples
#' 
#' data_obj = dataTransform(data_obj);
#' 
dataTransform <- function(obj_in,method = "log10")
{
	if (method == "log10"){
		fpkm_transform = log10(obj_in$normed_fpkm+1)
	}

	if (method == "log10_10x"){
		fpkm_transform = log10(obj_in$normed_fpkm+0.1)
	}  
return(append(obj_in,list("fpkm_transformed" = fpkm_transform)))
}
