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
    
  #### 1: reading input ####  

  fpkm_temp = obj_in$normed_fpkm;
  
  #### 2: choosing method ####

  if (method == "log10"){
    pseudo_count = 1;
    log_fpkm_temp = fpkm_temp;
    log_fpkm_temp[log_fpkm_temp<=pseudo_count] = pseudo_count;
    log_fpkm_temp = log10(log_fpkm_temp);
    fpkm_transform = log_fpkm_temp;
  }

  #### 3: writing output ####  
  obj_out = append(obj_in,
              list("fpkm_transformed" = fpkm_transform
                   )
              )
  
  return(obj_out)
  
}