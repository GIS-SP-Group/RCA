#' Normalize data
#' 
#' Normalize single cell RNAseq data. 
#' 
#' @param obj_in data object.
#' @param method can choose between "no_norm" (Default, no normalization would be performed) and "scQ" (an adjusted quantile normalization for single cell data).
#' @return data object.
#' @export
#' @examples
#' 
#' data_obj = cellNormalize(data_obj);
#' 
cellNormalize <- function(obj_in,method="no_norm")
{
  
  #### 1: reading input ####  
  
  fpkm_temp = obj_in$fpkm;
  
  #### 2: choosing method ####
  if (method == "scQ"){
    
    normed_fpkm = scQ(fpkm_temp,(min(apply(fpkm_temp,2,function(x) sum(x>0)))-5));

  }
  
  if (method == "no_norm"){
    normed_fpkm = fpkm_temp;
  }
  
  #### 3: writing output ####  
  obj_out = append(obj_in,
              list("normed_fpkm" = normed_fpkm)
              ); 
  return(obj_out)
  
}