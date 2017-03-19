#' Filt genes
#' 
#' remove genes with too low expression levels in most cells.
#' 
#' @param obj_in data object
#' @return data object
#' @export
#' @examples
#' 
#' data_obj = geneFilt(data_obj);
#' 
geneFilt <- function(obj_in, method = "default")
{
  #### 1: reading input ####  
  fpkm_temp = obj_in$fpkm_raw;
  
  #### 2: evaluating method ####
  if (method == "default"){
    gene_filter = apply(fpkm_temp,1,function(x) sum(x>1e-3))>=2;
    fpkm_filtered = fpkm_temp[gene_filter,,drop=FALSE];
  }
  #### 3: writing output ####  
  obj_out = append(obj_in,
              list("fpkm"=fpkm_filtered,
                   "geneFilter"=gene_filter)
  )
  
  return(obj_out)
}