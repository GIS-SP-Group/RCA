#' Compute a umap projection into 2 or 3 dimensions
#'
#' @param rca.obj RCA object
#' @param nDIMS number of UMAP dimensions
#' @param nThreads Number of threads to use for UMAP
#' @param nThreadsSGD Number of threads used for stochastic gradient descent (WARNING Results are not reproducible if this is > 1)
#' @param approxPow approximate the power function in stochastic gradient descent (Default FALSE)
#' @export
#'
computeUMAP <- function(rca.obj, nDIMS=2, nThreads=1, nThreadsSGD=1, approxPow=FALSE) {

    ### Extract projection data from RCA object
    projection = as.matrix(rca.obj$projection.data)
    clusterColorList = rca.obj$clustering.out$dynamicColorsList
    
    # Compute UMAP projection from cell type projection
    if ((nDIMS < 2) | (nDIMS > 3)){
	print("Error: nDIMS must be set to 2 or 3")
    	return(rca.obj)}
    else{
        umap.projection <- umap::umap(t(projection),n_components=nDIMS, random_state=1,transform_state=1, n_threads = nThreads, n_sgd_threads = nThreadsSGD, approx_pow = approxPow)
    }

    # Store UMAP layout in data frame for plotting
    umap.df <- as.data.frame(umap.projection$layout)
    if (nDIMS==2){
        colnames(umap.df) <- c("UMAP1", "UMAP2")
    }
    else{
        colnames(umap.df) <- c("UMAP1", "UMAP2","UMAP3")
    }
    rca.obj$umap.coordinates <- umap.df
    return(rca.obj)
}
