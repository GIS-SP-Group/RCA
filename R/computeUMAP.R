#' Compute a umap projection into 2 or 3 dimensions
#'
#' @param rca.obj RCA object
#' @param nDIMS number of UMAP dimensions
#' @param nThreads Number of threads to use for UMAP
#' @export
#'
computeUMAP <- function(rca.obj, nDIMS=2, nThreads=1) {

    ### Extract projection data from RCA object
    projection = as.matrix(rca.obj$projection.data)
    clusterColorList = rca.obj$clustering.out$dynamicColorsList
    
    # Compute UMAP projection from cell type projection
    if ((nDIMS < 2) | (nDIMS > 3)){
	print("Error: nDIMS must be set to 2 or 3")
    	return(rca.obj)}
    else{
        umap.projection <- umap::umap(t(projection),n_components=nDIMS, random_state=1,transform_state=1, n_threads = nThreads)
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
