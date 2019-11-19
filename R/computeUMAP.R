#' Compute a umap projection into 2 or 3 dimensions
#'
#' @param rca.obj RCA object
#' @param nDIMS
#' @export
#'

computeUMAP <- function(rca.obj, nDIMS=2) {

    ### Extract projection data from RCA object
    projection = as.matrix(rca.obj$projection.data)
    clusterColorList = rca.obj$clustering.out$dynamicColorsList

    ### Check if package dependencies are available; if not, download from CRAN and require those packages
    # umap
    if (!require(umap))
        install.packages("umap", repos = "http://cran.us.r-project.org")
    require(umap)

    # Compute UMAP projection from cell type projection
    if ((nDIMs < 2) | (nDIMS > 3))
	print("Error: nDIMS must be set to 2 or 3");
    	return(rca.obj)
    else{
        umap.projection <- umap(t(projection),n_components=nDIMS)
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
