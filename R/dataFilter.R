#' @objective To filter the dataset and remove genes and cells that could be of poor quality using user-defined thresholds
#'
#' @param data data matrix (genes x cells)
#' @param nGene.thresholds numeric vector with lower and upper nGene thresholds
#' @param nUMI.thresholds numeric vector with lower and upper nUMI thresholds
#' @param percent.mito.thresholds numeric vector with lower and upper pMito thresholds
#' @param min.cell.exp minimum number of cells a gene must be expressed in
#' @return filtered data matrix
#' @export
#' @examples
#'
#' data_obj = dataFilter(data, nGene.thresholds = c(100, Inf), nUMI.thresholds = c(1000, Inf), percent.mito.thresholds = c(0.0, 0.2), min.cell.exp = 10);
#'
dataFilter <- function(data, nGene.thresholds = c(100, NULL), nUMI.thresholds = c(1000, NULL), percent.mito.thresholds = c(0.0, 0.2), min.cell.exp = 0.001*ncol(data)) {
    
    ### Gene filter ###
    
    if(!is.null(min.cell.exp)) {
        
        # Select genes with expression in a minimum number of cells
        geneExpVec <- rowSums(data>0)
        filt.genes <- rownames(data)[which(geneExpVec > min.cell.exp)]
        
    } else {
        filt.genes <- rownames(data)
    }
    
    ### Cell filter - nGene ###
    
    # Assign filtered data
    filt.cells <- colnames(data)
    
    if(!is.null(nGene.thresholds)) {
        
        # Compute nGene vector
        nGeneVec <- colSums(data>0)
        
        if(is.numeric(nGene.thresholds) & (length(nGene.thresholds) == 2)) {
            
            nGene.filt.cells <- filt.cells[which((nGeneVec >= nGene.thresholds[1]) & (nGeneVec <= nGene.thresholds[2]))]
        } else if(is.numeric(nGene.thresholds) & (length(nGene.thresholds) == 1)) {
            nGene.filt.cells <- filt.cells[which(nGeneVec >= nGene.thresholds)]
        } else {
            warning("nGene.thresholds was not of the appropriate format. Please enter a numeric vector with lower and upper thresholds.")
            nGene.filt.cells <- filt.cells
        }
    } else {
        warning("No nGene.thresholds provided.")
        nGene.filt.cells <- filt.cells
    }
    
    ### Cell filter - nUMI ###
    
    if(!is.null(nUMI.thresholds)) {
        
        # Compute nUMI vector
        nUMIVec <- colSums(data)
        
        if(is.numeric(nUMI.thresholds) & (length(nUMI.thresholds) == 2)) {
            
            nUMI.filt.cells <- filt.cells[which((nUMIVec >= nUMI.thresholds[1]) & (nUMIVec <= nUMI.thresholds[2]))]
        } else if(is.numeric(nGene.thresholds) & (length(nUMI.thresholds) == 1)) {
            nUMI.filt.cells <- filt.cells[which(nUMIVec >= nUMI.thresholds)]
        } else {
            warning("nUMI.thresholds was not of the appropriate format. Please enter a numeric vector with lower and upper thresholds.")
            nUMI.filt.cells <- filt.cells
        }
    } else {
        warning("No nUMI.thresholds provided.")
        nUMI.filt.cells <- filt.cells
    }
    
    ### Cell filter - Mitochondrial Rate ###
    
    if(!is.null(percent.mito.thresholds)) {
        
        # Select mito genes
        mito.genes = grep(pattern = "^MT-", x = rownames(data), value = T)
        
        # Compute percent.mito vector
        pMitoVec <- colSums(data[mito.genes, ])/colSums(data)
        
        
        if(is.numeric(percent.mito.thresholds) & (length(percent.mito.thresholds) == 2)) {
            
            pMito.filt.cells <- filt.cells[which((pMitoVec >= percent.mito.thresholds[1]) & (pMitoVec <= percent.mito.thresholds[2]))]
        } else if(is.numeric(percent.mito.thresholds) & (length(percent.mito.thresholds) == 1)) {
            pMito.filt.cells <- filt.cells[, which(pMitoVec >= percent.mito.thresholds)]
        } else {
            warning("percent.mito.thresholds was not of the appropriate format. Please enter a numeric vector with lower and upper thresholds.")
            pMito.filt.cells <- filt.cells
            
        }
    } else {
        warning("No percent.mito.thresholds provided.")
        pMito.filt.cells <- filt.cells
    }
    
    # Select only the cells that satisfy all 3 cell filtering criteria
    filt.cells <- intersect(intersect(nGene.filt.cells, nUMI.filt.cells), pMito.filt.cells)
    
    # Combine results from cell and gene filtering
    filt.data <- data[filt.genes, filt.cells]
    
    # Return
    return(filt.data)
}