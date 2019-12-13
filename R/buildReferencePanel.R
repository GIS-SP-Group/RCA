#' Function to build a reference panel from your bulk RNA sequencing data
#' 
#' @param bulk.rna.data Bulk RNA sequencing data - ideally TPM normalized.
#' @param celltype.vec Vector of cell type names of each bulk RNA sequencing sample.
#' @param folder.path path to store reference panel
#' @param filename file name of reference panel (should end with .rds)
#' 
#' @export
#' 

buildReferencePanel <- function(bulk.rna.data, celltype.vec, folder.path = ".", filename = "my_reference_panel.rds") {
    
    if(ncol(bulk.rna.data) != length(celltype.vec)) {
        stop("Data doesn't match cell type vector.")
    }
    
    ## Name bulk RNA replicates using cell type vector
    celltype.vec <- gsub(pattern = "_", replacement = "\\.", x = celltype.vec)
    
    unique.ct <- unique(celltype.vec)
    ct.list <- lapply(unique.ct, function(ct){0})
    names(ct.list) <- unique.ct
    for(i in 1:length(celltype.vec)) {
        ct.list[[celltype.vec[i]]] <- ct.list[[celltype.vec[i]]] + 1
        
        colnames(bulk.rna.data)[i] <- paste0(celltype.vec[i], "_", ct.list[[celltype.vec[i]]])
    }
    
    
    ## Identify genes with at least 0.5 cpm in at least 2 samples
    thresh <- bulk.rna.data > 0.5
    keep <- rowSums(thresh) >= 2
    
    # Subset the rows of countdata to keep the more highly expressed genes
    filt.bulk.rna.data <- bulk.rna.data[keep,]
    
    library(edgeR)
    
    ## Convert to an edgeR object
    dgeObj <- DGEList(filt.bulk.rna.data, group = celltype.vec)
    
    # Estimate common and tag-wise dispersion for the pair of clusters
    dgeObj <- estimateCommonDisp(dgeObj)
    dgeObj <- estimateTagwiseDisp(dgeObj)
    
    # Calculate norm factors for data
    dgeObj <- calcNormFactors(dgeObj, method = "none")
    
    etList <- list()
    for(i in 1:(length(unique.ct)-1)) {
        for(j in (i+1):length(unique.ct)) {
            
            # Perform exact test
            etObj <- exactTest(object = dgeObj, pair = c(i, j))
            etObj$table$FDR <- p.adjust(p = etObj$table$PValue, method = "fdr")
            
            etList[paste0(i, "_", j)] <- etObj
        }
    }
    
    # Get list of DE genes
    de.list <- lapply(etList, function(etTable) {
        rownames(etTable)[which((abs(etTable$logFC) > 7.5) & (etTable$FDR < 0.0001))]
    })
    
    # Get union of DE genes
    de.union <- unique(unlist(de.list))
    
    # Define reference panel
    ref.panel <- filt.bulk.rna.data[de.union, ]
    
    # Average out replicates
    average.ref.panel <- matrix(0, nrow(ref.panel), ncol = length(unique.ct))
    rownames(average.ref.panel) <- rownames(ref.panel)
    colnames(average.ref.panel) <- unique.ct
    for(ct in unique.ct) {
        if(sum(grepl(pattern = paste0("^", ct, "_"), x = colnames(ref.panel))) > 1) {
            average.ref.panel[, ct] <- rowMeans(ref.panel[, grep(pattern = paste0("^", ct, "_"), x = colnames(ref.panel), value = T)])
        } else {
            average.ref.panel[, ct] <- ref.panel[, grep(pattern = paste0("^", ct, "_"), x = colnames(ref.panel), value = T)]
        }
        
    }
    
    # Save reference panel
    saveRDS(object = average.ref.panel, file = paste0(folder.path, "/", filename))
}