#' Function to build a reference panel from your bulk RNA sequencing data
#' 
#' @param bulk.rna.data Bulk RNA sequencing data - ideally TPM normalized.
#' @param celltype.vec Vector of cell type names of each bulk RNA sequencing sample.
#' @param fc.thrs fold change threshold
#' @param fdr.thrs FDR threshold
#' @param gene.nomenclature Type of gene nomenclature. "ENS" for Ensembl IDs, "SYMBOL" for gene symbols. Default is "SYMBOL".
#' @param species "HUMAN" for human data, "MOUSE" for mouse data. Default is "HUMAN".
#' @param folder.path path to store reference panel
#' @param filename file name of reference panel (should end with .rds)
#' 
#' @export
#' 

buildReferencePanel <- function(bulk.rna.data, celltype.vec = colnames(bulk.rna.data), fc.thrs = 7.5, fdr.thrs = 1e-02, gene.nomenclature = "SYMBOL", species = "HUMAN", folder.path = ".", filename = "my_reference_panel.rds") {
    
    if(ncol(bulk.rna.data) != length(celltype.vec)) {
        stop("Data doesn't match cell type vector.")
    }
    
    # Name bulk RNA replicates using cell type vector
    celltype.vec <- gsub(pattern = "_", replacement = "\\.", x = celltype.vec)
    
    unique.ct <- unique(celltype.vec)
    ct.list <- lapply(unique.ct, function(ct){0})
    names(ct.list) <- unique.ct
    for(i in 1:length(celltype.vec)) {
        ct.list[[celltype.vec[i]]] <- ct.list[[celltype.vec[i]]] + 1
        
        colnames(bulk.rna.data)[i] <- paste0(celltype.vec[i], "_", ct.list[[celltype.vec[i]]])
    }
    
    library(edgeR)
    
    # Convert to an edgeR object
    dgeObj <- DGEList(bulk.rna.data, group = celltype.vec)
    
    
    ## Identify genes with at least 1 cpm in at least 2 samples
    cpm.bulk.rna.data <- cpm(dgeObj)
    thresh <- cpm.bulk.rna.data > 1
    keep <- rowSums(thresh) >= 2
    
    # Filter genes based on expression
    dgeObj <- dgeObj[keep, ]
    
    # Estimate common and tag-wise dispersion for the pair of clusters
    dgeObj <- estimateCommonDisp(dgeObj)
    dgeObj <- estimateTagwiseDisp(dgeObj)
    
    # Calculate norm factors for data
    dgeObj <- calcNormFactors(dgeObj, method = "TMM")
    
    etList <- list()
    for(i in 1:(length(unique.ct)-1)) {
        for(j in (i+1):length(unique.ct)) {
            
            # Perform exact test
            etObj <- exactTest(object = dgeObj, pair = c(i, j))
            etObj$table$FDR <- p.adjust(p = etObj$table$PValue, method = "fdr")
            print(paste0(i, ", ", j))
            
            etList[paste0(i, "_", j)] <- etObj
        }
    }
    
    # Get list of DE genes
    de.list <- lapply(etList, function(etTable) {
        rownames(etTable)[which((abs(etTable$logFC) > fc.thrs) & (etTable$FDR < fdr.thrs))]
    })
    
    # Get union of DE genes
    de.union <- unique(unlist(de.list))
    
    # Define reference panel
    ref.panel <- cpm.bulk.rna.data[de.union, ]
    
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
    
    # For Ensembl gene IDs
    if(gene.nomenclature == "ENS") {
        ensembl_genes <- grep(pattern = "^ENS", x = rownames(average.ref.panel), value = T)
        
        # if Ensembl gene version numbers exist, remove them
        ensembl_genes <- gsub(pattern = "\\..*$", replacement = "", x = ensembl_genes)
        
        if(length(ensembl_genes) != length(unique(ensembl_genes))) {
            stop("Non-unique gene names in input data.")
        }
        
        rownames(average.ref.panel) <- ensembl_genes
        
        if(species == "HUMAN") {
            gene.id.df <- ensembldb::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys= ensembl_genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
        } else if(species == "MOUSE") {
            gene.id.df <- ensembldb::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79, keys= ensembl_genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
        }
        
        average.ref.panel <- average.ref.panel[gene.id.df$GENEID, ]
        
        # Removing duplicated gene symbols
        average.ref.panel <- average.ref.panel[-which(duplicated(gene.id.df$SYMBOL)), ]
        rownames(average.ref.panel) <- gene.id.df$SYMBOL[-which(duplicated(gene.id.df$SYMBOL))]
        
    }
    
    # Save reference panel
    saveRDS(object = average.ref.panel, file = paste0(folder.path, "/", filename))
    
    # Return object
    return(average.ref.panel)
}