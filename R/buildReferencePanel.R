#' Internal function to build a reference panel from your bulk RNA sequencing data and a list of marker genes
#' Bulk RNA sequencing data must contain gene symbols - we recommend conversion prior to use this function, e.g. using Biomart.
#'
#' @param bulk.rna.data Bulk RNA sequencing data - ideally TPM normalized
#' @param de.genes list of marker genes
#' @param species "HUMAN" for human data, "MOUSE" for mouse data. Default is "HUMAN".
#' @return reference panel boject
#'

buildPanel <- function(bulk.rna.data, de.genes,  species = "HUMAN") {

    # Get union of DE genes
    de.union <- base::unique(de.genes)

    # Define reference panel
    ref.panel <- bulk.rna.data[de.union, ]
    unique.ct <- base::unique(base::as.vector(base::sapply(base::colnames(bulk.rna.data),function(x){return(base::strsplit(x,"_")[[1]][1])})))

    # Average out replicates
    average.ref.panel <- base::matrix(0, base::nrow(ref.panel), ncol = base::length(unique.ct))
    base::rownames(average.ref.panel) <- base::rownames(ref.panel)
    base::colnames(average.ref.panel) <- unique.ct
    for(ct in unique.ct) {
        if(base::sum(base::grepl(pattern = base::paste0("^", ct, "_"), x = base::colnames(ref.panel))) > 1) {
            average.ref.panel[, ct] <- base::rowMeans(ref.panel[, base::grep(pattern = base::paste0("^", ct, "_"), x = base::colnames(ref.panel), value = T)])
        } else {
            average.ref.panel[, ct] <- ref.panel[, base::grep(pattern = base::paste0("^", ct, "_"), x = base::colnames(ref.panel), value = T)]
        }

    }

    # Return object
    return(average.ref.panel)
}


#' Function to build a reference panel.
#' Celltype and replicate information must be included in the colnames of the bulk.rna.data.object, separted by an underscore ("_"), e.g.: Bnaive_R1.
#'
#' @param bulk.rna.data Bulk RNA sequencing data
#' @param fc.thrs.general absolulte log fold change threshold applied for coarse marker detection (default 6)
#' @param fc.thrs.specific absolulte log fold change threshold applied for detailed marker detection (default 2)
#' @param fdr.thrs FDR threshold (default 0.01)
#' @param cut_height to cut the clustering of the reference data sets (default 1.1)
#' @param species "HUMAN" for human data, "MOUSE" for mouse data. Default is "HUMAN".
#' @param filename file name of reference panel (should end with .rds)
#' @param verbose Generate debug output and figures (default FALSE)
#' @return reference panel as a data.frame
#' @export
#'

buildReferencePanel <- function(bulk.rna.data, fc.thrs.general = 6, fc.thrs.specific =2, fdr.thrs = 0.01, cut_height = 1.1, species = "HUMAN", filename = "my_reference_panel.rds", verbose=FALSE) {
	#Split cell type from replicate information
	cellTypes<-base::as.vector(base::sapply(base::colnames(bulk.rna.data),function(x){return(base::strsplit(x,"_")[[1]][1])}))

	#Generating pseudo bulk across replicates
	pseudo_pseudo_bulk<-base::c()
	for (ct in base::unique(cellTypes)){
		pseudo_pseudo_bulk<-base::cbind(pseudo_pseudo_bulk,base::apply(bulk.rna.data[,base::which(cellTypes==ct)],1,mean))
	}
	base::colnames(pseudo_pseudo_bulk)<-base::unique(cellTypes)

	#Generate PCA of the pseudo bulk generated above
	pca_pseudo<-stats::prcomp(base::t(pseudo_pseudo_bulk),scale. = T)
	totalVar<-(base::sum(pca_pseudo$sdev*pca_pseudo$sdev))
	varExplained<-(pca_pseudo$sdev*pca_pseudo$sdev)/totalVar
	if(verbose){
		pca_heatmap<-gplots::heatmap.2(stats::cor(base::t(pca_pseudo$x[,base::c(1:(base::min(base::which(varExplained < 0.01)) - 1))]), method="spearman"), col="bluered", margins=base::c(10,10), trace = NULL, tracecol = NULL)
		grDevices::png("Heatmap_PCA_Clustering.png", width=1280, height=1280)
		pca_heatmap
		grDevices::dev.off()
	}


	bulk.rna.data_cluster <- stats::hclust(stats::as.dist(1 - stats::cor(base::t(pca_pseudo$x[,base::c(1:(base::min(base::which(varExplained < 0.01)) - 1))]),method="spearman")))
	bulk.rna.data_cluster_assignment<-stats::cutree(bulk.rna.data_cluster,h=cut_height)


	#Loop through the cluster identified in the PCA to identify marker genes within those clusters using the fc.thrs.specific
	SubMat = base::list()
	markerGenesSpecific <- base::c()
	for(cluster in base::unique(bulk.rna.data_cluster_assignment)){
		if(verbose)
		    base::print(cluster)
		current_CellTypes <- base::unique(cellTypes[base::which(cellTypes %in% base::names(bulk.rna.data_cluster_assignment)[base::which(bulk.rna.data_cluster_assignment == cluster)])]			      )
		if(verbose)
		    base::print(current_CellTypes)

		SubMat[[cluster]] <- bulk.rna.data[, base::which(cellTypes %in% base::names(bulk.rna.data_cluster_assignment)[base::which(bulk.rna.data_cluster_assignment == cluster)])]
		if (base::length(current_CellTypes) > 1){
	        	n_unique_cellTypes <- base::length(base::unique(current_CellTypes))
		        tmpGeneList <- base::c()
		for (i in base::c(1:(n_unique_cellTypes - 1))){
			for (j in base::c((i + 1):n_unique_cellTypes)){
				if(verbose){
				    base::print(current_CellTypes[i])
				    base::print(current_CellTypes[j])
				}
			        selected_cellTypes <- bulk.rna.data[, base::union(base::which(cellTypes == current_CellTypes[i]), base::which(cellTypes == current_CellTypes[j]))]
			        dgeObj <- edgeR::DGEList(selected_cellTypes, group = base::as.vector(base::sapply(base::colnames(selected_cellTypes), function(x){return(base::strsplit(x, "_")[[1]][1])})))
			        dgeObj <- edgeR::estimateCommonDisp(dgeObj)
			        dgeObj <- edgeR::estimateTagwiseDisp(dgeObj)
			        etObj <- edgeR::exactTest(object = dgeObj)
			        etObj$table$FDR <- stats::p.adjust(p = etObj$table$PValue, method = "fdr")
				etObj$table <- etObj$table[base::which(etObj$table$FDR < fdr.thrs),]
				etObj$table <- etObj$table[base::which(base::abs(etObj$table$logFC) > fc.thrs.specific),]
				tmpGeneList <- base::unique(c(tmpGeneList, base::row.names(etObj$table)[base::order(etObj$table$logFC,decreasing = T)]))
				if(verbose)
				    base::print(etObj$table)
				}
			}
		markerGenesSpecific <- base::unique(base::union(markerGenesSpecific, tmpGeneList))
		}
	}

	#Find general markers between all cell types with a strong fold change.
	ct <- base::unique(cellTypes)
	n_unique_cellTypes <- base::length(ct)
	tmpGeneList <- base::c()
	markerGenesLenient <- base::c()
	for (i in base::c(1:(n_unique_cellTypes - 1))){
		for (j in base::c((i + 1):n_unique_cellTypes)){
			if(verbose){
			    base::print(ct[i])
			    base::print(ct[j])
			}
		        selected_cellTypes<-bulk.rna.data[, base::union(base::which(cellTypes == ct[i]), base::which(cellTypes == ct[j]))]
		        dgeObj <- edgeR::DGEList(selected_cellTypes, group = base::as.vector(base::sapply(base::colnames(selected_cellTypes), function(x){return(base::strsplit(x, "_")[[1]][1])})))
		        dgeObj <- edgeR::estimateCommonDisp(dgeObj)
		        dgeObj <- edgeR::estimateTagwiseDisp(dgeObj)
		        etObj <- edgeR::exactTest(object = dgeObj)
		        etObj$table$FDR <- stats::p.adjust(p = etObj$table$PValue, method = "fdr")
			etObj$table <- etObj$table[base::which(etObj$table$FDR<fdr.thrs),]
			etObj$table <- etObj$table[base::which(base::abs(etObj$table$logFC) > fc.thrs.general),]
			tmpGeneList <- base::unique(base::c(tmpGeneList,base::row.names(etObj$table)[base::order(etObj$table$logFC, decreasing = T)]))
			if(verbose)
			    base::print(etObj$table)
			}
		}
	markerGenesLenient<- base::unique(base::union(markerGenesLenient, tmpGeneList))

	#Generate the panel considering the union of specific and general marker genes
	Panel_Selection <- buildPanel(bulk.rna.data, base::unique(base::union(markerGenesSpecific, markerGenesLenient)), species = "HUMAN")

	if(verbose){
		heatmap_panel <- gplots::heatmap.2(stats::cor(Panel_Selection, method=  "kendal"), col = "bluered", margins = base::c(10,10))
		grDevices::png("Heatmap_Reference_Panel.png", width = 1280, height = 1280)
		heatmap_panel
		grDevices::dev.off()
	}

	base::saveRDS(Panel_Selection, file = filename)

	return(Panel_Selection)
}

