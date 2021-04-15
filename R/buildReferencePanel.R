#' Internal function to build a reference panel from your bulk RNA sequencing data and a list of marker genes
#'
#' @param bulk.rna.data Bulk RNA sequencing data - ideally TPM normalized
#' @param de.genes list of marker genes
#' @param gene.nomenclature Type of gene nomenclature. "ENS" for Ensembl IDs, "SYMBOL" for gene symbols. Default is "SYMBOL".
#' @param species "HUMAN" for human data, "MOUSE" for mouse data. Default is "HUMAN".
#' @return reference panel boject
#'

buildPanel <- function(bulk.rna.data, de.genes, gene.nomenclature = "SYMBOL", species = "HUMAN") {

    # Get union of DE genes
    de.union <- unique(de.genes)

    # Define reference panel
    ref.panel <- bulk.rna.data[de.union, ]
    unique.ct <- unique(as.vector(sapply(colnames(bulk.rna.data),function(x){return(strsplit(x,"_")[[1]][1])})))

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
#' @param gene.nomenclature Type of gene nomenclature. "ENS" for Ensembl IDs, "SYMBOL" for gene symbols. Default is "SYMBOL".
#' @param species "HUMAN" for human data, "MOUSE" for mouse data. Default is "HUMAN".
#' @param filename file name of reference panel (should end with .rds)
#' @param verbose Generate debug output and figures (default FALSE)
#' @return reference panel as a data.frame
#' @export
#'

buildReferencePanel <- function(bulk.rna.data, fc.thrs.general = 6, fc.thrs.specific =2, fdr.thrs = 0.01, cut_height = 1.1, gene.nomenclature = "SYMBOL", species = "HUMAN", filename = "my_reference_panel.rds", verbose=FALSE) {
	#Split cell type from replicate information
	cellTypes<-as.vector(sapply(colnames(bulk.rna.data),function(x){return(strsplit(x,"_")[[1]][1])}))

	#Generating pseudo bulk across replicates
	pseudo_pseudo_bulk<-c()
	for (ct in unique(cellTypes)){
		pseudo_pseudo_bulk<-cbind(pseudo_pseudo_bulk,apply(bulk.rna.data[,which(cellTypes==ct)],1,mean))
	}
	colnames(pseudo_pseudo_bulk)<-unique(cellTypes)

	#Generate PCA of the pseudo bulk generated above
	pca_pseudo<-prcomp(t(pseudo_pseudo_bulk),scale. = T)
	totalVar<-(sum(pca_pseudo$sdev*pca_pseudo$sdev))
	varExplained<-(pca_pseudo$sdev*pca_pseudo$sdev)/totalVar
	if(verbose){
		pca_heatmap<-gplots::heatmap.2(cor(t(pca_pseudo$x[,c(1:(min(which(varExplained<0.01))-1))]),method="spearman"),col="bluered",margins=c(10,10),trace = NULL,tracecol = NULL)
		png("Heatmap_PCA_Clustering.png",width=1280,height=1280)
		pca_heatmap
		dev.off()
	}


	bulk.rna.data_cluster<-hclust(as.dist(1-cor(t(pca_pseudo$x[,c(1:(min(which(varExplained<0.01))-1))]),method="spearman")))
	bulk.rna.data_cluster_assignment<-cutree(bulk.rna.data_cluster,h=cut_height)


	#Loop through the cluster identified in the PCA to identify marker genes within those clusters using the fc.thrs.specific
	SubMat=list()
	markerGenesSpecific<-c()
	for(cluster in unique(bulk.rna.data_cluster_assignment)){
		if(verbose)
			print(cluster)
		current_CellTypes<-unique(cellTypes[which(cellTypes%in%names(bulk.rna.data_cluster_assignment)[which(bulk.rna.data_cluster_assignment==cluster)])]			      )
		if(verbose)
		    	print(current_CellTypes)

		SubMat[[cluster]]<-bulk.rna.data[,which(cellTypes%in%names(bulk.rna.data_cluster_assignment)[which(bulk.rna.data_cluster_assignment==cluster)])]
		if (length(current_CellTypes)>1){
	        	n_unique_cellTypes<-length(unique(current_CellTypes))
		        tmpGeneList<-c()
		for (i in c(1:(n_unique_cellTypes-1))){
			for (j in c((i+1):n_unique_cellTypes)){
				if(verbose){
					print(current_CellTypes[i])
				        print(current_CellTypes[j])
				}
			        selected_cellTypes<-bulk.rna.data[,union(which(cellTypes==current_CellTypes[i]),which(cellTypes==current_CellTypes[j]))]
			        dgeObj <- edgeR::DGEList(selected_cellTypes, group = as.vector(sapply(colnames(selected_cellTypes),function(x){return(strsplit(x,"_")[[1]][1])})))
			        dgeObj <- edgeR::estimateCommonDisp(dgeObj)
			        dgeObj <- edgeR::estimateTagwiseDisp(dgeObj)
			        etObj <- edgeR::exactTest(object = dgeObj)
			        etObj$table$FDR <- stats::p.adjust(p = etObj$table$PValue, method = "fdr")
				etObj$table<-etObj$table[which(etObj$table$FDR<fdr.thrs),]
				etObj$table<-etObj$table[which(abs(etObj$table$logFC)>fc.thrs.specific),]
				tmpGeneList<-unique(c(tmpGeneList,row.names(etObj$table)[order(etObj$table$logFC,decreasing = T)]))
				if(verbose)
					print(etObj$table)
				}
			}
		markerGenesSpecific<-unique(union(markerGenesSpecific,tmpGeneList))
		}
	}

	#Find general markers between all cell types with a strong fold change.
	ct<-unique(cellTypes)
	n_unique_cellTypes<-length(ct)
	tmpGeneList<-c()
	markerGenesLenient<-c()
	for (i in c(1:(n_unique_cellTypes-1))){
		for (j in c((i+1):n_unique_cellTypes)){
			if(verbose){
				print(ct[i])
			        print(ct[j])
			}
		        selected_cellTypes<-bulk.rna.data[,union(which(cellTypes==ct[i]),which(cellTypes==ct[j]))]
		        dgeObj <- edgeR::DGEList(selected_cellTypes, group = as.vector(sapply(colnames(selected_cellTypes),function(x){return(strsplit(x,"_")[[1]][1])})))
		        dgeObj <- edgeR::estimateCommonDisp(dgeObj)
		        dgeObj <- edgeR::estimateTagwiseDisp(dgeObj)
		        etObj <- edgeR::exactTest(object = dgeObj)
		        etObj$table$FDR <- stats::p.adjust(p = etObj$table$PValue, method = "fdr")
			etObj$table<-etObj$table[which(etObj$table$FDR<fdr.thrs),]
			etObj$table<-etObj$table[which(abs(etObj$table$logFC)>fc.thrs.general),]
			tmpGeneList<-unique(c(tmpGeneList,row.names(etObj$table)[order(etObj$table$logFC,decreasing = T)]))
			if(verbose)
				print(etObj$table)
			}
		}
	markerGenesLenient<-unique(union(markerGenesLenient,tmpGeneList))

	#Generate the panel considering the union of specific and general marker genes
	Panel_Selection<-buildPanel(bulk.rna.data,unique(union(markerGenesSpecific,markerGenesLenient)), gene.nomenclature = "SYMBOL", species = "HUMAN")

	if(verbose){
		heatmap_panel<-gplots::heatmap.2(cor(Panel_Selection,method="kendal"),col="bluered",margins =c(10,10))
		png("Heatmap_Reference_Panel.png",width=1280,height=1280)
		heatmap_panel
		dev.off()
	}

	saveRDS(Panel_Selection,file=filename)

	return(Panel_Selection)
}

