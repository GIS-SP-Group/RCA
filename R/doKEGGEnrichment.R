#' Perform KEGG enrichment analysis using cluster profiler.
#'
#' @param rca.obj An RCA2 object
#' @param annotation An ontology dataset that can be obtained from bioconductor. Used for ID matching.
#' @param org Supported organism in KEGG. An overview is listed in 'http://www.genome.jp/kegg/catalog/org_list.html'. Default: hsa (human)
#' @param key one of "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot', Default: KEGG.
#' @param p.Val p-value cut-off, default 0.05
#' @param q.Val q-value cut-off, default 0.2
#' @param p.Adjust.Method p-value adjustment method to be used, default BH
#' @param gene.label.type Type of gene.labels used, default SYMBOL
#' @param filename postfix of the plots generated, GoEnrichment.pdf
#' @param background.set.threshold, minimum expression threshold used for mean gene-expression. Either a numerical value or one of the following thresholds computed on the mean gene-expression values across all genes: Min, 1stQ, Mean, Median, 3rdQ. Default: 1stQ
#' @param n.Cells.Expressed Alternative threshold to filter genes. Keep only genes that are expressed in at least n.Cells.Expressed cells
#' @param cluster.ID ID of a cluster for which the GO enrichment should be computed. If this is not provided, enrichment will be computed for all clusters. Default: NULL
#' @param deep.split Deep.split to be used if hierachical clustering was used to cluster the projection
#' @return RCA object.
#' @export
#'
doEnrichKEGG<-function(rca.obj,
		     annotation=NULL,
		     org="hsa"
		     key="kegg",
		     p.Val=0.05,
		     q.Val=0.2,
		     p.Adjust.Method="BH",
		     gene.label.type="SYMBOL",
		     filename="KEGG_Enrichment.pdf",
		     background.set.threshold="1stQ",
		     n.Cells.Expressed=NULL,
		     cluster.ID=NULL,
		     deep.split=NULL){
	#check annotation provided
	if (is.null(annotation)){
		print("No annotation provided. Download from bioconductor, e.g. org.Hs.eg.db for homo sapiens") 
		stop()
	}
  
	#check background.se.thresholds
	if (!(is.null(background.set.threshold)) & !(is.null(n.Cells.Expressed))){
		print("Only one threshold can be used.")
		stop()
	}
    
	#Check type of clustering
	if (class(rca.obj$clustering.out)!="hclust"){
		deep.split=1
	}else{
		if (is.null(deep.split)){
			print("Please specify the desired cluster split to be used")
			stop()
		}
	}
      
	#map cluster colors to numbers#
	clusters<-unique(rca.obj$clustering.out$dynamicColorsList[[deep.split]])
	map<-c(1:length(clusters))
        names(map)<-clusters
	#Determine clusters to be subjected to go test
 	if(is.null(cluster.ID)){
		allClusters<-map[unique(rca.obj$clustering.out$dynamicColorsList[[1]])]
	}else{
		allClusters<-map[cluster.ID]
	    }
	###Loop through all clusters
	for (cluster in allClusters){
		if (is.null(cluster.ID)){
			print(paste0("Performing KEGG enrichment for cluster ",names(allClusters)[cluster]))
		} else {
			print(paste0("Performing KEGG enrihment for cluster ",names(allClusters)[1]))
			}
		clusterGenes<-as.character(rca.obj$DE.genes$Top.DE.genes$Gene[which(rca.obj$DE.genes$Top.DE.genes$Cluster==cluster)])
		clusterGenes<-str_to_upper(clusterGenes)
		###Generate background set using a threshold based on the mean expression of genes across all cells.
		if (!(is.null(background.set.threshold))){
			if(is.null(cluster.ID)){
				clusterMeanExp<-apply(rca.obj$data[,which(rca.obj$clustering.out$dynamicColorsList[[1]]==names(allClusters)[cluster])],1,mean)
			}else{
				clusterMeanExp<-apply(rca.obj$data[,which(rca.obj$clustering.out$dynamicColorsList[[1]]==names(allClusters)[1])],1,mean)
			}
			if (is.numeric(background.set.threshold)){
				backgroundGeneNames<-row.names(rca.obj$data)[which(clusterMeanExp>background.set.threshold)]
			}else{
				labels<-c("Min","1stQ","Median","Mean","3rdQ")
				index<-which(labels == background.set.threshold)
				if (isEmpty(index)){
					print("The provided value for the background.set.threshold is not valid")
					stop()
					}
				backgroundGeneNames<-row.names(rca.obj$data)[which(clusterMeanExp>summary(clusterMeanExp)[index])]
				}
		}
		     
		###Generate background set using a threshold based on the mean expression of genes across all cells.
		if (!(is.null(n.Cells.Expressed))){
			geneExpVec <- Matrix::rowSums(rca.obj$raw.data>0)
			backgroundGeneNames<-row.names(rca.obj$data)[which(geneExpVec > min.cell.exp)]
			}
		###Relabel genes in case they are not in ENTREZ gene ID format already
		if (gene.label.type != "ENTREZID"){
			backgroundENTREZ<-clusterProfiler::bitr(backgroundGeneNames, fromType = gene.label.type, toType = c("ENTREZID"), OrgDb = annotation)
			gene.df <- clusterProfiler::bitr(clusterGenes, fromType = gene.label.type, toType = c("ENTREZID"), OrgDb = annotation)
		} else {
			backgroundENTREZ <- data.frame(ENTREZID=backgroundGeneNames)
			gene.df <- data.frame(ENTREZID=clusterGenes)
			}
		###Perfom actual enrichemt test
		ggo <- clusterProfiler::enrichKEGG(gene = gene.df$ENTREZID, 
						 organism = org, 
						 keyType = key, 
						 pAdjustMethod = p.Adjust.Method, 
						 pvalueCutoff  = p.Val, 
						 qvalueCutoff  = q.Val, 
						 readable = TRUE, 
						 universe = backgroundENTREZ$ENTREZID)
		
		if (!is.null(ggo)){
			if (dim(as.data.frame(ggo))[1] != 0){
				if (is.null(cluster.ID)){
				#Generate and save barplot
				ggplot2::ggsave(paste0("barplot_",names(allClusters)[cluster],"_",filename),barplot(ggo),width=15,height=8,units="in")
				#Generate and save dotplot
				ggplot2::ggsave(paste0("Dotplot_",names(allClusters)[cluster],"_",filename),clusterProfiler::dotplot(ggo),width=15,height=8,units="in")
				} else {
				#Generate and save barplot
				ggplot2::ggsave(paste0("barplot_",names(allClusters)[1],"_",filename),barplot(ggo),width=15,height=8,units="in")
				#Generate and save dotplot
				ggplot2::ggsave(paste0("Dotplot_",names(allClusters)[1],"_",filename),clusterProfiler::dotplot(ggo),width=15,height=8,units="in")
		                }
		        }
	      }
      }  
}
