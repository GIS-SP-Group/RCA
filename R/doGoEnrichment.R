#' Perform GO enrichment analysis using cluster profiler.
#'
#' @param rca.obj An RCA2 object
#' @param annotation An ontology dataset that can be obtained from bioconductor
#' @param ontology Either MF, BP, or CC, default BP
#' @param p.Val p-value cut-off, default 0.05
#' @param q.Val q-value cut-off, default 0.2
#' @param p.Adjust.Method p-value adjustment method to be used, default BH
#' @param gene.label.type Type of gene.labels used, default SYMBOL
#' @param filename postfix of the plots generated, GoEnrichment.pdf
#' @param background.set, ALL indicates that all genes are considered, CLUSTER indicates that only genes with the cluster of interest are considered (default: ALL).
#' @param background.set.threshold, minimum expression threshold used for mean gene-expression. Either a numerical value or one of the following thresholds computed on the mean gene-expression values across all genes within a considered cluster: Min, 1stQ, Mean, Median, 3rdQ. Default: NULL
#' @param n.Cells.Expressed Alternative threshold to filter genes. Keep only genes that are expressed in at least n.Cells.Expressed cells
#' @param cluster.ID ID of a cluster for which the GO enrichment should be computed. If this is not provided, enrichment will be computed for all clusters. Default: NULL
#' @param deep.split Deep.split to be used if hierachical clustering was used to cluster the projection
#' @return RCA object.
#' @export
#'
doEnrichGo<-function(rca.obj,
		     annotation=NULL,
		     ontology="BP",
		     p.Val=0.05,
		     q.Val=0.2,
		     p.Adjust.Method="BH",
		     gene.label.type="SYMBOL",
		     filename="GoEnrichment.pdf",
		     background.set="ALL",
		     background.set.threshold=NULL,
		     n.Cells.Expressed=NULL,
		     cluster.ID=NULL,
		     deep.split=NULL){
	#check annotation provided
	if (base::is.null(annotation)){
	    base::print("No annotation provided. Download from bioconductor, e.g. org.Hs.eg.db for homo sapiens")
		stop()
	} else {
	    base::print(utils::citation(annotation))
	}

	#check background.se.thresholds
	if (!(base::is.null(background.set.threshold)) & !(base::is.null(n.Cells.Expressed))){
	    base::print("Only one threshold can be used.")
		stop()
	}

	#Check type of clustering
	if (base::class(rca.obj$clustering.out)!="hclust"){
		deep.split=1
	}else{
		if (base::is.null(deep.split)){
		    base::print("Please specify the desired cluster split to be used")
			stop()
		}
	}

	#map cluster colors to numbers#
	clusters<-base::unique(rca.obj$clustering.out$dynamicColorsList[[deep.split]])
	map<-base::c(1:base::length(clusters))
	base::names(map)<-clusters
	#Determine clusters to be subjected to go test
 	if(base::is.null(cluster.ID)){
		allClusters<-map
	}else{
		allClusters<-map[cluster.ID]
	    }
	###Loop through all clusters
	for (cluster in allClusters){
		if (base::is.null(cluster.ID)){
		    base::print(base::paste0("Performing Go enrichment for cluster ",base::names(allClusters)[cluster]))
		} else {
		    base::print(base::paste0("Performing Go enrichment for cluster ",base::names(allClusters)[1]))
			}
		clusterGenes<-base::as.character(rca.obj$DE.genes$Top.DE.genes$Gene[base::which(rca.obj$DE.genes$Top.DE.genes$Cluster==cluster)])
		clusterGenes<-stringr::str_to_upper(clusterGenes)
		###Generate background set using a threshold based on the mean expression of genes across all cells.
		backgroundGeneNames=base::row.names(rca.obj$data)
		if (!(base::is.null(background.set.threshold))){
			if (background.set=="CLUSTER"){
				if(base::is.null(cluster.ID)){
					clusterMeanExp<-base::apply(rca.obj$data[,base::which(rca.obj$clustering.out$dynamicColorsList[[1]]==base::names(allClusters)[cluster])],1,mean)
				}else{
					clusterMeanExp<-base::apply(rca.obj$data[,base::which(rca.obj$clustering.out$dynamicColorsList[[1]]==base::names(allClusters)[1])],1,mean)
				}
			}else{
				clusterMeanExp<-base::apply(rca.obj$data,1,mean)
			}

			if (base::is.numeric(background.set.threshold)){
					backgroundGeneNames<-base::row.names(rca.obj$data)[base::which(clusterMeanExp>background.set.threshold)]
			}else{
				labels<-base::c("Min","1stQ","Median","Mean","3rdQ")
				index<-base::which(labels == background.set.threshold)
				if (S4Vectors::isEmpty(index)){
				    base::print("The provided value for the background.set.threshold is not valid")
					stop()
					}
				backgroundGeneNames<-base::row.names(rca.obj$data)[base::which(clusterMeanExp>base::summary(clusterMeanExp)[index])]
			}
		}


		###Generate background set using a threshold based on the mean expression of genes across all cells.
		if (!(base::is.null(n.Cells.Expressed))){
			geneExpVec <- Matrix::rowSums(rca.obj$raw.data>0)
			backgroundGeneNames<-base::row.names(rca.obj$data)[base::which(geneExpVec > n.Cells.Expressed)]
			}
		###Relabel genes in case they are not in ENTREZ gene ID format already
		if (gene.label.type != "ENTREZID"){
			backgroundENTREZ<-clusterProfiler::bitr(backgroundGeneNames, fromType = gene.label.type, toType = base::c("ENTREZID"), OrgDb = annotation)
			gene.df <- clusterProfiler::bitr(clusterGenes, fromType = gene.label.type, toType = base::c("ENTREZID"), OrgDb = annotation)
		} else {
			backgroundENTREZ <- base::data.frame(ENTREZID=backgroundGeneNames)
			gene.df <- base::data.frame(ENTREZID=clusterGenes)
			}
		###Perfom actual enrichemt test
		ggo <- clusterProfiler::enrichGO(gene = gene.df$ENTREZID,
						 OrgDb = annotation,
						 ont = ontology,
						 pAdjustMethod = p.Adjust.Method,
						 pvalueCutoff  = p.Val,
						 qvalueCutoff  = q.Val,
						 readable = TRUE,
						 universe = backgroundENTREZ$ENTREZID)

		if (!base::is.null(ggo)){
			if (base::dim(base::as.data.frame(ggo))[1] != 0){
				if (base::is.null(cluster.ID)){
				#Generate and save barplot
				ggplot2::ggsave(base::paste0("barplot_",base::names(allClusters)[cluster],"_",ontology,"_",filename),
				                graphics::barplot(ggo),width=15,height=8,units="in")
				#Generate and save dotplot
				ggplot2::ggsave(base::paste0("Dotplot_",base::names(allClusters)[cluster],"_",ontology,"_",filename),
				                clusterProfiler::dotplot(ggo),width=15,height=8,units="in")
				#Generate and save goplot
				ggplot2::ggsave(base::paste0("gotplot_",base::names(allClusters)[cluster],"_",ontology,"_",filename),
				                clusterProfiler::goplot(ggo),width=15,height=8,units="in")
				} else {
				#Generate and save barplot
				ggplot2::ggsave(base::paste0("barplot_",base::names(allClusters)[1],"_",ontology,"_",filename),
				                graphics::barplot(ggo),width=15,height=8,units="in")
				#Generate and save dotplot
				ggplot2::ggsave(base::paste0("Dotplot_",base::names(allClusters)[1],"_",ontology,"_",filename),
				                clusterProfiler::dotplot(ggo),width=15,height=8,units="in")
				#Generate and save goplot
				ggplot2::ggsave(base::paste0("gotplot_",base::names(allClusters)[1],"_",ontology,"_",filename),
				                clusterProfiler::goplot(ggo),width=25,height=8,units="in")
		                }
		        }
	      }
      }
}
