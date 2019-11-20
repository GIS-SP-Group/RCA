#' Plot bar plots showing the composition of RCA clusters
#'
#' @param rca.obj data matrix (genes x cells)
#' @param folderpath path to save heatmap to
#' @param filename file name of saved heatmap
#' @export
#'

plotRCAClusterComposition <- function(rca.obj, folderpath = ".", filename = "RCA_Heatmap.pdf") {
    ### Extract projection data and clustering result from RCA object
    heatmapIn = as.matrix(rca.obj$projection.data)
    cellTree = rca.obj$clustering.out$cellTree
    clusterColorList = rca.obj$clustering.out$dynamicColorsList[[1]]

    ### Check if package dependencies are available; if not, download from CRAN and require those packages
    # dplyr Package
    if (!require(dplyr))
        install.packages("dplyr", repos = "http://cran.us.r-project.org")
    require(dplyr)

    # ggplot2 Package
    if (!require(ggplot2))
        install.packages("ggplot2", repos = "http://cran.us.r-project.org")
    require(ggplot2)

    # ggplot2 Package
    if (!require(gridExtra))
        install.packages("gridExtra", repos = "http://cran.us.r-project.org")
    require(gridExtra)

    enrichmentAll<-c()
    for(type in unique(clusterColors)){
	index=which(clusterColors==type)
        enrichmentAll<-rbind(enrichmentAll,(cbind(type,table(unlist(rca.obj$cell.Type.Estimate)[index]))))
    }
    enrichmentAll<-data.frame(cbind(row.names(enrichmentAll),enrichmentAll))
    colnames(enrichmentAll)<-c("CT","Cluster","Count")
    rownames(enrichmentAll)<-c(1:dim(enrichmentAll)[1])
    enrichmentAll$Count<-as.numeric(as.character(enrichmentAll$Count))
    totalCounts<-data.frame(count(enrichmentAll,wt=Count,Cluster))
    enrichmentAll<-left_join(enrichmentAll,totalCounts,by="Cluster")
    enrichmentAll<-cbind(enrichmentAll,enrichmentAll$Count/enrichmentAll$n*100)
    colnames(enrichmentAll)[5]<-"Ratio"
    enrichmentAll$Ratio<-as.numeric(as.character(enrichmentAll$Ratio))
    
    if (require(randomcoloR)){
       dColors<-distinctColorPalette(length(unique(enrichmentAll$CT)))
       ratioPlot<-ggplot2::ggplot(enrichmentAll,aes(x=Cluster,y=Ratio,fill=CT))+geom_bar(stat="identity")+theme_bw(15)+ylab("Percentage")+coord_flip()+ggtitle("a)")+theme(legend.position = "none")+scale_fill_manual(values=dColors)
       countPlot<-ggplot2::ggplot(enrichmentAll,aes(x=Cluster,y=Count,fill=CT))+geom_bar(stat="identity")+theme_bw(15)+ylab("Count")+coord_flip()+ggtitle("b)")+scale_fill_manual(values=dColors)+labs(fill="Cell type")+guides(fill=guide_legend(ncol=1))
    }
    else{
       ratioPlot<-ggplot2::ggplot(enrichmentAll,aes(x=Cluster,y=Ratio,fill=CT))+geom_bar(stat="identity")+theme_bw(15)+ylab("Percentage")+coord_flip()+ggtitle("a)")+theme(legend.position = "none")
       countPlot<-ggplot2::ggplot(enrichmentAll,aes(x=Cluster,y=Count,fill=CT))+geom_bar(stat="identity")+theme_bw(15)+ylab("Count")+coord_flip()+ggtitle("b)")+labs(fill="Cell type")+guides(fill=guide_legend(ncol=1))
    }
   pdf(paste0(folderpath,"Heatmap_Composition_"filename,width=15,height=7))
   grid.arrange(ratioPlot,countPlot,widths=c(1,1.4),nrow=1)
   dev.off()
}
