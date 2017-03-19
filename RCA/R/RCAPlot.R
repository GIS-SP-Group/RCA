#' Visualize the result of RCA clustering
#' 
#' Visualize the result of RCA via scatter (PCA) plot and heat maps.
#' 
#' @param rca_obj input data object.
#' @param point_cex scaling parameter for point size in scatter plots. Default value is 1.
#' @param cluster_color_labels label of cell identities in colors based on external information. Default value is NULL.
#' @return Two types of png format plots, including scatter plot in PCA space and heat map with columns representing samples and rows representing Reference Component Projections.
#' @export
#' @examples
#' 
#' To color cells with default RCA clusters:
#' RCAplot(data_obj); 
#' 
#' To color cells with external color labels: 
#' RCAPlot(data_obj,cluster_color_labels = color_to_use);
#' Note: the "color_to_use" vector needs to be of the length of total number of input cells, and ordered according to the columns of the input expression data frame (e.g. FPKM_raw). 
#'          
RCAPlot <- function(rca_obj,point_cex=1,cluster_color_labels=NULL){
  
  if (!require("WGCNA")) install.packages("WGCNA",repos = "http://cran.us.r-project.org") 
  if (!require("gplots")) install.packages("gplots",repos = "http://cran.us.r-project.org")
#  if (!require("colorRamp")) install.packages("colorRamp",repos = "http://cran.us.r-project.org") 
  require(gplots)
  require(WGCNA)
#  require(colorRamp)
    data(sysdata, envir=environment())
    
    c = rca_obj$dynamicColors;
    pr = prcomp(t(scale(rca_obj$fpkm_for_clust))); 
    pc_projection = as.data.frame(pr$x);
    cell_projection = pc_projection[,1:2];
    pch_to_use = 21;
    cex_to_use = point_cex;

  png("RCAplot_PCA_RCA_clusters.png")
  main_name = "PCA of cell clusters in RCA space";
  #par(mar=c(1,1,1,1)*7)  
  plot(cell_projection[,1],cell_projection[,2],
       type="p",pch = pch_to_use,
       xlab = colnames(cell_projection)[1],ylab = colnames(cell_projection)[2],
       col = c,lwd = 1,bg = c,
       main = main_name, cex = cex_to_use, cex.main = 2,
       font.axis = 1, font.lab = 1, font.main = 1
      );
  dev.off();
  
  
  if (!is.null(cluster_color_labels)){
    
  png("RCAplot_PCA_external_labels.png")
  main_name = "PCA of cell clusters in RCA space";
  #par(mar=c(1,1,1,1)*7)  
  plot(cell_projection[,1],cell_projection[,2],
       type="p",pch = pch_to_use,
       xlab = colnames(cell_projection)[1],ylab = colnames(cell_projection)[2],
       col = cluster_color_labels,lwd = 1,bg = cluster_color_labels,
       main = main_name, cex = cex_to_use, cex.main = 2,
       font.axis = 1, font.lab = 1, font.main = 1
      );
  dev.off();
  }

  png("RCAplot_heatmap_RCA_clusters.png",width = 3000, height = 3000, units = "px", pointsize = 15);
  color_scheme =     colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100);
  heatmap.2(as.matrix(rca_obj$fpkm_for_clust),col=color_scheme,
          Colv=as.dendrogram(rca_obj$cellTree),
          ColSideColors=c,
           scale="none", margins=c(5,20),
           trace="none",
           key = TRUE,
           keysize = 0.5,
           cexCol = 1,cexRow =1,
           labCol = ""
           )  
  dev.off();
  
  if (!is.null(cluster_color_labels)){  
  png("RCAplot_heatmap_external_labels.png",width = 3000, height = 3000, units = "px", pointsize = 15);
  color_scheme =     colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                       "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100);
  heatmap.2(as.matrix(rca_obj$fpkm_for_clust),col=color_scheme,
          Colv=as.dendrogram(rca_obj$cellTree),
          ColSideColors=cluster_color_labels,
           scale="none", margins=c(5,20),
           trace="none",
           key = TRUE,
           keysize = 0.5,
           cexCol = 1,cexRow =1,
           labCol = ""
           )  
  dev.off();
  }


  }
 