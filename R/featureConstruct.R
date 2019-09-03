#'Helper Function
fuzzy_transfer <- function(x,min_val){
			upper_thresh = max(log10(50),median(x[x>min_val]));
			y=x;
			y[x==min_val]=0;
			y[x>=upper_thresh]=1;
			y[x<upper_thresh & x>min_val]=(x[x<upper_thresh & x>min_val]-min_val)/upper_thresh;
			return(y);
		}


#' Compute Reference Component features for clustering analysis
#' 
#' @param obj_in data object.
#' @param method: Either "GlobalPanel"(default), "ColonEpitheliumPanel", "SelfProjection", or the filename (including path) to any custom panel
#' @param corMeth: Any of the correlation measures supported by R, defaults to pearson
#' @param power: power to raise up to for the RCA features before clusterin, default is 4
#' @param scale: True if the data should be scaled, False otherwise
#' @return data object
#' @export
#' @examples
#' 
#' data_obj = featureConstruct(data_obj);
#' 
#' 
featureConstruct <- function(obj_in,method="GlobalPanel",corMeth="pearson",power=4,scale=T){
data(sysdata, envir=environment())
if (method == "GlobalPanel") {
	data11 = list()
	for (i in 1:length(sysdata[[1]])) {
		panel = sysdata[[1]][[i]]
		clustering_panel = hclust(as.dist(1-cor(panel)), method = "average")
		query_genes <- paste("", lapply(strsplit(gsub("^.*?_", "",  rownames(obj_in$fpkm_transformed)), "_ENS"), "[[",1), sep = "")
		shared_genes <- intersect(query_genes, rownames(panel))
		matched_panel_genes = rownames(obj_in$fpkm_transformed)[match(shared_genes, query_genes)]
		subset_panel = panel[shared_genes, ]
		subset_data = obj_in$fpkm_transformed[matched_panel_genes, , drop = FALSE]
		subset_panel[subset_panel <= (sysdata$at)[i]] = (sysdata$at)[i]
		cor_matrix <-cor(subset_panel,subset_data,method=corMeth)
		cor_matrix = abs(cor_matrix)^(power)*sign(cor_matrix)
		if (scale){
			cor_matrix = scale(cor_matrix, center = TRUE, scale = TRUE)
		}
	        data11[[i]] = cor_matrix
		}
	cor_matrix = as.data.frame(rbind(data11[[1]],data11[[2]]))
	}else{
		if (method == "ColonEpitheliumPanel"){
			fc = apply(sysdata$ColonEpiPanel,1,function(x) x-median(x));
			fs = fc>1.5;
			fs1 = rownames(sysdata$ColonEpiPanel[apply(fs,1,function(x) sum(x))>0,])
			gl_intersect = intersect(rownames(fpkm_temp),fs1)
			cor_matrix = as.data.frame(cor(fpkm_temp[gl_intersect,],sysdata$ColonEpiPanel[gl_intersect,],corMeth))
			cor_matrix = abs(cor_matrix)^(power)*sign(cor_matrix)
			if (scale){
				cor_matrix = scale(cor_matrix, center = TRUE, scale = TRUE)
			}
		}
		else{
			if (method == "SelfProjection"){    
				min_val = min(fpkm_temp)
				fpkm_transfer = apply(fpkm_temp,1,function(x) fuzzy_transfer(x,min_val))
				fpkm_transfer1 = as.data.frame(t(fpkm_transfer))
				ones_vec = apply(fpkm_transfer1,1,function(x) sum(x==1))
				zeros_vec = apply(fpkm_transfer1,1,function(x) sum(x==0))
				fs_vec_fuzzy = ones_vec >= dim(fpkm_temp)[2]*1/15 & zeros_vec >= dim(fpkm_temp)[2]*1/15     
				cor_matrix = cor(fpkm_temp[fs_vec_fuzzy,,drop=FALSE],method=corMeth)
				if (scale){
					cor_matrix = scale(cor_matrix, center = TRUE, scale = TRUE)
				}
			}
			else{
				panel<-readRDS(method)
				distance_panel = as.dist(1 - cor(panel))
				clustering_panel = hclust(distance_panel, method = "average")
				query_genes <- paste("", lapply(strsplit(gsub("^.*?_", "",  rownames(obj_in$fpkm_transformed)), "_ENS"), "[[",1), sep = "")
				shared_genes <- intersect(query_genes, rownames(panel))
				matched_panel_genes = rownames(obj_in$fpkm_transformed)[match(shared_genes, query_genes)]
				subset_panel = panel[shared_genes, ]
				subset_data = obj_in$fpkm_transformed[matched_panel_genes, , drop = FALSE]
				cor_matrix <-cor(subset_panel,subset_data,method=corMeth)
				cor_matrix = abs(cor_matrix)^(power)*sign(cor_matrix)
				if (scale){
					cor_matrix = scale(cor_matrix, center = TRUE, scale = TRUE)
				}
				cor_matrix = as.data.frame(cor_matrix)	
				}
			}
		}
return(append(obj_in,list("fpkm_for_clust" = cor_matrix)))
}
