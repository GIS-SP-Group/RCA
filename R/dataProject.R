#' Compute Reference Component features for clustering analysis
#'
#' @param rca.obj RCA object.
#' @param method: Either "GlobalPanel"(default), "ColonEpitheliumPanel", or "Custom"
#' @param customPath: directory path (including filename) to any custom panel stored in RDS format. Only used if method == "Custom".
#' @param corMeth: Any of the correlation measures supported by R, defaults to pearson
#' @param power: power to raise up to for the RCA features before clustering, default is 4
#' @param scale: True if the data should be scaled, False otherwise
#' @return RCA object.
#' @export
#' @examples
#'
#' data_projection = dataProject(data)
#'
dataProject <- function(rca.obj, method = "GlobalPanel", customPath = NULL, corMeth = "pearson", power = 4, scale = T) {
    
    # Extract data
    sc_data <- rca.obj$data
    
    # dplyr
    if (!require(dplyr))
        install.packages("dplyr", repos = "http://cran.us.r-project.org")
    require(dplyr)
    
    # Load reference panel data from environment
    data(sysdata, envir = environment())
    
    # If panel for correlation is GlobalPanel
    if (method == "GlobalPanel") {
        
        # Initialise variable to store projection data from the two fragments of the Global Panel
        projection_list = list()
        
        # For each fragment of the Global Panel
        for (i in 1:length(sysdata[[1]])) {
            
            # Initialise panel 
            panel = sysdata[[1]][[i]]
            
            # Select genes that are shared by the input data and the panel 
            shared_genes <- intersect(rownames(sc_data), rownames(panel))
            
            # Reduce the panel and input data to the shared genes
            subset_panel = panel[shared_genes, ]
            subset_data = sc_data[shared_genes, , drop = FALSE]
            
            # For values in the panel below the minimum threshold, set those values to threshold
            subset_panel[subset_panel <= (sysdata$at)[i]] = (sysdata$at)[i]
            
            # Compute projection of input data with the panel fragment
            projection_fragment <- cor(subset_panel, subset_data, method = corMeth)
            
            # Raise the projection fragment to power
            projection_fragment = abs(projection_fragment) ^ (power) * sign(projection_fragment)
            
            # If scaling is required
            if (scale) {
                
                # Scale
                projection_fragment = scale(projection_fragment, center = TRUE, scale = TRUE)
            }
            
            # Store projection data of fragment of Global Panel
            projection_list[[i]] = projection_fragment
        }
        
        # Combine the projection result of multiple Global Panel fragments
        projection = do.call("rbind", projection_list)
        
    }
    # If panel for correlation is ColonEpitheliumPanel
    else if (method == "ColonEpitheliumPanel") {
        
        # Scale panel by median
        fc = apply(sysdata$ColonEpiPanel, 1, function(x) x - median(x))
        
        fs = fc > 1.5
        
        fs1 = rownames(sysdata$ColonEpiPanel[apply(fs, 1, function(x)
            sum(x)) > 0,])
        gl_intersect = intersect(rownames(fpkm_temp), fs1)
        projection = as.data.frame(cor(fpkm_temp[gl_intersect,], sysdata$ColonEpiPanel[gl_intersect,], corMeth))
        projection = abs(projection) ^ (power) * sign(projection)
        if (scale) {
            projection = scale(projection,
                               center = TRUE,
                               scale = TRUE)
        }
    }
    # If no provided method is chosen, it is assumed that the user wishes to use a custom panel
    else {
        
        # Load panel from path provided
        panel <- readRDS(customPath)
        
        # Select genes that are shared by the input data and the panel
        shared_genes <- intersect(rownames(sc_data), rownames(panel))
        
        # Reduce the panel and input data to the shared genes
        subset_panel = panel[shared_genes, ]
        subset_data = sc_data[shared_genes, , drop = FALSE]
        
        # Compute projection of input data with the panel
        projection <- cor(subset_panel, subset_data, method = corMeth)
        
        # Raise the projection to power
        projection = abs(projection) ^ (power) * sign(projection)
        
        # If scaling is required
        if (scale) {
            
            # Scale
            projection = scale(projection,
                               center = TRUE,
                               scale = TRUE)
        }
        
        # Store projection result of custom panel as data frame
        projection = as.data.frame(projection)
    }
    
    # Assign projection result to RCA object
    rca.obj$projection.data <- projection
    
    ### Return RCA object
    
    return(rca.obj)
}
