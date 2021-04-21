#' Compute pairwise DE genes for supervised clustering result.
#'
#' @param rca.obj RCA object.
#' @param logFoldChange Log fold change required to call gene DE.
#' @param method Denotes which test to use. Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
#'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
#'  (McDavid et al., Bioinformatics, 2013)
#'  \item{"roc"} : Identifies 'markers' of gene expression using ROC analysis.
#'  For each gene, evaluates (using AUC) a classifier built on that gene alone,
#'  to classify between two groups of cells. An AUC value of 1 means that
#'  expression values for this gene alone can perfectly classify the two
#'  groupings (i.e. Each of the cells in cells.1 exhibit a higher level than
#'  each of the cells in cells.2). An AUC value of 0 also means there is perfect
#'  classification, but in the other direction. A value of 0.5 implies that
#'  the gene has no predictive power to classify the two groups. Returns a
#'  'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative differentially
#'  expressed genes.
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#'  }
#' @param mean.Exp Minimum mean expression of a gene to be considered in the DE gene calculation
#' @param deepsplit If hclust was used for clustering, the desired deepsplit can be specified here.. Values can range from 0 to 4. Default is 1.
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.25
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param random.seed Random seed for downsampling. default is 1
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param pseudocount.use Pseudocount to add to averaged expression values when calculating logFC. 1 by default.
#' @param p.adjust.methods correction method for calculating qvalue. default is BH (or FDR)
#' @param top.genes.per.cluster Number of top DE genes to be considered per cluster
#' @param pairwise Flag indicating whether DE genes should be compared derived in pairwise manner or 1 cluster vs all others (Default).
#' @param nCores Number of cores to used for parallel computation (default 1).

#' @return RCA object.
#' @export
#'
dataDE <- function(rca.obj,
                   logFoldChange = 1.5,
                   method = "wilcox",
                   mean.Exp = 0.5,
                   deepsplit = 1,
                   min.pct = 0.25,
                   min.diff.pct = -Inf,
                   random.seed = 1,
                   min.cells.group = 3,
                   pseudocount.use = 1,
                   p.adjust.methods =  "BH",
                   top.genes.per.cluster = 10,
		           pairwise=FALSE,
		           nCores=1) {
    clusteri <- group1 <- group2 <- avg_logFC <- Cluster <- NULL
    `%dopar%` <- foreach::`%dopar%`
    `%>%` <- magrittr::`%>%`
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)
    df <- base::c()
    if (base::is.null(mean.Exp)){
    temp.exp = base::exp(x = rca.obj$data)
    temp.exp.row = Matrix::rowMeans(temp.exp)
    temp.exp.row = base::sort(temp.exp.row, decreasing = T)
    temp.exp.row = temp.exp.row[6:base::length(temp.exp.row)]
    MeanExprsThrs = base::mean(temp.exp.row)}
    else{
	    MeanExprsThrs = mean.Exp
    }
    ############################
    #hclust used for clustering#
    ############################
    if (base::class(rca.obj$clustering.out$cellTree) == "hclust") {
        clusters <- rca.obj$clustering.out$dynamicColorsList[[deepsplit]]
        total.clus <- base::length(base::unique(clusters))
        remap <- c(1:total.clus)
        base::names(remap) <-
            base::unique(rca.obj$clustering.out$dynamicColorsList[[deepsplit]])
        clusters <- base::as.numeric(base::as.character(remap[clusters]))
    } else{
    #############################
    #Graph based clustering used#
    #############################
        clusters <- rca.obj$clustering.out$dynamicColorsList[[1]]
        total.clus <- base::length(base::unique(clusters))
        remap <- base::c(1:total.clus)
        base::names(remap) <- base::unique(rca.obj$clustering.out$dynamicColorsList[[1]])
        clusters <- base::as.numeric(base::as.character(remap[clusters]))
    }
    ###########################
    #Compute pairwise DE genes#
    ###########################
    if (pairwise) {
      df <- foreach::foreach(clusteri = 1:(total.clus - 1), .combine = rbind) %dopar% {
	 tmp <- base::c()
       	 for (clusterj in (clusteri + 1):(total.clus)) {
            cells.1 <- base::colnames(rca.obj$data)[base::which(clusters == clusteri)]
            cells.2 <- base::colnames(rca.obj$data)[base::which(clusters == clusterj)]
            marker.genes = ComputePairWiseDE(
                object = rca.obj$data,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = NULL,
                logfc.threshold = logFoldChange,
                test.use = method,
                min.pct = min.pct,
                min.diff.pct = min.diff.pct,
                verbose = TRUE,
                only.pos = FALSE,
                max.cells.per.ident = Inf,
                random.seed = random.seed,
                min.cells.group = min.cells.group,
                pseudocount.use = pseudocount.use,
                MeanExprsThrs = MeanExprsThrs,
                p.adjust.methods = p.adjust.methods)
            if (!(base::is.null(marker.genes))) {
                if (base::colnames(marker.genes)[1] != 'myAUC') {
                    marker.genes = marker.genes[marker.genes$p_val_adj < 0.05, ]
                }
                if (base::nrow(marker.genes) > 0) {
                    marker.genes$group1 = clusteri
                    marker.genes$group2 = clusterj
                    marker.genes$gene = base::rownames(marker.genes)
                    tmp <- base::rbind(tmp, marker.genes)
                }
            }
        }
	tmp
      }
    }else{
      df <- foreach::foreach(clusteri = 1:(total.clus),.combine = rbind) %dopar% {
            cells.1 <- base::colnames(rca.obj$data)[base::which(clusters == clusteri)]
            cells.2 <- base::colnames(rca.obj$data)[base::which(clusters != clusteri)]
            marker.genes = ComputePairWiseDE(
                object = rca.obj$data,
                cells.1 = cells.1,
                cells.2 = cells.2,
                features = NULL,
                logfc.threshold = logFoldChange,
                test.use = method,
                min.pct = min.pct,
                min.diff.pct = min.diff.pct,
                verbose = TRUE,
                only.pos = FALSE,
                max.cells.per.ident = Inf,
                random.seed = random.seed,
                min.cells.group = min.cells.group,
                pseudocount.use = pseudocount.use,
                MeanExprsThrs = MeanExprsThrs,
                p.adjust.methods = p.adjust.methods
    		)
            if (!(base::is.null(marker.genes))) {
                if (base::colnames(marker.genes)[1] != 'myAUC') {
                    marker.genes = marker.genes[marker.genes$p_val_adj < 0.05, ]
                }
                if (base::nrow(marker.genes) > 0) {
                    marker.genes$group1 = clusteri
                    marker.genes$gene = base::rownames(marker.genes)
#                    df = rbind(df, marker.genes)
                }
		marker.genes
            }
	}
    }

    #Determine top x DE genes per Cluster#
    ######################################
    if (pairwise) {

    mC1 <- df %>% dplyr::group_by(group1, group2) %>% dplyr::top_n(n = top.genes.per.cluster, wt = avg_logFC)
    markers2 <- df
    markers2$avg_logFC <- (-1) * (markers2$avg_logFC)
    mC2 <- markers2 %>% dplyr::group_by(group2, group1) %>% dplyr::top_n(n = top.genes.per.cluster, wt = avg_logFC)
    topMarkers <- base::data.frame(base::rbind(
            base::cbind(Cluster = mC1$group1, Gene = mC1$gene),
            base::cbind(Cluster = mC2$group2, Gene = mC2$gene)
        ))
    topMarkers <- topMarkers %>% dplyr::group_by(Cluster) %>% dplyr::distinct(.keep_all = T)
    topMarkers <- topMarkers[base::order(topMarkers$Cluster, decreasing = F), ]
    df$group1<-base::names(remap)[df$group1]
    df$group2<-base::names(remap)[df$group2]
    }else{
    	topMarkers <- base::data.frame(Cluster = df$group1, Gene = df$gene, avg_logFC = df$avg_logFC)
    	topMarkers <- topMarkers %>% dplyr::group_by(Cluster) %>% dplyr::top_n(n = top.genes.per.cluster, wt = avg_logFC) %>% dplyr::distinct(.keep_all = T)
    	topMarkers <- topMarkers[base::order(topMarkers$Cluster, decreasing = F), ]
	df$group1<-base::names(remap)[df$group1]
    }
    rca.obj$DE.genes <- base::list(All.DE.genes = df, Top.DE.genes = topMarkers)
    return(rca.obj)
}


#'@author Ignasius Joanito (Modified from Seurat FindMarkers)

#' @title Find markers (differentially expressed genes) between two group of cells.
#'
#'@param object dataMatrix of genes (rows) x cells (columns) expression matrix (log normalized value)
#'@param cells.1 Vector of cell names belonging to group 1
#'@param cells.2 Vector of cell names belonging to group 2
#'@param features Genes to test. Default is NULL which mean to use all genes
#'@param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is log(1.5)
#'@param test.use Denotes which test to use. Available options are:
#' \itemize{
#'  \item{"wilcox"} : Identifies differentially expressed genes between two
#'  groups of cells using a Wilcoxon Rank Sum test (default)
#'  \item{"bimod"} : Likelihood-ratio test for single cell gene expression,
#'  (McDavid et al., Bioinformatics, 2013)
#'  \item{"roc"} : Identifies 'markers' of gene expression using ROC analysis.
#'  For each gene, evaluates (using AUC) a classifier built on that gene alone,
#'  to classify between two groups of cells. An AUC value of 1 means that
#'  expression values for this gene alone can perfectly classify the two
#'  groupings (i.e. Each of the cells in cells.1 exhibit a higher level than
#'  each of the cells in cells.2). An AUC value of 0 also means there is perfect
#'  classification, but in the other direction. A value of 0.5 implies that
#'  the gene has no predictive power to classify the two groups. Returns a
#'  'predictive power' (abs(AUC-0.5) * 2) ranked matrix of putative differentially
#'  expressed genes.
#'  \item{"t"} : Identify differentially expressed genes between two groups of
#'  cells using the Student's t-test.
#'  }
#'@param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.25
#'@param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#'@param only.pos Only return positive markers (FALSE by default)
#'@param verbose Print a progress bar once expression testing begins
#'@param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf)
#'@param random.seed Random seed for downsampling. default is 1
#'@param min.cells.group Minimum number of cells in one of the groups
#'@param pseudocount.use Pseudocount to add to averaged expression values when calculating logFC. 1 by default.
#'@param MeanExprsThrs a minimum expression threshold of average cluster expression for a gene to be considered a DE gene.
#' the mean expression value is in the linear scale!
#'@param p.adjust.methods correction method for calculating qvalue. default is BH (or FDR)
#'@return returnObj
#list CompGeneList: List of genes for which the DE statistical test was performed for each pairwise cluster comparison
#list qValueList: list of q-values from the DE statistical test for genes where test was performed for each pairwise cluster comparison
#list logFCList: list of log-fold-changes from the DE statistical test for genes where test was performed for each pairwise cluster comparison=
#character vector deGeneUnion: union of top ndeg DE genes from up- and down- regulated set for all pairwise comparison
#numeric matrix deCountMatrix: matrix with number of DE genes for each pairwise cluster comparison
#list deGeneRegulationList: list of up and down regulated DE genes for each pairwise cluster comparison ordered by log-fold-change
#list log2FCDEList: list of log2-fold-change for up and down regulated DE genes for each pairwise cluster comparison - corresponds to same order as in deGeneRegulationList
#list qValueDEList: list of q-values for up and down regulated DE genes for each pairwise cluster comparison - corresponds to same order as in deGeneRegulationList
#list upregulatedDEGeneList: list of cluster-specific upregulated DE genes
#list downregulatedDEGeneList: list of cluster-specific downregulated DE genes

#' @export
ComputePairWiseDE <-  function(object,
                               cells.1 = NULL,
                               cells.2 = NULL,
                               features = NULL,
                               logfc.threshold = 1.5,
                               test.use = "wilcox",
                               min.pct = 0.25,
                               min.diff.pct = -Inf,
                               verbose = TRUE,
                               only.pos = FALSE,
                               max.cells.per.ident = Inf,
                               random.seed = 1,
                               min.cells.group = 3,
                               pseudocount.use = 1,
                               MeanExprsThrs = 0,
                               p.adjust.methods = "BH") {
    `%||%` <- rlang::`%||%`
    ## for Wilcox test
    WilcoxDETest <-
        function(data.use, cells.1, cells.2, verbose = TRUE) {
            group.info <- base::data.frame(row.names = base::c(cells.1, cells.2))
            group.info[cells.1, "group"] <- "Group1"
            group.info[cells.2, "group"] <- "Group2"
            group.info[, "group"] <- base::factor(x = group.info[, "group"])
            data.use <- data.use[, base::rownames(x = group.info), drop = FALSE]
            my.sapply <- base::ifelse(
                test = verbose && future::nbrOfWorkers() == 1,
                yes = pbapply::pbsapply,
                no = future.apply::future_sapply
            )
            p_val <- my.sapply(
                X = 1:nrow(x = data.use),
                FUN = function(x) {
                    return(stats::wilcox.test(data.use[x,] ~ group.info[, "group"])$p.value)
                }
            )
            return(base::data.frame(p_val, row.names = base::rownames(x = data.use)))
        }

    ## List of DE Functions
    ## for bimod
    # Likelihood ratio test for zero-inflated data
    # Identifies differentially expressed genes between two groups of cells using
    # the LRT model proposed in McDavid et al, Bioinformatics, 2013
    bimodLikData <- function(x, xmin = 0) {
        x1 <- x[x <= xmin]
        x2 <- x[x > xmin]
        xal <- Seurat::MinMax(
            data = base::length(x = x2) / base::length(x = x),
            min = 1e-5,
            max = (1 - 1e-5)
        )
        likA <- base::length(x = x1) * base::log(x = 1 - xal)
        if (base::length(x = x2) < 2) {
            mysd <- 1
        } else {
            mysd <- stats::sd(x = x2)
        }
        likB <- base::length(x = x2) *
            base::log(x = xal) +
            base::sum(stats::dnorm(
                x = x2,
                mean = base::mean(x = x2),
                sd = mysd,
                log = TRUE
            ))
        return(likA + likB)
    }
    DifferentialLRT <- function(x, y, xmin = 0) {
        lrtX <- bimodLikData(x = x)
        lrtY <- bimodLikData(x = y)
        lrtZ <- bimodLikData(x = base::c(x, y))
        lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
        return(stats::pchisq(
            q = lrt_diff,
            df = 3,
            lower.tail = F
        ))
    }
    DiffExpTest <-
        function(data.use, cells.1,  cells.2, verbose = TRUE) {
            my.sapply <- base::ifelse(
                test = verbose && future::nbrOfWorkers() == 1,
                yes = pbapply::pbsapply,
                no = future.apply::future_sapply
            )
            p_val <- base::unlist(x = my.sapply(
                X = 1:base::nrow(x = data.use),
                FUN = function(x) {
                    return(DifferentialLRT(
                        x = base::as.numeric(x = data.use[x, cells.1]),
                        y = base::as.numeric(x = data.use[x, cells.2])
                    ))
                }
            ))
            to.return <-
                base::data.frame(p_val, row.names = base::rownames(x = data.use))
            return(to.return)
        }
    ## for ROC test
    DifferentialAUC <- function(x, y) {
        prediction.use <- ROCR::prediction(
            predictions = base::c(x, y),
            labels = base::c(base::rep(x = 1, base::length(x = x)), base::rep(x = 0, base::length(x = y))),
            label.ordering = 0:1
        )
        perf.use <-
            ROCR::performance(prediction.obj = prediction.use, measure = "auc")
        auc.use <- base::round(x = perf.use@y.values[[1]], digits = 3)
        return(auc.use)
    }
    AUCMarkerTest <-
        function(data1, data2, mygenes, print.bar = TRUE) {
            myAUC <- base::unlist(x = base::lapply(
                X = mygenes,
                FUN = function(x) {
                    return(DifferentialAUC(x = base::as.numeric(x = data1[x,]),
                                           y = base::as.numeric(x = data2[x,])))
                }
            ))
            myAUC[base::is.na(x = myAUC)] <- 0
            iterate.fxn <-
                base::ifelse(test = print.bar,
                       yes = pbapply::pblapply,
                       no = base::lapply)
            avg_diff <- base::unlist(x = iterate.fxn(
                X = mygenes,
                FUN = function(x) {
                    return(Seurat::ExpMean(x = base::as.numeric(x = data1[x,])) - Seurat::ExpMean(x = base::as.numeric(x = data2[x,])))
                }
            ))
            toRet <- base::data.frame(base::cbind(myAUC, avg_diff), row.names = mygenes)
            toRet <- toRet[base::rev(x = base::order(toRet$myAUC)),]
            return(toRet)
        }
    MarkerTest <-
        function(data.use, cells.1, cells.2, verbose = TRUE) {
            to.return <- AUCMarkerTest(
                data1 = data.use[, cells.1, drop = FALSE],
                data2 = data.use[, cells.2, drop = FALSE],
                mygenes = base::rownames(x = data.use),
                print.bar = verbose
            )
            to.return$power <- base::abs(x = to.return$myAUC - 0.5) * 2
            return(to.return)
        }
    # Differential expression testing using Student's t-test
    # Identify differentially expressed genes between two groups of cells using the Student's t-test
    DiffTTest <-
        function(data.use, cells.1, cells.2, verbose = TRUE) {
            my.sapply <- base::ifelse(
                test = verbose && future::nbrOfWorkers() == 1,
                yes = pbapply::pbsapply,
                no = future.apply::future_sapply
            )
            p_val <- base::unlist(x = my.sapply(
                X = 1:base::nrow(data.use),
                FUN = function(x) {
                    stats::t.test(x = data.use[x, cells.1], y = data.use[x, cells.2])$p.value
                }
            ))
            to.return <-
                base::data.frame(p_val, row.names = base::rownames(x = data.use))
            return(to.return)
        }


    features <- features %||% base::rownames(x = object)
    # error checking
    if (base::length(x = cells.1) == 0) {
        stop("Cell group 1 is empty - identity of group 1 need to be defined ")
    }
    else if (base::length(x = cells.2) == 0) {
        stop("Cell group 2 is empty - identity of group 2 need to be defined ")
        return(NULL)
    }
    else if (base::length(x = cells.1) < min.cells.group) {
        base::print("The considered cluster for group 1 consists of too few cells")
        return(NULL)
    }
    else if (base::length(x = cells.2) < min.cells.group) {
        base::print("The considered cluster for group 2 consists of too few cells")
        return(NULL)
    }
    else if (base::any(!cells.1 %in% base::colnames(x = object))) {
        bad.cells <-
            base::colnames(x = object)[base::which(x = !base::as.character(x = cells.1) %in% base::colnames(x = object))]
        stop(
            "The following cell names provided to cells.1 are not present: ",
            base::paste(bad.cells, collapse = ", ")
        )
    } else if (base::any(!cells.2 %in% base::colnames(x = object))) {
        bad.cells <-
            base::colnames(x = object)[base::which(x = !base::as.character(x = cells.2) %in% base::colnames(x = object))]
        stop(
            "The following cell names provided to cells.2 are not present: ",
            base::paste(bad.cells, collapse = ", ")
        )
    }

    # feature selection (based on percentages)
    thresh.min <- 0
    pct.1 <- base::round(
        x = Matrix::rowSums(x = object[features, cells.1, drop = FALSE] > thresh.min) /
            base::length(x = cells.1),
        digits = 16
    )
    pct.2 <- base::round(
        x = Matrix::rowSums(x = object[features, cells.2, drop = FALSE] > thresh.min) /
            base::length(x = cells.2),
        digits = 16
    )
    object.alpha <- base::cbind(pct.1, pct.2)
    base::colnames(x = object.alpha) <- c("pct.1", "pct.2")
    alpha.min <- base::apply(X = object.alpha,
                       MARGIN = 1,
                       FUN = max)
    base::names(x = alpha.min) <- base::rownames(x = object.alpha)
    features <- base::names(x = base::which(x = alpha.min > min.pct))
    if (base::length(x = features) == 0) {
        base::print("No features pass min.pct threshold")
        return(NULL)
    }
    alpha.diff <-
        alpha.min - base::apply(X = object.alpha,
                          MARGIN = 1,
                          FUN = min)
    features <- base::names(x = base::which(x = alpha.min > min.pct &
                                    alpha.diff > min.diff.pct))
    if (base::length(x = features) == 0) {
        base::print("No features pass min.diff.pct threshold")
        return(NULL)
    }
    #
    # feature selection (based on average difference)
    mean.fxn <- function(x) {
        return(base::log(x = base::mean(x = base::expm1(x = x)) + pseudocount.use))
    }
    object.1 <- base::apply(X = object[features, cells.1, drop = FALSE],
                      MARGIN = 1,
                      FUN = mean.fxn)
    object.2 <- base::apply(X = object[features, cells.2, drop = FALSE],
                      MARGIN = 1,
                      FUN = mean.fxn)
    total.diff <- (object.1 - object.2)

    # feature selection (based on mean exprssion threshold)
    features = base::names(x = base::which(
        x = base::expm1(object.1) > MeanExprsThrs  |
            base::expm1(object.2) > MeanExprsThrs
    ))
    if (base::length(x = features) == 0) {
        base::print("No features pass log mean exprssion threshold")
        return(NULL)
    }

    # feature selection (based on logfc threshold)
    features.diff <- if (only.pos) {
        base::names(x = base::which(x = total.diff > logfc.threshold))
    } else {
        base::names(x = base::which(x = base::abs(x = total.diff) > logfc.threshold))
    }
    features <- base::intersect(x = features, y = features.diff)
    if (base::length(x = features) == 0) {
        base::print("No features pass logfc.threshold threshold")
        return(NULL)
    }

    # sampling cell for DE computation
    if (max.cells.per.ident < Inf) {
        base::set.seed(seed = random.seed)
        # Should be cells.1 and cells.2?
        if (base::length(x = cells.1) > max.cells.per.ident) {
            cells.1 <- base::sample(x = cells.1, size = max.cells.per.ident)
        }
        if (base::length(x = cells.2) > max.cells.per.ident) {
            cells.2 <- base::sample(x = cells.2, size = max.cells.per.ident)
        }
    }

    # perform DE
    de.results <- base::switch(
        EXPR = test.use,
        'wilcox' = WilcoxDETest(
            data.use = object[features, base::c(cells.1, cells.2), drop = FALSE],
            cells.1 = cells.1,
            cells.2 = cells.2,
            verbose = verbose
        ),
        'bimod' = DiffExpTest(
            data.use = object[features, base::c(cells.1, cells.2), drop = FALSE],
            cells.1 = cells.1,
            cells.2 = cells.2,
            verbose = verbose
        ),
        'roc' = MarkerTest(
            data.use = object[features, base::c(cells.1, cells.2), drop = FALSE],
            cells.1 = cells.1,
            cells.2 = cells.2,
            verbose = verbose
        ),
        't' = DiffTTest(
            data.use = object[features, base::c(cells.1, cells.2), drop = FALSE],
            cells.1 = cells.1,
            cells.2 = cells.2,
            verbose = verbose
        ),
        stop("Unknown test: ", test.use)
    )
    diff.col <- "avg_logFC"
    de.results[, diff.col] <- total.diff[base::rownames(x = de.results)]
    de.results <-
        base::cbind(de.results, object.alpha[base::rownames(x = de.results), , drop = FALSE])

    if (only.pos) {
        de.results <- de.results[de.results[, diff.col] > 0, , drop = FALSE]
    }

    if (test.use == "roc") {
        de.results <-
            de.results[base::order(-de.results$power,-de.results[, diff.col]),]
    } else {
        de.results <-
            de.results[base::order(de.results$p_val,-de.results[, diff.col]),]
        de.results$p_val_adj = stats::p.adjust(p = de.results$p_val,
                                        method = p.adjust.methods)
    }

    return(de.results)
}
