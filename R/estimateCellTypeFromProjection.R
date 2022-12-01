#' Estimate the most likely cell type from the projection to the reference panel
#'
#' @param rca.obj RCA object.
#' @param confidence a parameter indicating the difference between z-scores. If the difference is below this threshold, the cell type will be set to unknown (default NULL).
#' @param ctRank parameter indicating whether a relative rank coloring for each cell type shall be computed (default FALSE).
#' @param cSCompute parameter indicating wheter the confidence score should be computed for each cell (default FALSE).
#' @return RCA object.
#'
#' @examples
#' RCA.pbmcs <- createRCAObject(RCAv2::pbmc_small_counts)
#' RCA.pbmcs <- dataLogNormalise(RCA.pbmcs)
#' RCA.pbmcs <- dataProject(RCA.pbmcs, method = "GlobalPanel_CellTypes")
#' RCA.pbmcs <- dataClust(RCA.pbmcs)
#' RCA.pbmcs <- estimateCellTypeFromProjection(RCA.pbmcs)
#' @export
#'
estimateCellTypeFromProjection <- function(rca.obj, confidence = NULL, ctRank = F, cSCompute = F) {

        # Returns the likeliest cell type of a cell with respect to a confidence threshold
        cTIdf <- function(x, confidence) {
            temp <- x
            tempMax <- base::max(temp)
            index <- base::which(temp == tempMax)
            temp <- temp[-index]
            deltaMax <- base::max(temp) / tempMax
            if (deltaMax < confidence)
                return(base::names(x)[index])
            else
                return("Unkown")
        }

        # Returns the likeliest cell type of a cell neglecting any confidence value.
        cTIdfWU <- function(x) {
            return(base::names(x)[base::which(x == base::max(x))])
        }

        # Returns a alpha value for each cell, depending on the confidence score for the cell's cell type annotation among all possible cell types.
        cTIdfAlpha <- function(x) {
            temp <- x
            tempMax <- base::max(temp)
            index <- base::which(temp == tempMax)
            temp <- temp[-index]
            deltaMax <- base::max(temp) / tempMax
            return(1.0 - base::abs(deltaMax))
        }

        # Returns a color for each cell in a grey to 'cell type base color' color scheme, indicating the relative confidence of the annotation for a particular cell among all other cells of the same cell type
        cTIdfConfCol <- function(x, index, bC) {
            colorVec <- grDevices::colorRampPalette(base::c("grey", bC))(50)
            maxVal <- base::max(x[, index])
            maxIndex <- base::which(x[, index] == maxVal)
            cellTypeOrder <- base::order(x[maxIndex, ])
            ratio <-
                base::max(1, (base::which(cellTypeOrder == index) / base::length(cellTypeOrder)) * 50)
            result <- colorVec[ratio]
            return(result)
        }

        # Compute cell type assignments and confidence Scores (alpha values for transparency and relative color scale).
        cellTypes <- base::list()
        confidenceScore <- base::list()
        relativeColorRank <- base::list()
        for (i in base::c(1:base::dim(rca.obj$projection.data)[2])) {
            if (base::is.null(confidence)) {
                cellTypes <- base::c(cellTypes, cTIdfWU(rca.obj$projection.data[, i]))
            }
            else{
                cellTypes <-
                    base::c(cellTypes,
                      cTIdf(rca.obj$projection.data[, i], confidence))
            }
        }
        if (cSCompute) {
            for (i in base::c(1:base::dim(rca.obj$projection.data)[2])) {
                confidenceScore <-
                    base::c(confidenceScore,
                      cTIdfAlpha(rca.obj$projection.data[, i]))
            }
            rca.obj$cScore <- confidenceScore
        } else{
            rca.obj$cScore <- base::list()
        }
        if (ctRank) {
            myColors <-
                randomcoloR::distinctColorPalette(base::length(base::unique(cellTypes)))
            base::names(myColors) <- base::unique(cellTypes)
            baseColors <- myColors[base::unlist(cellTypes)]
            rca.obj$baseColors <- base::list(Colors = baseColors)
            for (i in c(1:dim(rca.obj$projection.data)[2])) {
                relativeColorRank <-
                    base::c(
                        relativeColorRank,
                        cTIdfConfCol(rca.obj$projection.data, i, baseColors[i])
                    )
            }
            rca.obj$rRank <- relativeColorRank
        } else{
            rca.obj$rRank <- base::list()
            rca.obj$baseColors <- base::list()
        }


        # Assign projection result to RCA object
        rca.obj$cell.Type.Estimate.per.cell <- cellTypes


        # Return RCA object
        return(rca.obj)
    }
