#' Create RCA Object from 10X data
#'
#' @param dataDir Directory containing 10X data
#' @param cellrangerVersion Version of cellranger used to generate data. Default is 3.0.
#' @param min.barcode.umi Minimum UMIs needed for a barcode to be considered. Default is 100.
#' @return RCA object.
#' @export
#'
createRCAObjectFrom10X <- function(dataDir, cellrangerVersion = 3.0, min.barcode.umi = 100) {

    # Check if data directory exists, otherwise stop and throw error
    if (!dir.exists(paths = dataDir)) {
        stop("10X Directory provided does not exist")
    }

    # Matrix
    if (!require(Matrix))
        install.packages("Matrix", repos = "http://cran.us.r-project.org")
    require(Matrix)

    # Load file paths
    barcodeFilePath <- file.path(dataDir, "barcodes.tsv.gz")

    # For earlier versions of cellranger, the filenames would be different
    if(cellrangerVersion < 3.0) {
        barcodeFilePath <- file.path(dataDir, "barcodes.tsv")
        featureFilePath <- file.path(dataDir, "genes.tsv")
        matrixFilePath <- file.path(dataDir, "matrix.mtx")

    } else {
        barcodeFilePath <- file.path(dataDir, "barcodes.tsv.gz")
        featureFilePath <- file.path(dataDir, "features.tsv.gz")
        matrixFilePath <- file.path(dataDir, "matrix.mtx.gz")
    }

    # Stop if any of the files don't exist
    if(!file.exists(barcodeFilePath)) {
        stop("Barcode file missing")
    }
    if(!file.exists(featureFilePath)) {
        stop("Gene/Feature file missing")
    }
    if(!file.exists(matrixFilePath)) {
        stop("Matrix file missing")
    }

    # Load 10X barcodes and feature names
    cellNames <- readLines(barcodeFilePath)
    featureNames <- read.table(featureFilePath, header = FALSE, stringsAsFactors = FALSE)[[2]]

    # Load 10X data
    rawData <- readMM(file = matrixFilePath)
    rownames(rawData) <- featureNames
    colnames(rawData) <- cellNames

    # Filter barcodes by minimum barcode UMI
    nUMIVec <- Matrix::colSums(rawData)
    filt.cells <- which(nUMIVec >= min.barcode.umi)
    rawData <- rawData[, filt.cells]

    # Create RCA object using RCAConstruct and the raw data provided
    rca.obj <- RCAConstruct$new(raw.data = rawData, data = NULL)

    # Return RCA object
    return(rca.obj)
}
