#' Create RCA Object from 10X data
#'
#' @param dataDir Directory containing 10X data
#' @param min.barcode.umi Minimum UMIs needed for a barcode to be considered. Default is 100.
#' @return RCA object.
#' @export
#'
createRCAObjectFrom10X <- function(dataDir, min.barcode.umi = 100) {

    # Check if data directory exists, otherwise stop and throw error
    if (!base::dir.exists(paths = dataDir)) {
        stop("10X Directory provided does not exist")
    }

    # Set barcodes file path
    if(base::file.exists(base::file.path(dataDir, "barcodes.tsv"))) {
        barcodeFilePath <- base::file.path(dataDir, "barcodes.tsv")
    } else if (base::file.exists(base::file.path(dataDir, "barcodes.tsv.gz"))) {
        barcodeFilePath <- base::file.path(dataDir, "barcodes.tsv.gz")
    } else {
        stop("Barcodes file not found.")
    }

    # Set genes/features file path
    if(base::file.exists(base::file.path(dataDir, "genes.tsv"))) {
        featureFilePath <- base::file.path(dataDir, "genes.tsv")
    } else if (base::file.exists(base::file.path(dataDir, "genes.tsv.gz"))) {
        featureFilePath <- base::file.path(dataDir, "genes.tsv.gz")
    } else if (base::file.exists(base::file.path(dataDir, "features.tsv"))) {
        featureFilePath <- base::file.path(dataDir, "features.tsv")
    } else if (base::file.exists(base::file.path(dataDir, "features.tsv.gz"))) {
        featureFilePath <- base::file.path(dataDir, "features.tsv.gz")
    } else {
        stop("Genes/Features file not found.")
    }

    # Set matrix file path
    if(base::file.exists(base::file.path(dataDir, "matrix.mtx"))) {
        matrixFilePath <- base::file.path(dataDir, "matrix.mtx")
    } else if (base::file.exists(base::file.path(dataDir, "matrix.mtx.gz"))) {
        matrixFilePath <- base::file.path(dataDir, "matrix.mtx.gz")
    } else {
        stop("Matrix file not found.")
    }

    # Load 10X barcodes and feature names
    cellNames <- base::readLines(barcodeFilePath)
    featureNames <- utils::read.table(featureFilePath, header = FALSE, stringsAsFactors = FALSE)[[2]]

    # Load 10X data
    rawData <- Matrix::readMM(file = matrixFilePath)
    base::rownames(rawData) <- featureNames
    base::colnames(rawData) <- cellNames

    # Filter barcodes by minimum barcode UMI
    nUMIVec <- Matrix::colSums(rawData)
    filt.cells <- base::which(nUMIVec >= min.barcode.umi)
    rawData <- rawData[, filt.cells]

    # Create RCA object using RCAConstruct and the raw data provided
    rca.obj <- methods::new(Class="RCA",raw.data = rawData)

    # Return RCA object
    return(rca.obj)
}
