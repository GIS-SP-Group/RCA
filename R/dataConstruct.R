#' Construct data object
#' 
#' Construct a data object for the analysis
#' 
#' @param input_data A data frame of expression values (FPKM, TPM, UMI counts ...), with rows representing genes and columns representing cells. Note the current version of RCA only accepts gene names in the following format: "GenomeLocation_HGNCGeneName_EnsembleID", from which the "HGNCGeneName" is extracted for RCA analysis. For input data with only HGNC names, the users need to attach two strings to the HGNC names to make them into the "XXXX_HGNCGeneNames_YYYY" format. 
#' @return The data object for the following RCA analysis.
#' @export
#' @examples
#' input_data = read.csv("fpkm_data.csv",row.names=1);
#' data_obj = dataConstruct(input_data);
#' 
dataConstruct <- function(input_data)
{
  obj = list("fpkm_raw" = input_data
  );
  return(obj)
}