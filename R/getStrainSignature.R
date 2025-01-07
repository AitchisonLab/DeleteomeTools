#' Get data for a deletion strain's differentially-expressed genes
#'
#' @param strain Name of deletion strain as one-element character vector
#' @param minAbsLog2FC The absolute log2 fold-change threshold when selecting differentially-expressed genes (use 0 to include all genes)
#' @param pDEGs The p-value threshold to use when selecting differentially-expressed genes (use 1 to include all genes)
#' @param consoleMessages Whether to output console messages
#'
#' @return A data frame with systematic gene ids, gene names, M values (log2 fold-changes vs. WT), and p-values for a given deletion strain's transcriptome. Only microarrayed genes that meet M value and p-value cutoffs are included.
#' @export
#'
#' @examples
#' # getStrainSignature(strain = "nup170", minAbsLog2FC = 0, pDEGs = 0.05)
#'
getStrainSignature <- function(strain = "",
                               minAbsLog2FC = 0,
                               pDEGs = 0.05,
                               consoleMessages=T){

  if( ! is.character(strain)){
    message("ERROR: Please enter a character string for the strain parameter")
    return(invisible(NULL))
  }

  if( ! checkNumeric("minAbsLog2FC", minAbsLog2FC, minValue = 0) | ! checkNumeric("pDEGs", pDEGs, minValue = 0, maxValue = 1.0)) return(invisible(NULL))
  if( ! checkLogical("consoleMessages", consoleMessages)) return(invisible(NULL))

  if(consoleMessages){
    message(paste0("Getting log2 fold-change values for ", strain, " from deleteome..."))
  }

  McolSuffix <- "_wt"
  Mcolname <- paste0(strain,"_del_vs",McolSuffix)
  pcolname <- paste0(strain,"_del_vs_wt_2")

  if(Mcolname %in% names(deldata) & pcolname %in% names(deldata)){
    selected <- deldata[abs(as.numeric(deldata[,Mcolname]))>=minAbsLog2FC & as.numeric(deldata[,pcolname])<=pDEGs,c("systematicName","geneSymbol",Mcolname,pcolname)]
    selected <- selected[selected$geneSymbol!=toupper(strain),]
    selected[,3] <- as.numeric(selected[,3])
    selected[,4] <- as.numeric(selected[,4])
  }
  else{
    if(consoleMessages){
      message(c("Could not find expression and/or p-values for ", Mcolname, " and ", pcolname))
    }
    temp <- data.frame(matrix(nrow=0,ncol=4))
    names(temp) <- c("systematicName","geneSymbol",Mcolname,pcolname)
    return(temp)
  }
  return(selected)

}
