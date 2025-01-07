
#' Find deletion strains with similar expression profiles to a query strain using reciprocal correlation
#'
#' @param strain Name of Deleteome strain to analyze (query strain/query gene)
#' @param outputDir Directory in which to save similarity analysis results
#' @param minAbsLog2FC Log2 fold-change cutoff used to classify genes as differentially-expressed (the absolute value of log2 fold-change must be higher than minAbsLog2FC)
#' @param pDEGs P-value cutoff used to identify differentially-expressed genes
#' @param pCor P-value cutoff used to identify statistically significant correlation tests
#' @param quantileCutoff Quantile cutoff for selecting the Deleteome strain matches with highest confidence
#' @param returnTestValues  If true, a data frame containing the correlation test results against each similar Deleteome strain is returned. If false, only the names of the similar strains are returned.
#' @param showMessages Whether to output progress messages
#'
#' @export
#'
#' @examples
#' # getSimilarStrainsByReciprocalCorrelation(strain = "nup170",
#' #                                          outputDir = "...output directory path...",
#' #                                          minAbsLog2FC = 0,
#' #                                          pDEGs = 0.05,
#' #                                          pCor = 0.05,
#' #                                          quantileCutoff = 0.1)
#'
getSimilarStrainsByReciprocalCorrelation <- function( strain="",
                                                      outputDir = NA,
                                                      minAbsLog2FC=0,
                                                      pDEGs=0.05,
                                                      pCor=0.05,
                                                      quantileCutoff=0.1,
                                                      returnTestValues=F,
                                                      showMessages = F){

  message("\n\nFinding Deleteome strains transcriptionally similar to ", strain, " deletion strain by reciprocal correlation...")

  conditions <- getAllStrainNames()

  if( is.character(strain)){
    if( ! strain %in% conditions){
      message(paste0("Could not identify similar strains: ", strain, " is not a strain in the Deleteome"))
      return(invisible(NULL))
    }
  }
  else{
    message(paste0("ERROR: Please enter a character string for the strain parameter. Use getAllStrainNames() to see a list of valid strain names."))
    return(invisible(NULL))
  }

  outputDir <- file.path(outputDir)

  if( ! is.character(outputDir)){
    message("ERROR: Please enter a character string for the directory in which to save results of similarity analysis.")
    return(invisible(NULL))
  }
  else if( ! dir.exists(outputDir)){
    message("ERROR: The specified output directory ", outputDir, " does not exist.")
    return(invisible(NULL))
  }

  # Check numerical parameters
  if( ! all( c(checkNumeric("minAbsLog2FC", minAbsLog2FC, minValue = 0.0),
               checkNumeric("pDEGs", pDEGs, minValue = 0.0, maxValue = 1.0),
               checkNumeric("pCor", pCor, minValue = 0.0, maxValue = 1.0),
               checkNumeric("quantileCutoff", quantileCutoff, minValue = 0.0, maxValue = 1.0)))){

    return(invisible(NULL))
  }

  # Check logical parameters
  if( ! checkLogical("returnTestValues", returnTestValues)) return(invisible(NULL))
  if( ! checkLogical("showMessages", showMessages)) return(invisible(NULL))



  # Get the query strain's profile (signature)
  strainSignature <- getStrainSignature(strain, minAbsLog2FC, pDEGs, consoleMessages = showMessages)

  if(dim(strainSignature)[1] <= 2){
    message(c("\nERROR: ", strain," had ", dim(strainSignature)[1],
              " differentially-expressed genes based on entered log2 fold-change and p-value thresholds, but a minimum of 3 is required for reciprocal correlation analysis."))
    return(invisible(NULL))
  }

  # get the direct correlations between the strain and all other deleteome strains
  allCorrResults <- data.frame(Deletion=as.character(),
                               CorrCoefficient=as.numeric(),
                               Pvalue=as.numeric(),
                               stringsAsFactors=F)
  h = 1

  for(cond in conditions){

    if(cond %in% c("wt_matA","wt_by4743","wt_ypd")) next() # Skip the WT control experiments in the deleteome

    allconddata <- getStrainSignature(cond, 0, 1, consoleMessages = showMessages)
    intrsct <- intersect(strainSignature$systematicName, allconddata$systematicName) # need to do this b/c condition profile won't include the KO'd gene, which might be in the strain's profile

    if(length(intrsct) > 2){
      # Perform correlation test
      correl <- cor.test(allconddata[allconddata$systematicName %in% intrsct, 3], strainSignature[strainSignature$systematicName %in% intrsct, 3])
      Rval <- as.vector(correl["estimate"][[1]])
      pval <- as.vector(correl["p.value"][[1]])
    }
    else{
      # Not enough overlap in signatures for the strain
      Rval <- NA
      pval <- NA
    }

    allCorrResults[h,] <- c(cond, Rval, pval)
    h = h+1
  }

  allCorrResults$Pvalue.FDR <- as.numeric(p.adjust(allCorrResults$Pvalue, method = "BH"))  # FDR-adjustment

  allCorrResultsNoNA <- allCorrResults[ ! is.na(allCorrResults$CorrCoefficient) & ! is.na(allCorrResults$Pvalue), ]
  allCorrResultsNoNA$CorrCoefficient <- as.numeric(allCorrResultsNoNA$CorrCoefficient)
  allCorrResultsNoNA$Pvalue <- as.numeric(allCorrResultsNoNA$Pvalue)
  allCorrResultsNoNA$Pvalue.FDR <- as.numeric(allCorrResultsNoNA$Pvalue.FDR)

  allSigCorrResults <- allCorrResultsNoNA[allCorrResultsNoNA$CorrCoefficient > 0, ] # limit to positive correlations

  # Add quantile level for FDR-adjusted p-values
  allSigCorrResults <- addQuantileLevel(allSigCorrResults)

  pvalcutoff <- quantile(allSigCorrResults$Pvalue.FDR, quantileCutoff, type = 1) # get strains with p-values that were in the desired quantile
  # message("Quantile-based FDR P-value cutoff set to ", pvalcutoff)

  # Sometimes the pvalcutoff can be zero, which will incorrectly return no significant strains. Catch that case here.
  if(pvalcutoff == 0) allSigCorrResults <- allSigCorrResults[allSigCorrResults$Pvalue.FDR == pvalcutoff & allSigCorrResults$Pvalue.FDR < pCor,]
  else allSigCorrResults <- allSigCorrResults[allSigCorrResults$Pvalue.FDR < pvalcutoff & allSigCorrResults$Pvalue.FDR < pCor,]

  allSigCorrResults <- allSigCorrResults[order(allSigCorrResults$CorrCoefficient, decreasing=T),]

  # for each signifcantly correlated deletion strain, see if reciprocal correlation is also significant
  allRecipCorrResults = data.frame(matrix(ncol=3, nrow=length(allSigCorrResults$Deletion)))
  names(allRecipCorrResults) = c("Deletion", "CorrCoefficient","Pvalue")

  m = 1

  strainSignatureANY <- getStrainSignature(strain, 0, 1, consoleMessages = showMessages)

  for(cond in allSigCorrResults$Deletion){

    condprofile <- getStrainSignature(cond, minAbsLog2FC, pDEGs, consoleMessages = showMessages)
    intrsct <- intersect(strainSignatureANY$systematicName,condprofile$systematicName) # need to do this b/c condition profile won't include the KO'd gene, which might be in the strain's profile

    if(length(intrsct) > 2){
      correl <- cor.test(condprofile[condprofile$systematicName %in% intrsct,3],strainSignatureANY[strainSignatureANY$systematicName %in% intrsct,3])
      Rval <- as.vector(correl["estimate"][[1]])
      pval <- as.vector(correl["p.value"][[1]])
    }
    else{
      Rval <- NA
      pval <- NA
    }
    allRecipCorrResults[m,] <- c(cond, Rval, pval)
    m = m + 1
  }

  allRecipCorrResults$Pvalue.FDR <- p.adjust(allRecipCorrResults$Pvalue, method="BH")
  allRecipCorrResultsNoNA <- na.omit(allRecipCorrResults)

  # Find overlap between direct and reciprocal correlations. Limit results to significant, positive reciprocal correlations.
  allSigRecipCorrResults <- allRecipCorrResultsNoNA[as.numeric(allRecipCorrResultsNoNA$Pvalue.FDR) < pCor &
                                                      allRecipCorrResultsNoNA$CorrCoefficient > 0,]

  sigDirectRecipCorrNames <- intersect(allSigRecipCorrResults$Deletion, allSigCorrResults$Deletion) # get all KO's that had significant direct and reciprocal correlations
  sigDirectRecipCorrResults <- allSigCorrResults[allSigCorrResults$Deletion %in% sigDirectRecipCorrNames,]
  sigDirectRecipCorrResults <- sigDirectRecipCorrResults[order(sigDirectRecipCorrResults$CorrCoefficient, decreasing=T),] # sort by direct correlation Rval


  outputfilename <- paste0(outputDir, "/", strain, "_sigCorr_L2FC", minAbsLog2FC, "_pDEGs", pDEGs, "_pCor", pCor, "_quantile",quantileCutoff,".tsv")

  resultsToOutput <- sigDirectRecipCorrResults[sigDirectRecipCorrResults$Deletion != strain,
                                               c("Deletion", "CorrCoefficient", "Pvalue", "Pvalue.FDR", "Pvalue.FDR.quantile")]
  similarStrains <- resultsToOutput$Deletion

  write.table(resultsToOutput, outputfilename, row.names = F, col.names = T, sep="\t", quote=F)
  message("Results written to ", outputfilename)

  if(length(similarStrains)==0){
    message("Could not find any deletion strains with signatures that significantly overlapped with ", strain, " deletion")
  }

  if(returnTestValues) return(resultsToOutput)
  else return(similarStrains) # Returns list of Deleteome strains matching the input strain
}

