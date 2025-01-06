
#' Find deletion strains with similar expression profiles to a query strain using hypergeometric enrichment tests on signature genes
#'
#' @param strain Name of Deleteome strain to analyze (query strain/query gene)
#' @param outputdir Directory in which to save similarity analysis results
#' @param minAbsLog2FC Log2 fold-change cutoff used to classify genes as differentially-expressed (the absolute value of log2 fold-change must be higher than minAbsLog2FC)
#' @param pDEGs P-value cutoff used to identify differentially-expressed genes
#' @param pEnrich P-value cutoff used to identify statistically significant enrichment tests
#' @param quantileCutoff Quantile cutoff for selecting the Deleteome strain matches with highest confidence
#' @param returnTestValues  If true, a data frame containing the correlation test results against each similar Deleteome strain is returned. If false, only the names of the similar strains are returned.
#' @param showMessages Whether to output progress messages
#'
#' @export
#'
#' @examples
#' # getSimilarStrainsByEnrichment(strain = "nup170",
#' #                               outputdir = "...output directory path...",
#' #                               minAbsLog2FC = 0,
#' #                               pDEGs = 0.05,
#' #                               pEnrich = 0.05,
#' #                               quantileCutoff = 0.1)

getSimilarStrainsByEnrichment <- function(strain=NA,
                                          outputdir = NA,
                                          minAbsLog2FC=0,
                                          pDEGs=0.05,
                                          pEnrich=0.05,
                                          quantileCutoff=0.05,
                                          returnTestValues=F,
                                          showMessages = F){

  message("Getting deleteome matches based on enrichment tests...")

  conds <- getAllStrainNames()

  if( ! strain %in% conds){
    stop(paste0(strain, " is not a strain in the Deleteome"))
  }

  outputdir <- file.path(outputdir)

  if(is.na(outputdir)){
    message("ERROR: Please specify the directory in which to save results of similarity analysis")
    return(invisible(NULL))
  }
  else if( ! dir.exists(outputdir)){
    message("ERROR: The specified output directory ", outputdir, " does not exist.")
    return(invisible(NULL))
  }


  hypergs <- data.frame(Condition=as.character(),HyperGpval=as.numeric(), Pvalue.FDR=as.numeric(),
                        sampleSize=as.numeric(), sampleSuccesses=as.numeric(), popSize=as.numeric(),
                        popSuccesses=as.numeric(), stringsAsFactors = F)

  conds <- getAllStrainNames()

  nconds <- length(conds)

  ndeleteomegenes <- length(getStrainSignature(strain, 0, 1, consoleMessages = showMessages)$systematicName)

  # Get signature for query strain
  strainSignature <- getStrainSignature(strain, minAbsLog2FC, pDEGs, consoleMessages = showMessages)
  strainSignatureALL <- getStrainSignature(strain, 0, 1.0, consoleMessages = showMessages)

  if(dim(strainSignature)[1]==0){
    message(c("\nERROR: ", strain," had no significantly changed genes, based on M-value and p-value thresholds"))
    return(c())
  }

  # compute hypergeometric test to find strains whose signatures are enriched for the strain of interest's signature
  for(cond in conds){

    if(cond==strain | cond %in% c("wt_matA","wt_by4743","wt_ypd")){
      next() # Skip the query strain and WT control experiments in the deleteome
    }

    signature <- getStrainSignature(cond, minAbsLog2FC, pDEGs, consoleMessages = showMessages)
    signatureint <- computeDirectionalMatches(signature, strainSignature)

    samplesuccesses <- length(signatureint[signatureint$samedir>0,"systematicName"])
    samplesize <- length(signature$systematicName)
    popsuccesses <- length(strainSignature$systematicName)
    popsize <- ndeleteomegenes

    hyperp <- phyper(samplesuccesses-1, popsuccesses, (popsize-popsuccesses), samplesize, lower.tail=F) # the minus 1 is because probabilities are P[X>x] by default but we want P[X>=x]
    hyperpFDR <- 1 # set later

    if( ! is.na(hyperp)){
      hypergs <- rbind(hypergs,data.frame(Condition=cond, HyperGpval=hyperp, Pvalue.FDR=hyperpFDR,
                                          sampleSize=samplesize, sampleSuccesses=samplesuccesses,
                                          popSize=popsize, popSuccesses=popsuccesses, stringsAsFactors = F))
    }
  }

  hypergs$Pvalue.FDR <- p.adjust(hypergs$HyperGpval, method = "BH")

  hypergs <- addQuantileLevel(hypergs)

  pctcutoff <- quantile(hypergs$Pvalue.FDR, quantileCutoff, type = 1) # make sure that we are only using the top X% of p-values
  sighypergs <- hypergs[hypergs$Pvalue.FDR <= pEnrich & hypergs$Pvalue.FDR <= pctcutoff, ] # X% cutoff and must meet significance criteria

  similarStrains <- sighypergs[order(sighypergs$Pvalue.FDR, decreasing=F),"Condition"]

  if(length(similarStrains)==0){
    message("Could not find any deletion strains with signatures that significantly overlapped with ", strain, " deletion")
    return(c())
  }

  sigresults <- sighypergs[order(sighypergs$Pvalue.FDR, decreasing=F), c("Condition", "HyperGpval", "Pvalue.FDR", "Pvalue.FDR.quantile")]
  sigresultsfile <- paste0(outputdir, "/", strain, "_sigHyperG_L2FC",minAbsLog2FC,"_pDEGs",pDEGs,"_pEnrich",pEnrich,"_quantile",quantileCutoff,".tsv")
  write.table(sigresults, file = sigresultsfile, sep="\t", quote=F, row.names = F, col.names = T)
  message("Results written to ", sigresultsfile)

  resultsToOutput <- sigresults

  if(returnTestValues) return(resultsToOutput)
  else return(similarStrains) # Returns list of Deleteome strains matching the input strain
}
