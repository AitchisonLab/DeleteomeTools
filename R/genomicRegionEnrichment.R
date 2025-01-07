
#' Perform genomic region enrichment tests
#'
#' @param genePositions Table of genes and genomic positions (obtained using getGenePositions() if not specified)
#' @param systematicNamesBG Background set of systematic gene names to use for enrichment tests
#' @param systematicNames Set of systematic gene names to test for regional enrichment
#' @param includeMito Whether to include mitochondrial genes in the analysis
#' @param relativeTo Genomic region that will be tested for enrichment of input genes. Valid values are "telomere" or "centromere".
#' @param rangeInKB Window (in kilobases) for defining a gene as located in the region
#'
genomicRegionEnrichment <- function(genePositions=NULL,
                                    systematicNamesBG,
                                    systematicNames,
                                    includeMito=F,
                                    relativeTo="telomere",
                                    rangeInKB=25){

  if( ! checkLogical("includeMito", includeMito)) return(invisible(NULL))

  if( ! is.data.frame(genePositions) | is.null(dim(genePositions)[1])){
    message("Using default gene position data included in package.")
    genePositions <- getGenePositions(includeMito=includeMito)
  }

  if( ! is.character(systematicNamesBG)){
    message("\nERROR: Genomic region enrichment test not performed because input background gene set was not a character vector.")
    return(invisible(NULL))
  }

  if( ! is.character(systematicNames)){
    message("\nERROR: Genomic region enrichment test not performed because input gene set to test was not a character vector.")
    return(invisible(NULL))
  }

  dfname = "dist_from_telo"
  dfkbname = "dist_from_telo_in_kb"

  if(relativeTo=="centromere"){
    dfname = "dist_from_cent"
    dfkbname = "dist_from_cent_in_kb"
  }else if(relativeTo != "telomere"){
    message("\nERROR: Invalid value for relativeTo parameter. Please use either \"telomere\" or \"centromere\".")
    return(invisible(NULL))
  }

  if( ! checkNumeric("rangeInKB", rangeInKB, 0)) return(invisible(NULL))

  # Limit the background set of genes to only those in systematicNamesBG
  genePositions <- genePositions[which(genePositions$Geneid %in% systematicNamesBG),]

  res <- data.frame(systematicName=as.character(systematicNames), x=as.numeric(NA),stringsAsFactors=F)
  names(res)[2] <- dfkbname

  # Record distance from region for each gene in results
  for( geneid in systematicNames ){

    rowWithDist = which(genePositions[,1]==geneid)

    if(length(rowWithDist) == 1 ){
      df = genePositions[rowWithDist,dfname]
      res[res$systematicName==geneid, dfkbname] <- df/1000	 # convert from bases to kilobases
    }
    else{
      message(paste0("WARNING: number of matching rows in gene positions table for gene ", geneid, " did not equal 1"))
    }
  }

  message("Done collecting gene distances from region")
  res <- na.omit(res)

  # compute hypergeometrics to see if list of signficantly upregulated genes are enriched for genes near the region of interest
  samplesuccesses <- length(which(res[,dfkbname]<=rangeInKB))
  samplesize <- length(res[,dfkbname])

  popsuccesses <- length(which(genePositions[,dfname]<=rangeInKB*1000))
  popsize <- length(genePositions[,dfname])

  hyperp = phyper(samplesuccesses-1, popsuccesses, (popsize-popsuccesses), samplesize, lower.tail=F) # the minus 1 is because probabilities are P[X>x] by default but we want P[X>=x]

  message(c("\nPopulation size: ", popsize, " Population successes: ", popsuccesses, " Sample size: ", samplesize, " Sample successes: ", samplesuccesses))

  return(hyperp)
}



#' Convenience function to test for enrichment of differentially-expressed genes in the centromeric region.
#' @param genePositions Table of genes and genomic positions (obtained using getGenePositions() if not specified)
#' @param systematicNamesBG Background set of systematic gene names to use for enrichment tests
#' @param systematicNames Set of systematic gene names to test for centromeric enrichment
#' @param includeMito Whether to include mitochondrial genes in the analysis
#' @param rangeInKB Window (in kilobases) for defining a gene as located in the centromeric region
#'
#' @export
#'
centromericEnrichment <- function(genePositions=NULL,
                                  systematicNamesBG,
                                  systematicNames,
                                  includeMito=F,
                                  rangeInKB=25){

  return(genomicRegionEnrichment(genePositions=genePositions,
                                 systematicNamesBG=systematicNamesBG,
                                 systematicNames=systematicNames,
                                 includeMito=includeMito,
                                 relativeTo="centromere",
                                 rangeInKB=rangeInKB))
}



#' Convenience function to test for enrichment of differentially-expressed genes in the subtelomeric region
#' @param genePositions Table of genes and genomic positions (obtained using getGenePositions() if not specified)
#' @param systematicNamesBG Background set of systematic gene names to use for enrichment tests
#' @param systematicNames Set of systematic gene names to test for subtelomeric enrichment
#' @param includeMito Whether to include mitochondrial genes in the analysis
#' @param rangeInKB Window (in kilobases) for defining a gene as located in the subtelomeric region
#'
#' @export
#'
subtelomericEnrichment <- function(genePositions=NULL,
                                   systematicNamesBG,
                                   systematicNames,
                                   includeMito=F,
                                   rangeInKB=25){

  return(genomicRegionEnrichment(genePositions=genePositions,
                                 systematicNamesBG=systematicNamesBG,
                                 systematicNames=systematicNames,
                                 includeMito=includeMito,
                                 relativeTo="telomere",
                                 rangeInKB=rangeInKB))
}

