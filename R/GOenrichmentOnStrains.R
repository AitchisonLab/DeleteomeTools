
#' Perform GO enrichment on a list of Deleteome strain names. By default, uses the list of strains in the Deleteome as the background for the enrichment tests. Performs enrichment test over GO:BP, GO:MF and GO:CC sub-ontologies.
#'
#' @param strains A character vector of Deleteome strain names
#' @param padjThresh FDR-adjusted P-value cutoff for GO enrichment significance
#' @param useDeleteomeBackground gene names, not systematic names
#'
#' @export
#'
#' @examples
#' # GOenrichmentOnStrains(strains = c("ctf18", "ctf8", "dcc1"), padjThresh = 0.1)
#'
GOenrichmentOnStrains <- function(strains=c(),
                                  padjThresh=0.1,
                                  useDeleteomeBackground = T){

  if(length(strains) == 0 | ! is.character(strains)){
    message("\nERROR performing GO enrichment: the input character vector of gene names is empty or invalid.")
    return(invisible(NULL))
  }

  if( ! checkNumeric("padjThresh", padjThresh, minValue = 0, maxValue = 1) |
      ! checkLogical("useDeleteomeBackground", useDeleteomeBackground)) return(invisible(NULL))

  message("Performing GO enrichment tests...")

  suppressMessages(suppressWarnings(require(clusterProfiler)))
  require(clusterProfiler)
  require(org.Sc.sgd.db)

  strainorfmap <- getStrainNameORFmap()

  # Note that the strain names WT_BY4743, WT_MATA, and WT_YPD will have ORF column values of <NA>
  # add will be omitted in the GO analysis

  strains <- toupper(strains)

  # Give warning if we couldn't find an associated ORF for a gene name
  for(i in 1:dim(strainorfmap)[1]){
    thestrain <- strainorfmap[i,"Strain"]
    if(is.na(strainorfmap[i,"ORF"]) & thestrain %in% strains){
      strains <- strains[strains != thestrain]  # Remove from analysis
      message("WARNING: Could not find ORF for strain ", strainorfmap[i,"Strain"], ". Removed from GO analysis.")
    }
  }

  orfs <- unique(strainorfmap[strainorfmap$Strain %in% strains, "ORF"])

  if(useDeleteomeBackground){
    allstrainorfs <- as.character(na.omit(unique(strainorfmap$ORF)))
    xGOUni <- clusterProfiler::enrichGO(orfs, pvalueCutoff = padjThresh, minGSSize = 1, OrgDb=org.Sc.sgd.db::org.Sc.sgd.db, keyType = "ORF", pAdjustMethod = "BH", universe = allstrainorfs, ont="ALL")
  }
  else{
    xGOUni <- clusterProfiler::enrichGO(orfs, pvalueCutoff = padjThresh, minGSSize = 1, OrgDb=org.Sc.sgd.db::org.Sc.sgd.db, keyType = "ORF", pAdjustMethod = "BH", ont="ALL")
  }
  #detach("package:clusterProfiler", unload=TRUE)
  return(xGOUni[ , ])
}
