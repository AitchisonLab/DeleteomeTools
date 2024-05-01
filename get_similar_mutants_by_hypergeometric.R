scriptloc <- dirname(rstudioapi::getActiveDocumentContext()$path) # Get script location
setwd(scriptloc)

# load deleteome_tools.R script
source(paste0(scriptloc, "/deleteome_tools.R")) # Assumes that deleteome_tools.R is in same directory as this file

alldata <- getDeleteomeExpData() # Load all the deleteome expression data
# getAllStrainNames(alldata) # Use this function to output a list of all available mutant strains in the Deleteome

mutantname <- "nup170" # The tool will find mutant strains transcriptionally similar to the strain specified here

Mthresh <- 0.0 # log2 fold-change cutoff for identifying differentially-expressed genes in deletion strain
pDEGs <- 0.05 # p-value cutoff for identifying differentially-expressed genes in deletion strain
pEnrich <- 0.05  # p-value cutoff for identifying statistically-significant correlation tests
qthresh <- 0.05 # quantile cutoff for selecting similar deletion strains

# Find similar strains in Deleteome
matchingStrains <- getDeleteomeMatchesByEnrichment(   delData = alldata,
                                                      mutant=mutantname, 
                                                      minAbsLog2FC = Mthresh, 
                                                      pDEGs = pDEGs,
                                                      pEnrich = pEnrich,
                                                      quantileCutoff = qthresh
                                                      )

# Show correlated mutants
message("\nSimilar mutant strains:")
print(matchingStrains)

# Generate heatmap for significantly correlated deleteome mutants
# First get a data frame that contains all the log2 fold-change and p-value info for the mutant strain (AKA its "profile")
mutantProfile <- getProfileForDeletion(delData = alldata, 
                                       deletionname = mutantname, 
                                       Mthresh = Mthresh, 
                                       pDEGs = pDEGs)

# Generate heatmap for significantly similar deleteome profiles
hm1 <- makeHeatmapDeleteomeMatches(mutantname = mutantname, 
                                   mutantProfile, 
                                   matchingStrains, 
                                   fileprefix = "HyperG_matches", 
                                   titledesc="transcriptional similarity (signature enrichment)", 
                                   MthreshForTitle = Mthresh, 
                                   pDEGsForTitle = pDEGs,
                                   pMatchesForTitle = pEnrich,
                                   quantileForTitle = qthresh,
                                   imagewidth = 5000,
                                   printToFile = T)

# Generate heatmap for significantly correlated deleteome strains that only shows subtelomeric genes
# NOTE: if there are less than 2 subtelomeric genes in the strain's signature, 
# this will generate an error indicating insufficient rows or columns for the heatmap
hm2 <- makeHeatmapDeleteomeMatches(mutantname = mutantname, 
                                   mutantProfile, 
                                   matchingStrains, 
                                   fileprefix = "HyperG_matches_subtelo", 
                                   titledesc="transcriptional similarity (signature enrichment)", 
                                   MthreshForTitle = Mthresh, 
                                   pDEGsForTitle = pDEGs,
                                   pMatchesForTitle = pEnrich,
                                   quantileForTitle = qthresh,
                                   subteloGenesOnly = T, 
                                   imagewidth = 5000,
                                   printToFile = T)

# Generate heatmap using specific deleteome profiles
hm3 <- makeHeatmapDeleteomeMatches(mutantname = mutantname, 
                                   mutantProfile, 
                                   c("hmo1","rif1","sir4","ctf8","ctf18","dcc1"), 
                                   fileprefix = "Specific_mutants", 
                                   titledesc="manual selection", 
                                   MthreshForTitle = Mthresh, 
                                   pDEGsForTitle = pDEGs,
                                   pMatchesForTitle = pEnrich,
                                   quantileForTitle = qthresh,
                                   imagewidth = 2550,
                                   printToFile = T)


# Make a mountain lake plot (this also reports subtelomeric enrichment p-values for the mutant's up- and down-regulated genes)
makeGenomicPositionHistogram(delData = alldata, 
                             mutant = mutantname, 
                             Mthresh = Mthresh,
                             pDEGs = pDEGs,
                             xmax=770, ymax = 40, 
                             upcolor="#d53e4f", 
                             downcolor="#3288bd")


# Run GO enrichment on similar mutants
GOpadjcutoff = 0.05
GO <- doGOenrichmentOnDeleteomeMatches(delData = alldata, 
                                       genes = matchingStrains, 
                                       padjthresh = GOpadjcutoff)
GOoutputfile <- paste0(scriptloc,"/output/GO_enrichment/", mutantname, "_GOenrichmentResults_HyperG_GOpadj", GOpadjcutoff,".tsv")
write.table(GO, file=GOoutputfile, sep="\t", quote = F, row.names = F)
message("Wrote GO enrichment results to ", GOoutputfile)
