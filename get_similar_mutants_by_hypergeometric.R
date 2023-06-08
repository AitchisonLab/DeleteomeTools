scriptloc <- dirname(rstudioapi::getActiveDocumentContext()$path) # Get script location
setwd(scriptloc)

# load deleteome_tools.R script
source(paste0(scriptloc, "/deleteome_tools.R")) # Assumes that deleteome_tools.R is in same directory as this file

alldata <- getDeleteomeExpData() # Load all the deleteome expression data
# getAllStrainNames(alldata) # Use this function to output a list of all available mutant strains in the Deleteome

mutantname <- "nup170" # The tool will find mutant strains transcriptionally similar to the strain specified here

Mthresh <- 0 # log2 fold-change cutoff for identifying differentially-expressed genes in deletion strain
pthresh <- 0.05 # p-value cutoff for identifying differentially-expressed genes in deletion strain

selectedConditions <- getDeleteomeMatchesByEnrichment(
                        mutant=mutantname, 
                        minAbsLog2FC = Mthresh, 
                        pCutoff = pthresh, 
                        deleteomeData = alldata)

# Show correlated mutants
message("\nSimilar mutant strains:")
print(selectedConditions)

mutantProfile <- getProfileForDeletion(alldata, mutantname, Mthresh, pthresh)

# Generate heatmap for significantly correlated deleteome profiles
hm1 <- makeHeatmapDeleteomeMatches(mutantname=mutantname, 
                                   mutantProfile, 
                                   selectedConditions, 
                                   fileprefix = "sigHypergMutants", 
                                   titledesc="transcriptional similarity (hypergeometric analysis)", 
                                   MthreshForTitle = Mthresh, 
                                   pthreshForTitle = pthresh, 
                                   rowFontSize = 0.3, 
                                   printToFile = T)

# Generate heatmap for significantly correlated deleteome profiles (subtelomic genes only)
hm2 <- makeHeatmapDeleteomeMatches(mutantname=mutantname, 
                                   mutantProfile, 
                                   selectedConditions, 
                                   fileprefix = "sigHypergSubtelo", 
                                   titledesc="transcriptional similarity (hypergeometric analysis)", 
                                   MthreshForTitle = Mthresh, 
                                   pthreshForTitle = pthresh, 
                                   subteloGenesOnly = T, 
                                   printToFile = T, 
                                   rowFontSize=0.275)

# Generate heatmap using specific deleteome profiles
hm3 <- makeHeatmapDeleteomeMatches(mutantname=mutantname, mutantProfile, 
                                   c("hmo1","rif1","sir4","ctf8","ctf18","dcc1"), 
                                   fileprefix = "specificMutants", 
                                   titledesc="manual selection", 
                                   MthreshForTitle = Mthresh, 
                                   pthreshForTitle = pthresh, 
                                   imagewidth = 2000,
                                   printToFile = T)


# Make a mountain lake plot (this also reports subtelomeric enrichment p-values for the mutant's up- and down-regulated genes)
makeGenomicPositionHistogram(alldata = alldata, 
                             mutant = mutantname, 
                             Mthresh = Mthresh, 
                             xmax=770, ymax = 40, 
                             upcolor="#d53e4f", 
                             downcolor="#3288bd")

# Run GO enrichment on similar mutants
GOpadjcutoff = 0.05
GO <- doGOenrichmentOnDeleteomeMatches(alldata, 
                                       selectedConditions, 
                                       pthresh = GOpadjcutoff)
GOoutputfile <- paste0(scriptloc,"/output/GO_enrichment/",mutantname,"GOenrichmentResults_HyperGguilt_GOpadj",GOpadjcutoff,".tsv")
write.table(GO,file=GOoutputfile, sep="\t", quote = F, row.names = F)
message("Wrote GO enrichment results to ", GOoutputfile)
