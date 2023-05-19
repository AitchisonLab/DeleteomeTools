scriptloc <- dirname(rstudioapi::getActiveDocumentContext()$path) # Get script location
setwd(scriptloc)

# load deleteome_tools.R script
source(paste0(scriptloc, "/deleteome_tools.R")) # Assumes that deleteome_tools.R is in same directory as this file

alldata <- getDeleteomeExpData() # Load all the deleteome expression data

Mthresh <- 0 # log2 fold-change cutoff
pthresh <- 0.05 # p-value cutoff

mutantname <- "nup170"

selectedConditions <- getDeleteomeMatchesByEnrichment(
                        mutant=mutantname, 
                        minAbsLog2FC = Mthresh, 
                        pCutoff = pthresh, 
                        deleteomeData = alldata)

# Show correlated mutants
print(selectedConditions)

mutantProfile <- getProfileForDeletion(alldata, mutantname, Mthresh, pthresh)

# Generate heatmap for significantly correlated deleteome profiles
hm1 <- makeHeatmapDeleteomeMatches(mutantname=mutantname, 
                                   mutantProfile, 
                                   selectedConditions, 
                                   fileprefix = "sigHypergMutants", 
                                   titledesc="hypergeometric analysis", 
                                   MthreshForTitle = Mthresh, 
                                   pthreshForTitle = pthresh, 
                                   rowFontSize = 0.3, 
                                   printToFile = T)

# Generate heatmap for significantly correlated deleteome profiles (subtelomic genes only)
hm2 <- makeHeatmapDeleteomeMatches(mutantname=mutantname, 
                                   mutantProfile, 
                                   selectedConditions, 
                                   fileprefix = "sigCorrHypergSubtelo", 
                                   titledesc="hypergeometric analysis", 
                                   MthreshForTitle = Mthresh, 
                                   pthreshForTitle = pthresh, 
                                   subteloGenesOnly = T, 
                                   printToFile = T, 
                                   rowFontSize=0.275)

# Mountain lake plot
makeGenomicPositionHistogram(alldata = alldata, 
                             mutant = mutantname, 
                             Mthresh = Mthresh, 
                             upcolor="#d53e4f", 
                             downcolor="#3288bd")

# Do GO enrichment on deleteome matches
GOpadjcutoff = 0.05
GO <- doGOenrichmentOnDeleteomeMatches(alldata, 
                                       selectedConditions, 
                                       pthresh = GOpadjcutoff)
GOoutputfile <- paste0(scriptloc,"/output/GO_enrichment/",mutantname,"GOenrichmentResults_HyperGguilt_GOpadj",GOpadjcutoff,".tsv")
write.table(GO,file=GOoutputfile, 
            sep="\t", quote = F, row.names = F)
message("Wrote GO enrichment results to ", GOoutputfile)


