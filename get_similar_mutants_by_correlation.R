scriptloc <- dirname(rstudioapi::getActiveDocumentContext()$path) # Get script location
setwd(scriptloc)

# load deleteome_tools.R script
source(paste0(scriptloc, "/deleteome_tools.R")) # Assumes that deleteome_tools.R is in same directory as this file

alldata <- getDeleteomeExpData() # Load all the deleteome expression data

Mthresh <- 0 # log2 fold-change cutoff
pthresh <- 0.05 # p-value cutoff

mutantname <- "nup170"

selectedConditions <- getDeleteomeMatchesByReciprocalCorrelation(mutant=mutantname, 
                                                                 minAbsLog2FC = Mthresh, 
                                                                 pCutoff = pthresh, 
                                                                 deleteomeData = alldata)

# Show correlated mutants
print(selectedConditions)

# Generate heatmap for significantly correlated deleteome mutants
# First get a data frame that contains all the log2 fold-change and p-value info for the mutant strain 
# (AKA its "profile")
mutantProfile <- getProfileForDeletion(delData = alldata, 
                                       deletionname = "nup170", 
                                       Mthresh = Mthresh, 
                                       pthresh = pthresh)

hm1 <- makeHeatmapDeleteomeMatches(mutantname=mutantname, 
                                   mutantProfile, 
                                   selectedConditions, 
                                   fileprefix = "sigCorrMutants", 
                                   titledesc="reciprocal correlation", 
                                   MthreshForTitle = Mthresh, 
                                   pthreshForTitle = pthresh, 
                                   printToFile = T)

# Generate heatmap using subtelomic genes only
hm2 <- makeHeatmapDeleteomeMatches(mutantname=mutantname, 
                                   mutantProfile, 
                                   selectedConditions, 
                                   fileprefix = "sigCorrMutantsSubtelo", 
                                   titledesc="reciprocal correlation (only subtelo genes shown)",
                                   MthreshForTitle = Mthresh, 
                                   pthreshForTitle = pthresh, 
                                   subteloGenesOnly = T, 
                                   rowFontSize=0.275, 
                                   printToFile = T)

# Generate heatmap using specific deleteome profiles
hm3 <- makeHeatmapDeleteomeMatches(mutantname=mutantname, mutantProfile, 
                                   c("hmo1","rif1","sir4","ctf8","ctf18","dcc1"), 
                                   fileprefix = "specificMutants", 
                                   titledesc="specific deleteome profiles", 
                                   MthreshForTitle = Mthresh, pthreshForTitle = pthresh,
                                   printToFile = T)

# Make a mountain lake plot (this also reports subtelomeric enrichment p-values for the mutant's up- 
# and down-regulated genes)
makeGenomicPositionHistogram(alldata = alldata, 
                             mutant = mutantname, 
                             Mthresh = Mthresh, 
                             xmax=770, ymax = 40, 
                             upcolor="#d53e4f", 
                             downcolor="#3288bd", 
                             showupdownlabels = T)


# Run GO enrichment on similar mutants
GOpadjcutoff = 0.05
GO <- doGOenrichmentOnDeleteomeMatches(alldata, 
                                       selectedConditions, 
                                       pthresh = GOpadjcutoff)
GOoutputfile <- paste0(scriptloc,"/output/GO_enrichment/", mutantname, "GO_enrichmentResults_RecipCorr_GOpadj",GOpadjcutoff,".tsv")
write.table(GO, file=GOoutputfile, sep="\t", quote = F, row.names = F)
message("Wrote GO enrichment results to ", GOoutputfile)
