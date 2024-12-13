# ------ CONFIGURATION ------ #
# If using RStudio, this will automatically source the deleteome_tools.R codebase, which is assumed to be in same folder as this script
if(Sys.getenv("RSTUDIO")=="1"){ dtdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
}else{
  dtdir <- "[...]"  # If not using RStudio, manually enter path to deleteome_tools.R parent folder (not the file) here
}
srcfile <- paste0(dtdir, "/deleteome_tools.R")

while( ! file.exists(srcfile)){ 
  message("\nERROR: Could not find Deleteome-Tools codebase at ", srcfile, "\nPlease enter full path to deleteome_tools.R file:")
  srcfile <- gsub("^\"|^\'|\"$|\'$", "", readline())  # Get path from user
  dtdir <- dirname(srcfile)
}

message("Location of deleteome_tools.R set to ", srcfile)
setwd(dtdir)  # Set working directory
source(srcfile) #  Source the deleteome_tools.R file
# --------------------------- #


alldata <- getDeleteomeExpData() # Load Deleteome expression data

# This script will find mutant strains transcriptionally similar to the strains specified here.
# Use getAllStrainNames(alldata) to output a list of all available Deleteome mutant strains that can be processed
mutantnames <- c("nup170", "pex10")  #A character string, vector or list object can be used here

Mthresh <- 0.0 # Absolute log2 fold-change cutoff for identifying differentially-expressed genes in deletion strain
pDEGs <- 0.05 # P-value cutoff for identifying differentially-expressed genes in deletion strain
pEnrich <- 0.05  # P-value cutoff for identifying statistically-significant correlation tests
qthresh <- 0.05 # quantile cutoff for down-selecting similar deletion strains

# Find similar strains in Deleteome
for(mutantname in mutantnames){  # Iterate through query deletion strains
  
  # Find similar strains in Deleteome
  similarStrains <- getSimilarStrainsByEnrichment( delData = alldata,
                                                   mutant=mutantname, 
                                                   minAbsLog2FC = Mthresh, 
                                                   pDEGs = pDEGs,
                                                   pEnrich = pEnrich,
                                                   quantileCutoff = qthresh
                                                  )
  
  # Show correlated strains
  message("\nSimilar deletion strains:")
  print(similarStrains)
  
  # Generate heatmap for significantly correlated deleteome strains. First get a data frame that contains all
  # the log2 fold-change and p-value info for the query strain's set of differentially-expressed genes...
  mutantProfile <- getProfileForDeletion(delData = alldata, 
                                         mutant = mutantname, 
                                         minAbsLog2FC = Mthresh, 
                                         pDEGs = pDEGs)

  # ...then generate a heatmap showing expression values for the query strain and its significantly similar
  # deleteome strains. Rows are the deletion strain's differentially-expressed genes. Columns are the deletion
  # strains found to be similar to the query strain.
  hm1 <- makeHeatmapDeleteomeMatches(mutantname = mutantname, 
                                     mutantProfile = mutantProfile, 
                                     selectedConditions = similarStrains, 
                                     fileprefix = "SignatureEnrich_matches", 
                                     titledesc="transcriptional similarity (signature enrichment)", 
                                     MthreshForTitle = Mthresh, 
                                     pDEGsForTitle = pDEGs,
                                     pMatchesForTitle = pEnrich,
                                     quantileForTitle = qthresh,
                                     imagewidth = 5000,
                                     printToFile = T)

  # Make the same heatmap as above but only show expression values for SUBTELOMERIC genes. This demonstrates
  # how users can make heatmaps for comparing subtelomeric gene expression patterns between a query strain
  # and other Deleteome strains.
  # NOTE: if there are less than 2 subtelomeric genes in the strain's signature
  # this will report an error indicating insufficient rows or columns for the heatmap
  hm2 <- makeHeatmapDeleteomeMatches(mutantname = mutantname, 
                                     mutantProfile = mutantProfile, 
                                     selectedConditions = similarStrains, 
                                     fileprefix = "SignatureEnrich_matches_SUBTELO", 
                                     titledesc="transcriptional similarity (signature enrichment)", 
                                     MthreshForTitle = Mthresh, 
                                     pDEGsForTitle = pDEGs,
                                     pMatchesForTitle = pEnrich,
                                     quantileForTitle = qthresh,
                                     subteloGenesOnly = T, 
                                     imagewidth = 5000,
                                     printToFile = T)

  # Generate a heatmap as in the first heatmap examplebut show query strain and manually-selected strains 
  # in columns. This demonstrates how users can compare expression values of a query strain's 
  # differentially-expressed genes to corresponding values in any arbitrary Deleteome strain.
  hm3 <- makeHeatmapDeleteomeMatches(mutantname = mutantname, 
                                     mutantProfile = mutantProfile, 
                                     selectedConditions = c("hmo1","rif1","sir4","ctf8","ctf18","dcc1"), 
                                     fileprefix = "Specific_mutants", 
                                     titledesc="manual selection", 
                                     MthreshForTitle = Mthresh, 
                                     pDEGsForTitle = pDEGs,
                                     pMatchesForTitle = pEnrich,
                                     quantileForTitle = qthresh,
                                     imagewidth = 2550,
                                     printToFile = T)

  # Make a mountain lake plot showing distance of a strain's differentially-expressed genes from the closest
  # telomere. This also reports (in the console) subtelomeric enrichment p-values for the mutant's up- and 
  # down-regulated genes.
  makeGenomicPositionHistogram(delData = alldata,
                               mutant = mutantname, 
                               Mthresh = Mthresh, 
                               pDEGs = pDEGs,
                               ymax = 40,
                               printToFile = T)

  # Run GO enrichment on genes deleted in similar deletion strains
  GOpadjcutoff = 0.1
  GO <- doGOenrichmentOnDeleteomeMatches(delData = alldata, 
                                         genes = similarStrains, 
                                         padjthresh = GOpadjcutoff)
  
  # If GO enrichment was performed, output results to a file
  if(! is.null(GO)){
    GOoutputfile <- paste0(dtdir,"/output/GO_enrichment/", mutantname, "_GOresults_SignatureEnrich_GOpadj", GOpadjcutoff,".tsv")
    write.table(GO, file=GOoutputfile, sep="\t", quote = F, row.names = F)
    message("Wrote GO enrichment results to ", GOoutputfile)
  }
}
