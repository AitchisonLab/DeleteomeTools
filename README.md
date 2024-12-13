# Deleteome-Tools

The Deleteome compendium, a publicly-available collection of ~1,500 transcriptomes from single-gene deletion yeast strains, is described in 

[Kemmeren P, et al. Large-scale genetic perturbations reveal regulatory networks and an abundance of gene-specific repressors. Cell. 2014 Apr 24;157(3):740-52.](https://pubmed.ncbi.nlm.nih.gov/24766815/)
and can be downloaded from [https://deleteome.holstegelab.nl](https://deleteome.holstegelab.nl)


Our Deleteome-Tools software was originally developed to identify gene products that are functionally associated with the nucleoporin NUP170.
We used the data in the Deleteome to identify genes that, when deleted, alter the yeast transcriptome in a way that is similar to the alterations caused by a NUP170 deletion.
This work is described in 

[Kumar, et al. Nuclear pore complexes mediate subtelomeric gene silencing by regulating PCNA levels on chromatin. J Cell Biol (2023) 222 (9): e202207060](https://doi.org/10.1083/jcb.202207060)

For any single-gene deletion profiled in the Deleteome, our software can be used to identify other single-gene deletion strains that are transcriptomically similar. It offers two methods for assessing similarity between transcriptomic profiles. The first quantifies similarity by conducting correlation tests on log2 fold-change values of transcriptomic profiles. The second method employs hypergeometric tests to determine whether the set of significantly altered genes shared among transcriptomic profiles occurs more frequently than expected by chance. 

The first approach considers the magnitude of expression changes in differentially-expressed genes, while the second approach focuses on whether a gene is differentially expressed or not, based on user-defined thresholds. 

In our studies with NUP170, we have observed that the two approaches yield similar, complementary results. The correlation-based method tends to be the more conservative option.

## Getting started in RStudio

* Clone this repository to your location of choice.
* Open the _get_similar_mutants_by_correlation.R_ script in R/RStudio and run it (e.g. by using "Source" in RStudio). This will identify deletion strains in the Deleteome that are similar to an input list of "query" deletion strains using the correlation-based methodology described above. (Run _get_similar_mutants_by_enrichment.R_ to perform the same analysis using the hypergeometric-based alternative.)
* A table showing the ranked list of similar deletion strains is saved to the "output/mutant_similarity" folder within the repository.
* The example script will also generate heatmaps showing gene expression values across the similar Deleteome strains it identifies. These are saved in the "output/heatmaps" folder. 
* A plot showing numbers of significantly up- and down-regulated genes according to their distance from telomeres will also be shown. We have used these "mountain lake" plots to identify and illustrate subtelomeric silencing defects in the NUP170 deletion strain as well as other Deleteome strains.
* The script will also perform a Gene Ontology (GO) enrichment analysis on the collected set of genes deleted among the similar deletion strains it finds. For example, the _get_similar_mutants_by_correlation.R_ script identifies 39 strains similar to the NUP170 deletion strain. The GO analysis is performed on the collected set of 39 genes deleted among those strains. The list of all genes deleted across the Deleteome is used as the background for these enrichment tests. GO analysis results are saved in the "output/GO_enrichment" folder.

To change the strains analyzed in the example scripts, change the values of the "mutantnames" variable.
For example, to find deletion strains similar to the NUP188 and PEX5 deletion strains, set the variable to c("nup188", "pex10"):
```
mutantnames <- c("nup188", "pex5")
```
Users can view the full list of genes that have associated deletion strains in the Deleteome using the following code:

```
alldata <- getDeleteomeExpData() # Loads Deleteome data
getAllStrainNames(alldata) # Prints the gene deleted for each Deleteome strain
```

## Getting started in the base R environment

* Clone this repository to your location of choice.
* Base R users will need to specify the path to the downloaded Deleteome-Tools repository and set it as their working directory.
* If running the example _get_similar_mutants_by_correlation.R_ script, edit line 5 so the ```dtdir``` variable is set to the path to the downloaded repository folder.
```
dtdir <- "[...]" 
```
Where [...] is replaced with the full path to the repository. Users can then source the _get_similar_mutants_by_correlation.R_ script from within base R to run it.
```
source("FULL PATH TO EXAMPLE SCRIPT")
```

## Step-by-step console usage
The scripts _get_similar_mutants_by_correlation.R_ and _get_similar_mutants_by_enrichment.R_ in this repository contain examples of the commands used to perform strain similarity matching,  visualize results, and run GO enrichment analyses. Here we detail each major function illustrated in the examples in turn.
* First, set the working directory to the repository folder
```
dtdir <- "ENTER YOUR PATH TO REPOSITORY FOLDER"
setwd(dtdir)  # Set working directory
```
* Load the Deleteome-Tools codebase
```
source(paste0(dtdir, "/deleteome_tools.R")) #  Source the deleteome_tools.R file
```
* Load the Deleteome expression data
```
alldata <- getDeleteomeExpData() # Load Deleteome expression data
```

* Use reciprocal correlation to find strains in the Deleteome that are transcriptionally similar to a query gene/strain (nup170 in this example). Default values for cut-off parameters are shown.
```
similarStrains <- getSimilarStrainsByReciprocalCorrelation( delData = alldata,    # Deleteome expression data object
                                                            mutant = "nup170",    # Query gene/strain
                                                            minAbsLog2FC = 0.0,   # Absolute log2 fold-change cutoff for identifying differentially-expressed genes in deletion strain
                                                            pDEGs = 0.05,         # P-value cutoff for identifying differentially-expressed genes in deletion
                                                            pCor = 0.05,          # P-value cutoff for identifying statistically-significant correlation tests
                                                            quantileCutoff = 0.1  # Quantile cutoff for down-selecting significantly similar deletion strains
                                                          )
``` 
* Alternatively, signature enrichment can be used to find similar strains
```
similarStrainsErich <- getSimilarStrainsByEnrichment(  delData = alldata,     # Deleteome expression data object
                                                       mutant = "nup170",     # Query gene/strain
                                                       minAbsLog2FC = 0.0,    # Absolute log2 fold-change cutoff for identifying differentially-expressed genes in deletion strain
                                                       pDEGs = 0.05,          # P-value cutoff for identifying differentially-expressed genes in deletion
                                                       pEnrich = 0.05,        # P-value cutoff for identifying statistically-significant correlation tests
                                                       quantileCutoff = 0.05  # Quantile cutoff for down-selecting significantly similar deletion strains
                                                       )
```
* Users can print the list of available query genes/strains using this command
```
print(getAllStrainNames(alldata))
```
* Make a heatmap showing gene expression values for the query strain's differentially-expressed genes and corresponding values in similar strains
```
# First collect expression values for the query strain's differentially-expressed genes
mutantProfile <- getProfileForDeletion( delData = alldata,   # Deleteome expression data object
                                        mutant = "nup170",   # A Deleteome strain
                                        minAbsLog2FC = 0.0,  # Absolute log2 fold-change cutoff for identifying differentially-expressed genes in deletion strain
                                        pDEGs = 0.05         # P-value cutoff for identifying differentially-expressed genes in deletion strain
                                      )

# Make the heatmap
hm1 <- makeHeatmapDeleteomeMatches( mutantname = "nup170",  # Strain name to display in heatmap title
                                    mutantProfile = mutantProfile,  # Expression values for strain's differentially-expressed genes
                                    selectedConditions = similarStrains,  # THe set of other strains to show in heatmap (character vector)
                                    fileprefix = "Corr_matches",  # A filename prefix to use when writing heatmap to a file
                                    titledesc = "transcriptional similarity (reciprocal correlation)",  # Description of comparison method for use in heatmap title
                                    MthreshForTitle = 0.0, # Absolute log2 fold-change cutoff used in obtaining mutantProfile. Used in heatmap title.
                                    pDEGsForTitle = 0.05,  # P-value cutoff used in obtaining mutantProfile. Used in heatmap title.
                                    pMatchesForTitle = 0.05,  # If the selectedConditions are included by virtue of similarity analysis, the P-value cutoff used in that analysis. Used in heatmap title.
                                    quantileForTitle = 0.1,  # If the selectedConditions are included by virtue of similarity analysis, the P-value cutoff used in that analysis. Used in heatmap title.
                                    imagewidth = 5000,  # If wrwiting to file, width of image in pixels
                                    printToFile = T  # Whether to write the heatmap to a file or show in a new window
                                    )
```
* Make a Mountain Lake plot for the query gene/strain that shows the genomic position of the strain's differentially-expressed genes relative to closest telomere
```
makeGenomicPositionHistogram(  delData = alldata,  # Deleteome expression data object
                               mutant = "nup170",  # Deletion strain to visualize
                               Mthresh = 0.0,      # Absolute log2 fold-change cutoff for identifying differentially-expressed genes in deletion strain
                               pDEGs = 0.05,       # P-value cutoff for identifying differentially-expressed genes in deletion strain
                               ymax = 40,          # Maximum value for Y-axis
                               printToFile = T     # Whether to write the heatmap to a file or show in a new window
                            )
```

* Perform GO enrichment to functionally profile the genes deleted in strains similar to the query strain, thereby predicting function of query gene
```
GOresults <- doGOenrichmentOnDeleteomeMatches( delData = alldata,         # Deleteome expression data object
                                  genes = similarStrains,    # Set of genes to test for GO enrichment (the genes deleted in each strain found to be similar to the query strain)
                                  padjthresh = 0.1           # FDR-adjusted P-value cutoff for significant enrichment. Only enrichment tests meeting this cutoff are returned.
                                )
```

## Dependencies
Users will need to have the following R packages installed to run the example scripts:
* [org.Sc.sgd.db](https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html)
* [gplots](https://www.rdocumentation.org/packages/gplots/versions/3.1.3)
* [ggplot2](https://ggplot2.tidyverse.org)
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
