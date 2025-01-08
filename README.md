# DeleteomeTools

The Deleteome compendium, a publicly-available collection of ~1,500 transcriptomes from single-gene deletion yeast strains, is described in 

[Kemmeren P, et al. Large-scale genetic perturbations reveal regulatory networks and an abundance of gene-specific repressors. Cell. 2014 Apr 24;157(3):740-52.](https://pubmed.ncbi.nlm.nih.gov/24766815/)
and can be downloaded from [https://deleteome.holstegelab.nl](https://deleteome.holstegelab.nl)


Our DeleteomeTools software was originally developed to identify gene products that are functionally associated with the nucleoporin NUP170.
We used the data in the Deleteome to identify genes that, when deleted, alter the yeast transcriptome in a way that is similar to the alterations caused by a NUP170 deletion.
This work is described in 

[Kumar, et al. Nuclear pore complexes mediate subtelomeric gene silencing by regulating PCNA levels on chromatin. J Cell Biol (2023) 222 (9): e202207060](https://doi.org/10.1083/jcb.202207060)

For any single-gene deletion profiled in the Deleteome, our software can be used to identify other single-gene deletion strains that are transcriptomically similar. It offers two methods for assessing similarity between transcriptomic profiles. The first quantifies similarity by conducting correlation tests on log2 fold-change values of transcriptomic profiles. The second method employs hypergeometric tests to determine whether the set of significantly altered genes shared among transcriptomic profiles occurs more frequently than expected by chance. 

The first approach considers the magnitude of expression changes in differentially-expressed genes, while the second approach focuses on whether or not a gene is differentially expressed and the direction of its expression change. 

In our studies with NUP170, we have observed that the two approaches yield similar, complementary results. The correlation-based method tends to be the more conservative option.

## Installation
Make sure the _devtools_ or _remotes_ package is installed in your R environment. To install, use
```
install.packages("devtools")
```
or
```
install.packages("remotes")
```
Then, install DeleteomeTools via this GitHub site
```
library(devtools)  # or library(remotes)
install_github("AitchisonLab/DeleteomeTools")
```
If needed, users can overwrite a previous installation using
```
install_github("AitchisonLab/DeleteomeTools", force = T)
```
## Example usage
Load the DeleteomeTools package
```
library(DeleteomeTools)
```
### Identifying similar deletion strains
DeleteomeTools was primarily developed to identify deletion strains in the Deleteome whose gene expression changes vs. wild-type are similar to a strain where a gene of interest was deleted. We refer to the gene of interest as the "query gene" and its associated deletion strain as the "query strain". 

For example, users can identify deletion strains similar to the _nup170_ deletion strain using
```
sim <- getSimilarStrainsByReciprocalCorrelation(strain = "nup170",
                                                outputDir = "[enter output folder path]")
```
The _sim_ variable will contain the names of deletion strains found to be transcriptionally similar to the _nup170_ deletion strain. 
A table showing the ranked list of similar deletion strains is saved to the folder specified in the _outputDir_ parameter.

To perform the same analysis using the hypergeometric-based alternative mentioned above, use
```
sim <- getSimilarStrainsByEnrichment(strain = "nup170",
                                     outputDir = "[enter output folder path]")
```
These similarity analysis functions also include adjustable parameters that set statistical cutoffs for identifying the query strain's set of differentially-expressed genes (DEGs) as well as its similar strains. The default values are based on a systematic analysis of DeleteomeTools' performance for predicting gene function. They are set to optimize the trade-off between functional prediction sensitivity and total number of functional predictions. Users can adjust these parameters to make statisitcal cutoffs more stringent or permissive as desired.

To see the full list of strain names that can be used as query strains/genes, use
```
getAllStrainNames()
```
To iteratively perform similarity analyses over multiple query genes and store results in a list, you can use
```
results <- list()  # List for storing results
for(astrain in c("nup170", "pex5")){
  sim <- getSimilarStrainsByReciprocalCorrelation(astrain, outputDir = "[enter output folder path]")
  results[[mystrain]] <- sim
}
```
### Gene function prediction
After identifying the set of deletion strains that are similar to a query strain, users can functionally profile the set of genes deleted in those strains using GO enrichment analysis. 
The GO biological processes, cellular components, and molecular functions that are significantly enriched in these analysis indicate the query strain may be associated with those processes, components and functions. By default, the list of all genes deleted across Deleteome strains is used as the background for the enrichment analysis.

Perform GO enrichment on the genes deleted in strains found to be similar to a query strain using
```
gores <- GOenrichmentOnStrains(sim)
print(gores)
```
The GOenrichmentOnStrains() function also includes an adjustable parameter _padjThresh_ that sets the signficance cutoff for enriched GO terms.

## Visualization features

### Gene expression heatmaps
It may be useful for users to visualize gene expression values for a query strain alongside values from other strains such as those found to be transcriptionally similar. DeleteomeTools allows users to generate customizeable heatmaps showing such comparisons.

The following code will generate a heatmap of gene expression values for the _nup170_ deletion strain's differentially expressed genes (DEGs, rows) alongside values from strains that are transcriptionally similar (columns). In this example, default values for statistical cutoff parameters (to select the query strain's DEGs and similar strains) are used.
```
# Specify query strain
querystrain <- "nup170"

# Set statistical cutoffs (default values for correlation-based similarity analyses are used here)
fccut <- 0       # Absolute log2 fold-change cutoff for identifying differentially-expressed genes in query strain
pdegcut <- 0.05  # P-value cutoff for identifying differentially-expressed genes in query strain
psimcut <- 0.05  # P-value cutoff for identifying statistically-significant correlation tests
qcut <- 0.1      # Quantile cutoff for down-selecting significantly similar strains

# Run similarity analysis using reciprocal correlation
sim <- getSimilarStrainsByReciprocalCorrelation( strain = querystrain,  # Query gene/strain
                                                 outputDir = "[enter output folder path]",  # Path to folder where analysis results are saved
                                                 minAbsLog2FC = fccut,  # Absolute log2 fold-change cutoff for identifying differentially-expressed genes in query strain
                                                 pDEGs = pdegcut,       # P-value cutoff for identifying differentially-expressed genes in query strain
                                                 pCor = psimcut,        # P-value cutoff for identifying statistically-significant correlation tests
                                                 quantileCutoff = qcut   # Quantile cutoff for down-selecting significantly similar strains
                                                 )

# Get dataframe of expression data for query strain
querysig <- getStrainSignature(strain = querystrain, minAbsLog2FC = fccut, pDEGs = pdegcut)

# Make heatmap of expression values
hm1 <- makeHeatmapDeleteomeMatches( strain = querystrain,         # Strain name to display in heatmap title
                                    strainSignature = querysig,   # Expression values for strain's differentially-expressed genes
                                    otherStrains = sim,           # The set of other strains to show in heatmap (character vector)
                                    filePrefix = "Corr_matches",  # A filename prefix to use when writing heatmap to a file
                                    titleDesc = "transcriptional similarity (reciprocal correlation)",  # Description of comparison method for use in heatmap title
                                    minAbsLog2FCforTitle = fccut, # Absolute log2 fold-change cutoff used in obtaining mutantProfile. Only used for heatmap title.
                                    pDEGsForTitle = pdegcut,      # P-value cutoff used in obtaining mutantProfile. Only used for heatmap title.
                                    pMatchesForTitle = psimcut,   # If the selectedConditions are included by virtue of similarity analysis, the P-value cutoff used in that analysis. Only used for heatmap title.
                                    quantileForTitle = qcut,      # If the selectedConditions are included by virtue of similarity analysis, the quantile cutoff used in that analysis. Only used for heatmap title.
                                    imageWidth = 5000,            # If writing to file, width of image in pixels
                                    outputDir="[enter output folder path]",  # If writing to file, folder in which to save image
                                    printToFile = T               # Whether to write the heatmap to a file or show in a new window
                                    )
```
Instead of writing the heatmap to an image file, users also have the option to just show the heatmap in a new window by setting the _pathToFile_ parameter to FALSE. Users can also limit the expression values in heatmap rows to subtelomeric genes only by setting the _subteloGenesOnly_ parameter to TRUE.

### Mountain lake plots
DeleteomeTools was developed in the context of research on the yeast nucleoporin NUP170 which mediates subtelomeric silencing. Therefore, the package also includes a  feature for visualizing subtelomeric silencing defects/enhancements. These "mountain lake" plots consist of histograms indicating the genomic position of a deletion strain's significantly up- and down-regulated genes relative to the closest telomere. Alternatively, their position relative to the centromere can be shown.

Default values for statistical cutoffs and axis limits are used in this example
```
makeMountainLakePlot(strain = "nup170",   # Deletion strain to visualize
                     minAbsLog2FC = 0,    # Absolute log2 fold-change cutoff for identifying differentially-expressed genes in deletion strain
                     pDEGs = 0.05,        # P-value cutoff for identifying differentially-expressed genes in deletion strain
                     yMax = 40,           # Y-axis maximum value
                     outputDir = "[enter path to output folder]",  # If writing to file, folder to save image
                     printToFile = T)     # Whether to save image as file or open in new window
```

## Dependencies
Users will need to have the following R packages installed to run the set examples above:
* [org.Sc.sgd.db](https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html)
* [gplots](https://www.rdocumentation.org/packages/gplots/versions/3.1.3)
* [ggplot2](https://ggplot2.tidyverse.org)
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
