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

## Getting started

* Clone this repository to your location of choice.
* Open the _get_similar_mutants_by_correlation.R_ script in R/RStudio and run it (e.g. by using "Source" in RStudio). This will identify deletion strains in the Deleteome that are similar to the NUP170 deletion strain using the correlation-based methodology described above. (Run _get_similar_mutants_by_hypergeometric.R_ to perform the same analysis using the hypergeometric-based alternative.)
* A table showing the ranked list of similar deletion strains is saved to the "output/mutant_similarity" folder within the repository.
* The example script will also generate heatmaps showing gene expression values across the similar Deleteome strains it identifies. These are saved in the "output/heatmaps" folder. 
* A plot showing numbers of significantly up- and down-regulated genes according to their distance from telomeres will also be shown. We have used these "mountain lake" plots to identify and illustrate subtelomeric silencing defects in the NUP170 deletion strain as well as other Deleteome strains.
* The script will also perform a Gene Ontology (GO) enrichment analysis on the collected set of genes deleted among the similar deletion strains it finds. For example, the _get_similar_mutants_by_correlation.R_ script identifies 39 strains similar to the NUP170 deletion strain. The GO analysis is performed on the collected set of 39 genes deleted among those strains. The list of all genes deleted across the Deleteome is used as the background for these enrichment tests. GO analysis results are saved in the "output/GO_enrichment" folder.

To change the strain in the example scripts from NUP170 to a gene of interest, change the value of the "mutantname" variable.
For example, to find deletion strains similar to the PEX10 deletion strain, set the variable to "pex10":
```
mutantname <- "pex10"
```
Users can view the full list of genes that have associated deletion strains in the Deleteome using the following code:

```
alldata <- getDeleteomeExpData() # Loads Deleteome data
getAllStrainNames(alldata) # Prints the gene deleted for each Deleteome strain
```

## Dependencies
Users will need to have the following R packages installed to run the example scripts:
* [org.Sc.sgd.db](https://bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html)
* [gplots](https://www.rdocumentation.org/packages/gplots/versions/3.1.3)
* [ggplot2](https://ggplot2.tidyverse.org)
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
