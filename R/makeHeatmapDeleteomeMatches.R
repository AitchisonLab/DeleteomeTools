#' Make a heatmap showing expression data for a query strain's signature compared to other Deleteome strains
#'
#' @param strain A Deleteome strain whose expression signature will be shown. Only used for plot titles and file names.
#' @param strainSignature The transcriptional signature of the deletion strain
#' @param otherStrains List of additional Deleteome strains to include in heatmap
#' @param filePrefix String that can be added to beginning of heatmap file name. Appears after strain name.
#' @param titleDesc Rationale for inclusion of additional strains in heatmap. Used in title of heatmap.
#' @param minAbsLog2FCforTitle Absolute log2 fold-change threshold used to get the mutant strain's signature. Only used in heatmap title and file name.
#' @param pDEGsForTitle P-value threhshold used to get the mutant strain's signature (if applicable). Only used in heatmap title and file name.
#' @param pMatchesForTitle P-value threshold used to select similar mutant strains (if applicable)
#' @param quantileForTitle Quantile threshold used to select similar strains (if applicable). Only used in heatmap title and file name.
#' @param subteloGenesOnly Whether to only include subtelomeric genes in the rows of the heatmap
#' @param clusterColumns Whether to automatically cluster the columns of the heatmap
#' @param colFontSize Font size for column labels
#' @param showRowLabels Whether to show row labels
#' @param rowFontSize Font size for row labels
#' @param imageWidth If writing to file, the Width of the saved image in pixels
#' @param imageHeight If writing to file, the Width of the saved image in pixels
#' @param outputDir If writing to file, output directory to save image
#' @param printToFile Whether to write heatmap to a file or show in new window
#'
#' @export
#'
#' @returns Data matrix containing log2 fold-change expression data for heatmap
#'
#'
makeHeatmapDeleteomeMatches <- function(strain=NA,
                                        strainSignature=NA,
                                        otherStrains=NA,
                                        filePrefix=NA,
                                        titleDesc="",
                                        minAbsLog2FCforTitle=NA,
                                        pDEGsForTitle=NA,
                                        pMatchesForTitle=NA,
                                        quantileForTitle=NA,
                                        subteloGenesOnly=F,
                                        clusterColumns=T,
                                        colFontSize=1,
                                        showRowLabels=F,
                                        rowFontSize=0.275,
                                        imageWidth=3000,
                                        imageHeight=2400,
                                        outputDir=NA,
                                        printToFile=T){

  # Make heatmap comparing all significantly matched mutants to mutant of interest
  if( ! is.data.frame(strainSignature)){
    message("ERROR: Entry for parameter strainSignature must be a data frame. Use getStrainSignature() to create a valid data frame.")
  }

  if( ! is.character(otherStrains) | length(otherStrains)==0){
    message("\nERROR: Cannot make heatmap because the number of Deleteome strains to show alongside the ", strain, " deletion is 0. At least 1 is required.")
    return(invisible(NULL))
  }

  # Check logical parameters
  if( ! all( c(checkLogical("subteloGenesOnly", subteloGenesOnly),
               checkLogical("clusterColumns", clusterColumns),
               checkLogical("showRowLabels", showRowLabels),
               checkLogical("printToFile", printToFile)))){
    return(invisible(NULL))
  }

  if( ! all( c(checkNumeric("colFontSize", colFontSize, minValue = 0),
               checkNumeric("rowFontSize", rowFontSize, minValue = 0),
               checkNumeric("imageWidth", imageWidth, minValue = 1),
               checkNumeric("imageHeight", imageHeight, minValue = 1)))){
    return(invisible(NULL))
  }

  outputDir <- file.path(outputDir)
  if(is.na(outputDir) & printToFile){
    message("ERROR: Please specify the directory in which to save heatmap")
    return(invisible(NULL))
  }
  else if( ! dir.exists(outputDir) & printToFile){
    message("ERROR: The specified output directory ", outputDir, " does not exist.")
    return(invisible(NULL))
  }


  require(gplots)

  genesThatCanBeCompared <- strainSignature[ ! strainSignature$geneSymbol %in% toupper(otherStrains),]  # Omit strain from list of similar mutants

  rowtitleprefix <- "Genes"

  # If we're only considering subtelomeric genes...remove genes outside the region
  if(subteloGenesOnly){
    gps <- getGenePositions()
    subgps <- gps[gps$dist_from_telo<25000,"Geneid"]
    genesThatCanBeCompared <- genesThatCanBeCompared[genesThatCanBeCompared$systematicName %in% subgps,]
    rowtitleprefix <- "Subtelomeric genes"
  }

  if(dim(genesThatCanBeCompared)[1] <=2){
    message("\nERROR: Cannot make heatmap because the number of ", tolower(rowtitleprefix), " (row entries) to plot in the heatmap is ",
            dim(genesThatCanBeCompared)[1], ", but at least two are required.")
    return(invisible(NULL))
  }

  sigcorrheatmap <- data.frame(Gene=as.character(genesThatCanBeCompared$systematicName), stringsAsFactors=F)
  sigcorrheatmap[ , strain] <- as.numeric(strainSignature[strainSignature$systematicName %in% genesThatCanBeCompared$systematicName,3])

  # Collect transcriptional profiles for Deleteome strains that will be included in the heatmap
  for(cond in otherStrains){
    condprofileALL <- getStrainSignature(cond, 0, 1, consoleMessages = F)
    condhmdatadf <- condprofileALL[condprofileALL$systematicName %in% genesThatCanBeCompared$systematicName, c(2,3)]
    condhmdata <- condhmdatadf[,2]
    names(condhmdata) <- condhmdatadf[,1]
    sigcorrheatmap[,cond] <- as.numeric(condhmdata)
  }

  mybigmat <- data.matrix(sigcorrheatmap[,-1])
  rownames(mybigmat) <- sigcorrheatmap$Gene

  # Set colormap scale
  mybreaks <- seq(-1.5, 1.5, length.out=101)

  if(printToFile){
    hmfile <- paste0(outputDir, "/",strain,"_",filePrefix,"_Heatmap_l2FC",minAbsLog2FCforTitle,"_pDEGs",
                     pDEGsForTitle, "_pMatches",pMatchesForTitle,"_quantile",quantileForTitle,".png")
    png(file = hmfile, width=imageWidth, height = imageHeight, res=300, units="px")
  }
  else dev.new(width=10,height=8,noRStudioGD = TRUE)

  par(cex.main=0.75) ## Set size of main title and legend title

  rowlabels <- rownames(mybigmat)
  rightmar <- 5

  if( ! showRowLabels){
    rowlabels = rep("",dim(mybigmat)[1])
    rightmar <- 2
  }

  mycolnamesdeltaitalic <- colnames(mybigmat)
  mycolnamesdeltaitalic <- lapply(mycolnamesdeltaitalic, function(x) bquote(italic(.(x)*Delta)))
  names(mycolnamesdeltaitalic) <- colnames(mybigmat)
  mycolnamesdeltaitalic[[strain]] <- bquote(bolditalic(.(strain)*Delta))  # make the strain label bold

  mycolcolors <- rep("#4A4A4B", dim(mybigmat)[2])
  names(mycolcolors) <- colnames(mybigmat)
  mycolcolors[strain] <- "black"  # Make mutant name label black, not almost-black

  # Make the heatmap object
  clust <- heatmap.2(mybigmat,
                     Colv = clusterColumns, trace = "none", symbreaks = T,
                     xlab = "Deletion strain",
                     ylab = bquote(.(rowtitleprefix)*" in "*italic(.(strain)*Delta)*" signature"),
                     labRow = rowlabels,
                     labCol = as.expression(mycolnamesdeltaitalic),
                     colCol = mycolcolors,
                     srtCol = 45,
                     symm = F,
                     symkey = F,
                     scale = "none",
                     breaks = mybreaks,
                     margins = c(6, rightmar), # Makes sure the y-axis label is more flush with plot
                     cexCol = colFontSize,
                     cexRow = rowFontSize,
                     col = colorRampPalette(c("#0064FF","white","#DA0119"))(100),
                     key.title= "",
                     keysize = 1,
                     key.xlab = "Log2 fold-change vs. WT",
                     main = bquote(atop("Expression for "*italic(.(strain)*Delta)*" and other strains included by "*.(titleDesc),
                                        Cutoffs*": Abs. "*log[2]~fold*"-"*change*">="*.(minAbsLog2FCforTitle)*", "*DEG~italic("P-")*value*"="*.(pDEGsForTitle)*", "*Similarity~italic("P-")*value*"="*.(pMatchesForTitle)*", "*Quantile*"="*.(quantileForTitle)))
  )

  if(printToFile){
    dev.off()
    message("\nPrinted heatmap to ", hmfile)
  }

  return(mybigmat)
}
