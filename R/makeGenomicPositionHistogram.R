#' Make a mountain lake plot: a histogram showing position of differentially-expressed genes relative to nearest telomere.
#'
#' @param strain Name of Deleteome mutant strain to analyze
#' @param genePositions Table of genes and genomic positions (obtained using getGenePositions())
#' @param relativeTo Valid options are "telomere" or "centromere". If "telomere" specified, distance to closest telomere (right or left) is plotted.
#' @param rangeInKB Window (in kilobases) for defining a gene as either in the telomeric or centromeric region
#' @param Mthresh Log2 fold-change threshold for defining differentially-expressed genes
#' @param pDEGs P-value threshold for defining differentially-expressed genes
#' @param xMax X-axis maximum (in kilobases). If NA, uses maximum distance to closest telomere or to the centromere among all genes.
#' @param yMax Y-axis maximum
#' @param upColor Color to use for upregulated genes
#' @param downColor Color to use for downregulated genes
#' @param showUpDownLabels Whether to inclue "upregulated" and "downregulated" labels in plot
#' @param outputDir If writing to file, output directory to save image
#' @param printToFile Whether to write heatmap to a file or show in new window.
#'
#' @export
#'
# Mountain lake plots
makeGenomicPositionHistogram <- function(strain=NULL,
                                         genePositions=NULL,
                                         relativeTo="telomere",
                                         rangeInKB=25,
                                         Mthresh=0,
                                         pDEGs=0.05,
                                         xMax=NA,
                                         yMax=50,
                                         upColor="#d53e4f",
                                         downColor="#3288bd",
                                         showUpDownLabels=T,
                                         outputDir=NA,
                                         printToFile=T
){

  if(is.null(strain)) message("Cannot produce genomic position histogram: mutant strain ID is invalid.")

  if( ! relativeTo %in% c("telomere", "centromere")){
    message("\nERROR: please use either \"telomere\" or \"centromere\" for the \"relativeTo\" parameter. ", relativeTo, " is not valid." )
    return(invisible(NULL))
  }

  if( ! is.numeric(rangeInKB)){
    message("\nERROR: Input value for rangeInKB parameter must be numeric.")
    return()
  }

  outputDir <- file.path(outputDir)
  if(is.na(outputDir) & printToFile){
    message("ERROR: Please specify the directory in which to save plot")
    return(invisible(NULL))
  }
  else if( ! dir.exists(outputDir) & printToFile){
    message("ERROR: The specified output directory ", outputDir, " does not exist.")
    return(invisible(NULL))
  }

  message("Making genomic position histogram for mutant strain ", strain)

  require(ggplot2)


  if(is.null(genePositions)){
    message("Loading gene positions across deleteome. Excluding mitochondrial genes.")
    genePositions <- getGenePositions()
  }

  if(relativeTo=="telomere") relativeToCol <- "dist_from_telo"
  else if(relativeTo=="centromere") relativeToCol <- "dist_from_cent"

  # Use the following to show VdVFig3 plots of selected strains
  allDF <- genePositions[,relativeToCol]/1000 # get all distance from telomere numbers; filter out mitochondrial genes? If so: genePositions$Chr != 'Mito'

  strainSignature <- getStrainSignature(strain, Mthresh, pDEGs, consoleMessages = F)

  signatureUP <- strainSignature[strainSignature[,3] > Mthresh, ]
  signatureDOWN <- strainSignature[strainSignature[,3] <= -Mthresh, ]
  signatureALL <- getStrainSignature(strain, 0, 1.0, consoleMessages = F)

  upDF <- genePositions[genePositions$Geneid %in% signatureUP$systematicName, relativeToCol]
  downDF <- genePositions[genePositions$Geneid %in% signatureDOWN$systematicName, relativeToCol]

  upDF <- na.omit(upDF/1000)
  downDF <- na.omit(downDF/1000)

  message("Number of significantly upregulated genes: ", length(upDF))
  message("Number of significantly downregulated genes: ", length(downDF))

  mybreaks = seq(from=0, to=max(na.omit(genePositions[,relativeToCol]))/1000, by=5)

  if(relativeTo=="telomere") mybreaks <- c(mybreaks,770) # max dist from telo in genePositions is 765.197 kb. Need one more breakpoint for histogram to include it.
  else mybreaks <- c(mybreaks, 1085) # max dist from centromere in genePositions is 1081.042 kb. Need one more breakpoint for histogram to include it.

  regionUPhyperp <- genomicRegionEnrichment(genePositions,
                                            systematicNamesBG = signatureALL[signatureALL$systematicName %in% genePositions$Geneid, "systematicName"],
                                            systematicNames = signatureUP[signatureUP$systematicName %in% genePositions$Geneid, "systematicName"],
                                            relativeTo = relativeTo, rangeInKB = rangeInKB)
  message("Hypergeometric test P-value for upregulated subtelomeric genes: ", regionUPhyperp)

  regionDOWNhyperp <- genomicRegionEnrichment(genePositions, systematicNamesBG=signatureALL[signatureALL$systematicName %in% genePositions$Geneid, "systematicName"],
                                              systematicNames = signatureDOWN[signatureDOWN$systematicName %in% genePositions$Geneid, "systematicName"],
                                              relativeTo = relativeTo, rangeInKB = rangeInKB)
  message("Hypergeometric test P-value for downregulated subtelomeric genes: ", regionDOWNhyperp)

  # ggplot2 style
  upDFgg <- as.data.frame(upDF)
  names(upDFgg) <- "dist"
  downDFgg <- as.data.frame(downDF)
  names(downDFgg) <- "dist"
  updownDF <- rbind()
  allDFgg <- as.data.frame(allDF)
  names(allDFgg) <- "dist"

  relativeToLabel <- relativeTo
  if(relativeTo == "telomere") relativeToLabel <- "closest telomere"

  if(is.na(xMax)) xMax <- max(mybreaks)

  # Create plot with ggplot
  p <- ggplot() +
    geom_histogram(data = allDFgg[allDFgg$dist <= xMax, , drop = F], aes(x = dist, y = after_stat(count)/3), fill="gray", binwidth = 5, boundary=-5) +
    geom_histogram(data = allDFgg[allDFgg$dist <= xMax, , drop = F ], aes(x = dist, y = -after_stat(count)/3), fill="gray", binwidth = 5, boundary=-5) +
    geom_histogram(data = upDFgg[upDFgg$dist <= xMax, , drop = F], aes(x = dist, y = after_stat(count)), fill=upColor, binwidth = 5, boundary=-5) +
    geom_histogram(data = downDFgg[downDFgg$dist <= xMax, , drop = F], aes(x = dist, y = -after_stat(count)), fill= downColor, binwidth = 5, boundary=-5) +
    scale_x_continuous(breaks=seq(0, xMax, by=25)) +
    scale_y_continuous(labels = abs) +
    theme(plot.title = element_text(size=25, face="bold"), axis.text=element_text(size=11),
          axis.title=element_text(size = 18,face="plain"), axis.text.x = element_text(angle = 45, hjust=1)) +
    xlab(paste0("Distance from ",relativeToLabel," (kb)")) +
    ylab(paste0("Number of ORFs")) +
    ggtitle(bquote(italic(.(strain)*Delta)))

  # Show the "upregulated" and "downregulated" labels if desired
  if(showUpDownLabels)
    p <- p + annotate("text", x=(xMax-(xMax/4)), y=(yMax-(yMax/4)), label= "upregulated", color=upColor, size = 6) +
    annotate("text", x=(xMax-(xMax/4)), y=(-yMax+(yMax/4)), label= "downregulated", color=downColor, size = 6)

  if(! printToFile) dev.new(width=9,height=5.5, noRStudioGD = T)
  else{
    mlfile <- paste0(outputDir, "/",strain,"_MountainLake_l2FC", Mthresh,"_pDEGs", pDEGs, "_", relativeTo, "_range", rangeInKB, "KB.png")
    png(filename = mlfile,width=9, height=5.5, res=300, units="in")
  }

  print(p)

  if(printToFile){
    dev.off()
    message("\nPrinted mountain lake plot to ", mlfile)
  }
}

#' Make a mountain lake plot: a histogram showing position of differentially-expressed genes relative to nearest telomere.
#'
#' @param strain Name of Deleteome mutant strain to analyze
#' @param genePositions Table of genes and genomic positions (obtained using getGenePositions())
#' @param relativeTo Valid options are "telomere" or "centromere". If "telomere" specified, distance to closest telomere (right or left) is plotted.
#' @param rangeInKB Window (in kilobases) for defining a gene as either in the telomeric or centromeric region
#' @param Mthresh Log2 fold-change threshold for defining differentially-expressed genes
#' @param pDEGs P-value threshold for defining differentially-expressed genes
#' @param xMax X-axis maximum (in kilobases). If NA, uses maximum distance to closest telomere or to the centromere among all genes.
#' @param yMax Y-axis maximum
#' @param upColor Color to use for upregulated genes
#' @param downColor Color to use for downregulated genes
#' @param showUpDownLabels Whether to inclue "upregulated" and "downregulated" labels in plot
#' @param outputDir If writing to file, output directory to save image
#' @param printToFile Whether to write heatmap to a file or show in new window.
#'
#' @export
#'
makeMountainLakePlot <- makeGenomicPositionHistogram
