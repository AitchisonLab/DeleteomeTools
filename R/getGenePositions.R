#' Get genomic locations for genes
#'
#' @param includeMito Whether to include mitochondrial genes in the output
#'
#' @returns Data frame with distances of genes (in bases) from closest telomere and from centromere

getGenePositions <- function(includeMito=F){

  # Create table that associates gene ID with chromosome, as well as start and end positions on chromosome
  gpfile <- system.file("extdata", "genePositions.csv", package = "DeleteomeTools")
  centfile <- system.file("extdata", "YeastMine_CentromerePositions.tsv", package = "DeleteomeTools")
  chrfile <- system.file("extdata", "chrNameLength.txt", package = "DeleteomeTools")

  genePositions <- read.table(gpfile, sep=',',quote="",stringsAsFactors=F, header=T)
  centPositions <- read.table(centfile, sep="\t", quote="", stringsAsFactors = F, header=T)

  # Remove mitochondrial genes, if requested
  if( ! includeMito) genePositions <- genePositions[! grepl("Mito", genePositions$Chr),]

  # Column name for distances from telomere
  dftname = "dist_from_telo"

  # Columne name for distances from centromere
  dfcname <- "dist_from_cent"

  genePositions[,dftname] <- NA
  genePositions[,dfcname] <- NA

  # Create table of chromosome names and lengths
  chrNameLengths <- read.table(chrfile, header=FALSE, sep="\t")

  # Create table that associates gene with distance from telomere
  message("Computing gene distances from telomere and centromere")

  for(geneid in genePositions[,1]){

    rowNum = which(genePositions[,1] == geneid)

    # Get chromosome containing the gene. Convert from data frame element to character array.
    chrName = genePositions[rowNum, "Chr"]

    # Get gene start position
    geneStartPos = genePositions[rowNum, "Start"]
    geneEndPos = genePositions[rowNum, "End"]

    # Look up chromosome length
    chrLength = chrNameLengths[which(chrNameLengths[,1] == chrName),2]

    distFromCent <- NA

    if(chrName!="Mito"){
      # Look up position of centromere on chromosome
      chrCentStart <- as.integer(centPositions[centPositions$Description==paste0("Chromosome ",chrName, " centromere"),"Start"])
      chrCentEnd <- as.integer(centPositions[centPositions$Description==paste0("Chromosome ",chrName, " centromere"),"End"])

      # Compute distance from centromere
      if(geneEndPos < chrCentStart){
        distFromCent <- chrCentStart-geneEndPos
      }
      else if(geneStartPos > chrCentEnd){
        distFromCent <- geneStartPos-chrCentEnd
      }
      else{
        message("\nERROR: ",geneid," is neither ahead of nor behind centromere") # would only be thrown if a gene overlaps the centromere
        return(NULL)
      }
    }

    # Compute distance from telomere
    distFromChrStart = (geneStartPos - 1)
    distFromChrEnd = (chrLength - geneStartPos)
    distFromTelo = min(distFromChrStart, distFromChrEnd)
    genePositions[rowNum, dftname] <- distFromTelo
    genePositions[rowNum, dfcname] <- distFromCent
  }
  return(genePositions) # If including mitochondrial genes, then distance from centromere is included as NA
}
