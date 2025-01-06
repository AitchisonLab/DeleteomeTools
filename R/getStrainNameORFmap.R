
#' Getting mappings between Deleteome strain names and their open reading frames (ORFs)
#' @returns A data frame with three columns (Strain, ORF, MappingType)
#'
#' @export
getStrainNameORFmap <- function(){

  require(org.Sc.sgd.db)

  allstrainnames <- toupper(getAllStrainNames())

  allgenenames <- keys(org.Sc.sgd.db, keytype = "GENENAME")
  allorfs <- keys(org.Sc.sgd.db, keytype = "ORF")
  allaliases <- keys(org.Sc.sgd.db, keytype = "ALIAS")

  genenameorfmap <- suppressMessages(AnnotationDbi::select(org.Sc.sgd.db, keys = allgenenames, keytype = "GENENAME", columns = c("GENENAME", "ORF")))
  aliasorfmap <- suppressMessages(AnnotationDbi::select(org.Sc.sgd.db, keys = allstrainnames, keytype = "ALIAS", columns = c("ALIAS", "ORF")))

  strainorfmap <- data.frame(Strain=allstrainnames, ORF=NA)

  for(strain in allstrainnames){

    mapmethod <- ""
    orf <- c()

    # If strain name is in gene names, look up its ORF
    if(strain %in% genenameorfmap$GENENAME){
      orf <- na.omit(genenameorfmap[genenameorfmap$GENENAME==strain,"ORF"][1])  # 1-to-1 mapping exists for each strain name
      mapmethod <- "Gene name to ORF"
    }
    else if(strain %in% allorfs){  # Keep the gene name as-is
      orf <- strain
      mapmethod <- "ORF to ORF"
    }
    # If strain is not in gene names but it is in aliases, get ORF
    else if(strain %in% aliasorfmap$ALIAS){
      orf <- na.omit(aliasorfmap[aliasorfmap$ALIAS==strain,"ORF"])  # 1-to-many mapping exists here
      mapmethod <- "Gene alias to ORF"
    }

    if(length(orf) == 1){
      strainorfmap[strainorfmap$Strain==strain, "ORF"] <- orf
      strainorfmap[strainorfmap$Strain==strain, "MappingType"] <- mapmethod
    }
    else if(length(orf) > 1){
      message("\nERROR: Found multiple ORFs for strain ", strain)
      return(invisible(NULL))
    }
  }

  # Need to manually assign some deleteome strain names to their ORFs
  strainorfmap[strainorfmap$Strain=="ARG5_6","ORF"] <- "YER069W"
  strainorfmap[strainorfmap$Strain=="CYCC","ORF"] <- "YNL025C"  # SSN8
  strainorfmap[strainorfmap$Strain=="MF_ALPHA_1","ORF"] <- "YPL187W"
  strainorfmap[strainorfmap$Strain=="MF_ALPHA_2","ORF"] <- "YGL089C"
  strainorfmap[strainorfmap$Strain=="YAL044W_A","ORF"] <- "YAL044W-A"
  strainorfmap[strainorfmap$Strain=="YDR034W_B","ORF"] <- "YDR034W-B"
  strainorfmap[strainorfmap$Strain=="YIL014C_A","ORF"] <- "YIL014C-A"
  strainorfmap[strainorfmap$Strain=="YOL086W_A","ORF"] <- "YOL086W-A"

  return(strainorfmap)
}
