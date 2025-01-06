#' Get names of all deletion strains in the Deleteome
#'
#' @return A character vector of Deleteome strain names
#' @export
#'
#' @examples
#' # getAllStrainNames()
#'
getAllStrainNames <- function(){

  colheads <- gsub("_vs.*","",names(deldata))
  colheads <- gsub("_del","",colheads)
  colheads <- colheads[ ! colheads %in% c("systematicName","geneSymbol")] # exclude systematicName and geneSymbol columns

  return(sort(unique(colheads)))
}
