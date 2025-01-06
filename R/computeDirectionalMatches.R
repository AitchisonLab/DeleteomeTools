#' Used in hypergeometric analysis to identify genes that changed expression in the same direction
#'
#' @param signature1 Strain 1 signature
#' @param signature2 Strain 2 signature
#'
#' @returns Data frame with information about overlapping genes in the two signatures


computeDirectionalMatches <- function(signature1=NA, signature2=NA){

  signature1 <- signature1[order(signature1$systematicName), ]
  signature2 <- signature2[order(signature2$systematicName), ]

  intsct <- intersect(signature1$systematicName, signature2$systematicName)

  signature1int <- signature1[signature1$systematicName %in% intsct,]
  signature2int <- signature2[signature2$systematicName %in% intsct,]

  signature1int$mvalsToMatch <- signature2int[, 3]
  signature1int$samedir <- signature1int$mvalsToMatch*signature1int[, 3]

  return(signature1int)
}
