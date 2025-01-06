#' Add quantile level to data frame of deletion similarity results. Expects the input data frame to have a "Pvalue.FDR" column
#' @param results Correlation or enrichment test P-values
#' @returns results object with quantile level added as data frame column
addQuantileLevel <- function(results=data.frame()){

  if(dim(results)[1] > 0){
    cdf <- ecdf(results$Pvalue.FDR)
    results$Pvalue.FDR.quantile <- cdf(results$Pvalue.FDR)
  }
  else results$Pvalue.FDR.quantile <- numeric(0)

  return(results)
}
