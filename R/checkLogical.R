#' Function to check validity of logical parameters
#' @param parName Parameter name
#' @param parValue Input parameter value

checkLogical <- function(parName, parValue){

  if( ! is.logical(parValue)){
    message("ERROR: Please enter TRUE or FALSE for the parameter ", parName, ".")
    return(F)
  }

  return(T)
}
