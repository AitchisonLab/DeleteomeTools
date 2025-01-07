#' Function to check validity of numeric parameters
#' @param parName Parameter name
#' @param parValue Input parameter value
#' @param minValue Minimum valid value for parameter
#' @param maxValue Maximum valid value for parameter

checkNumeric <- function(parName = "", parValue = NULL, minValue = NULL, maxValue = NULL){

  if( ! is.numeric(parValue)){
    message("ERROR: expected numeric value for parameter ", parName, " but input was ", parValue)
    return(F)
  }
  else{
    if( ! is.null(minValue)){  # If min value specified, check if in range
      if(parValue < minValue){
      message("ERROR: value for parameter ", parName, " too low. The minimum allowed value is ", minValue)
      return(F)
      }
    }
    if( ! is.null(maxValue)){  # If max value specified, check if in range
      if(parValue > maxValue){
        message("ERROR: value for parameter ", parName, " too high. The maximum allowed value is ", maxValue)
        return(F)
      }
    }
  }
  return(T)
}
