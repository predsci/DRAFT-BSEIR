# This file provides extra details for the package man page

#' @details 
#' The primary function that a user will need is runbSEIR().  This function expects the user to input an incidence data.frame with columns 'date' and 'cases'.  For non-uniform reporting, the code assumes that the reporting period for the nth row runs from the date[n-1]+1 to date[n].  The 'date' column class should be "Date".  
#' @keywords internal
#' @useDynLib DRAFT, .registration = TRUE
"_PACKAGE"
#> [1] "_PACKAGE"

