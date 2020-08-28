#' Process the formula into parts we can use.
#'
#' @param formula
#' @return a list of formulas.
list_formulas <- function(formula, env){
  fstr <- strsplit(as.character(f), split = "\\|")[2:3]
  flist <- list(as.formula(fstr[[1]][1], env = env),
    as.formula(paste0(fstr[[1]][2], " ~ ", fstr[[2]]), env = env))
  flist
}
