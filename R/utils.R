#' Apply the functions in the list to the specified column
#' @param fun_list a list, whose elements are functions to be applied to covariates specified by \code{col_index} one by one.
#' @param col_index a vector, which specifies the indexes of covariates corresponding to \code{fun_list} in the data matrix.
#' @noRd
apply_fun <- function(x, fun_list, col_index){
  y <- 0
  for (i in 1:length(fun_list)) {
    y <- y + fun_list[[i]](x[, col_index[i]])
  }
  return(y)
}
