#' Calculate the mean of non-zero values
#'
#' Used as input for focal to handle edges padded with 0, and 0 values.
#'
#' @param x DESCRIPTION NEEEDED
#'
#' @return numeric mean value
#'
#' @author Louis-Etienne Robert
#' @export
meanNZ <- function(x) {
  a <- as.vector(x)
  if (sum(a) == 0) { #handling of neighborhood all 0
    out = 0
    return(out)
  } else {
    nonzero <- length(a[a > 0])
    out <- sum(a) / nonzero
  }

  return(out)
}
