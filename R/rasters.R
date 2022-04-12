#' Create dummy rasters based on raster to match
#'
#' @param rtm `RasterLayer` to use as template.
#' @param gaus logical. If `TRUE` (default), creates raster using Gaussian map.
#'             If `FALSE`, creates a raster with all pixel values equal to `1`.
#' @param rownumb number of rows for the new raster, if `rtm` not provided.
#' @param colnumb number of columns for the new raster, if `rtm` not provided.
#' @param xmn min x extent, if `rtm` not provided.
#' @param xmx max x extent, if `rtm` not provided.
#' @param ymn min y extent, if `rtm` not provided.
#' @param ymx max y extent, if `rtm` not provided.
#' @param sc spatial scale, passed to `gaussMap`.
#' @param vr spatial variance, passed to `gaussMap`.
#'
#' @return `RasterLayer`
#'
#' @author Louis-Etienne Robert and Alex Chubaty
#' @export
#' @importFrom raster extent extent<- raster
#' @importFrom SpaDES.tools gaussMap
createRaster <- function(rtm = NULL, gaus = TRUE, rownumb = 20, colnumb = 20,
                         xmn = -10, xmx = 10, ymn = -10, ymx = 10, sc= 100, vr = 5) {
  if (is.null(rtm)) {
    r <- raster(nrows = rownumb, ncols = colnumb, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
    if (isTRUE(gaus)) {
      r <- abs(round(gaussMap(r, scale = sc, var = vr)))
    } else {
      r[] <- 1
    }
  } else {
    r <- rtm
    if (isTRUE(gaus)) {
      r <- abs(round(gaussMap(r, scale = sc, var = vr)))
    } else {
      r[!is.na(r[])] <- 1
    }
  }

  return(r)
}
