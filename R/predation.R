#' Convert to enemy density
#'
#' Convert to density using the current budworm count (before scaling the budworm count).
#'
#' @param EnemyCount DESCRIPTION NEEDED
#' @param InsectCount DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom raster overlay
ConvertToEnemyDensity <- function(EnemyCount, InsectCount) {
  ## Initial equation:
  ##  sim$EnemyDensitySpring<-sim$EnemyCountSpring / sim$BudCountSpring

  ## Modified to take 0 budworm count into account
  EnemyDensity <- overlay(EnemyCount, InsectCount, fun = function(x, y) ifelse(y > 0, x/y, 0))

  return(EnemyDensity)
}

#' Apply density-independent loss
#'
#' Calculate density independent predation.
#' Apply a high scale spatial autocorrelated function to the budworm count to simulate
#' density-independent predation.
#' This simulates the effect of summer temperature (huge scale spatial autocorrelation) on ability
#' to escape density-independent predation from thermoregulatory predators.
#'
#' @param PreyCountRaster `RasterLayer`. DESCRIPTION NEEDED.
#'
#' @return egg count DESCRIPTION NEEDED
#'
#' @export
#' @importFrom raster raster
#' @importFrom SpaDES.tools gaussMap
ApplyDensityIndependentLoss <- function(PreyCountRaster) {
  ## Create a map that has the same extent/dimension as budworm count raster and assign values
  ## based on a Gaussian distribution of a scale of equal value of the initial raster row number
  ## to simulate the large scale spatial autocorrelation.

  ## TODO: parameter var need to be determined?
  Predation <- raster(PreyCountRaster)
  Predation <- gaussMap(Predation, scale = dim(Predation)[1], var = 10, method = "RMgauss")

  PreyCount <- PreyCountRaster * Predation

  ## Incorporate probability for prey to successfully complete development based on weather (Regniere 2012)
  ## Random normal map, need to be divided by 10 because variability is between 0-10 and using a
  ## random fraction to remove budworm.
  ## NOTE: THIS IS ARBITRARY

  PhenolLimit <- raster(PreyCountRaster)
  PhenolLimit <- gaussMap(PhenolLimit, scale = dim(PhenolLimit)[1], var = 10, method = "RMgauss")

  EggCount <- PreyCount * (PhenolLimit / 10)
  return(EggCount)
}
