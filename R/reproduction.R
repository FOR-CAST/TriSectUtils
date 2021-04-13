#' Calculate predator reproductive rate
#'
#' @param ParamB DESCRIPTION NEEEDED
#' @param ParamC DESCRIPTION NEEEDED
#' @param BackgrPred DESCRIPTION NEEEDED
#' @param PreyDensity DESCRIPTION NEEEDED
#' @param PredatorDensity DESCRIPTION NEEEDED
#' @param maxReproEnemy DESCRIPTION NEEEDED
#'
#' @return DESCRIPTION NEEEDED
#'
#' @export
CalculatePredatorReproductiveRate <- function(ParamB, ParamC, BackgrPred, PreyDensity,
                                              PredatorDensity, maxReproEnemy) {
  rYX <- 1 - exp(-1 * ParamB * PreyDensity - BackgrPred)
  rXY <- exp(-1 * ParamC * PredatorDensity)
  rTotal <- maxReproEnemy * rYX * rXY
  return(rTotal)
}

#' Calculate prey reproductive rate
#'
#' @param ParamA DESCRIPTION NEEEDED
#' @param ParamB DESCRIPTION NEEEDED
#' @param ParamC DESCRIPTION NEEEDED
#' @param Interference DESCRIPTION NEEEDED
#' @param PreyDensity DESCRIPTION NEEEDED
#' @param PredatorDensity DESCRIPTION NEEEDED
#' @param MatingEff DESCRIPTION NEEEDED
#' @param maxReproPrey DESCRIPTION NEEEDED
#' @param Fecundity DESCRIPTION NEEEDED
#'
#' @return DESCRIPTION NEEEDED
#'
#' @export
CalculatePreyReproductiveRate <- function(ParamA, ParamB, ParamC, Interference, PreyDensity,
                                          PredatorDensity, MatingEff, maxReproPrey, Fecundity) {
  rprimeYX <- exp(-1 * (ParamB + Interference) * PreyDensity^ParamA)
  rprimeXY <- exp(-1 * ParamC * PredatorDensity)
  rprimeTotal <- Fecundity * maxReproPrey * rprimeYX * rprimeXY * MatingEff
  return(rprimeTotal)
}

#' Add starvation to prey reproductive rate
#'
#' @param PctDefol DESCRIPTION NEEEDED
#' @param Constant DESCRIPTION NEEEDED
#' @param PreyReproRate DESCRIPTION NEEEDED
#'
#' @return DESCRIPTION NEEEDED
#'
#' @export
AddStarvationToPreyReproductiveRate <- function(PctDefol, Constant, PreyReproRate) {
  ## NOTE: Has to be divided by 100 because contrary to LANDIS code, pctdefol is the percentdefol
  ##       and not the proportion defol!
  rprimeZY <- Constant * (PctDefol / 100) + 1
  rprime2Total <- rprimeZY * PreyReproRate
  return(rprime2Total)
}
