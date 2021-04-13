#' Determine percent defoliation
#'
#' @param ParamA DESCRIPTION NEEDED
#' @param ParamB DESCRIPTION NEEDED
#' @param ParamC DESCRIPTION NEEDED
#' @param Interference DESCRIPTION NEEDED
#' @param PreyDensity DESCRIPTION NEEDED
#' @param PredatorDensity DESCRIPTION NEEDED
#' @param LarvaeConsump DESCRIPTION NEEDED
#' @param lambda DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
DeterminePctDefoliation <- function(ParamA, ParamB, ParamC, Interference, PreyDensity,
                                    PredatorDensity, LarvaeConsump, lambda) {
  rprimeYX <- exp(-1 * (ParamB + Interference) * PreyDensity^ParamA)
  rprimeXY <- exp(-1 * ParamC * PredatorDensity)

  ## mg convert to g * 0.001 factor
  PctDefol <- 100 * LarvaeConsump * 0.001 * PreyDensity * (lambda + (1 - lambda) * rprimeYX * rprimeXY)

  # Coerce larger value then 100% to a 100%
  PctDefol[PctDefol > 100] <- 100

  return(PctDefol)
}

#' Get mortality rate
#'
#' @param species DESCRIPTION NEEDED
#' @param age DESCRIPTION NEEDED
#' @param cumDefol DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
GetMortalityRate <- function(species, age, cumDefol) {
  ## TODO: parameterize function for species code wildCard of sim
  params_MR_AGE_BF <- c(0.1182, -0.003634, 0.00011120, -0.0001563, 0.000002064)
  params_MR_AGE_BS <- c(0.1509, -0.002707, 0.00007719, -0.0001869, 0.000002194)
  params_MR_AGE_RS <- c(0.1509, -0.002707, 0.00007719, -0.0001869, 0.000002194)
  params_MR_AGE_WS <- c(0.1005, -0.002187, 0.00006505, -0.0001375, 0.000001891)
  params_MR_AGE_NONHOST <- c(0.0, 0.0, 0.0, 0.0, 0.0)

  if (cumDefol < 0.35 || is.na(cumDefol)) {
    rate <- 0
    return(rate) # avoid extrapolation issues
  } else {
    ## get param set for specific species
    b <- switch(species,
      "Abies_bal" = params_MR_AGE_BF,
      "Pice_gla" = params_MR_AGE_WS,
      "Pice_mar" = params_MR_AGE_BS,
      "Pice_rub" = params_MR_AGE_RS,
      params_MR_AGE_NONHOST
    )

    ## cap maximal age to avoid extrapolation issues
    if (species == "Abies_bal" && age > 90) {
      age <- 90
    } else if (species == "Pice_gla" && age > 125) {
      age <- 125
    } else if (species == "Pice_mar" && age > 125) {
      age <- 125
    } else if (species == "Pice_rub" && age > 125) {
      age <- 125
    } else {
      age <- age
    }
    ## Periodic (5 year) Mortality Rate = Intercept + AGE + CD*AGE + CD^2 + CD^3
    rate <- b[1] + b[2] * age + b[3] * age * cumDefol + b[4] * cumDefol * cumDefol + b[5] * cumDefol * cumDefol * cumDefol

    ## Trap out of bounds extrapolated predictions.
    if (rate < 0) {
      rate <- 0
    }
    if (rate > 1) {
      rate <- 1
    }

    # Convert from periodic (5 year) to annual mortality rate using the reverse compound interest formula
    # (original fit was based on 5 year mortality rates, so we need to first predict periodic, then convert to annual).
    rate <- 1 - (1 - rate)^(1 / 5.0)

    return(rate)
  }
}
