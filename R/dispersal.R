globalVariables(c(
  ".N", ".SD", "V1", "V2", "V3", "V4"
))

#' CalcWinterMortalityAndDispersal
#'
#' Calculate overwintering loss and do a focal mean to simulate L2 dispersal.
#' Focal mean is done with equal weight of the 8 neighbor cells by default.
#'
#' @param WindowDistance DESCRIPTION NEEDED
#' @param Count DESCRIPTION NEEDED
#' @param SurvivingPercent DESCRIPTION NEEDED
#' @param LDisperse DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @author Louis-Etienne Robert
#' @export
#' @importFrom raster focal
CalcWinterMortalityAndDispersal <- function(WindowDistance, Count, SurvivingPercent,LDisperse) {
  #Debug
  #a<<-Count
  #b<<-SurvivingPercent

  ## Overwintering loss
  CountWinterMort <- Count * SurvivingPercent

  if (isTRUE(LDisperse)) {
    ## establish window matrix of cellvalues=1 where the number of rows/column = 2x+1 and
    ## total amount of cells is (2x+1)^2, x=DistanceofFocalWindow

    FocalWindow <- matrix(rep(1, (2*WindowDistance + 1)^2), nrow = 2*WindowDistance + 1, ncol = 2*WindowDistance + 1)

    ## Focal mean for Dispersal of L2. pad value for length of window matrix to account for edges.
    ## Then take the floor value since fraction of population is impossible
    CountWinterMort <- floor(focal(CountWinterMort, FocalWindow, meanNZ, pad = TRUE, padValue = 0))
  }

  ## TODO: handling edges in a toroidal geometry
  ##   Check wrap function in Spades.tools
  ## TODO: UNTESTED: Added an if trigger for local dispersal: CountWinterMort updating might not
  ##   work as expected with the trigger in place
  return(CountWinterMort)
}

#' \code{GetLDDHabitat}
#'
#' @param defol DESCRIPTION NEEDED
#' @param MinLDD DESCRIPTION NEEDED
#' @param HalfLDD DESCRIPTION NEEDED
#' @param MaxLDD DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
GetLDDHabitat <- function(defol, MinLDD, HalfLDD, MaxLDD) {
  # catch missing value #TODO verify consequences over whole model
  if (is.na(defol)) {
    LDDHabitat <- 0
    return(LDDHabitat)
  }
  if (defol <= MinLDD) {
    LDDHabitat <- 0
  } else if (defol > MinLDD && defol <= HalfLDD) {
    if (HalfLDD - MinLDD > 0) {
      m1 <- 0.5 / (HalfLDD - MinLDD)
    } else {
      m1 <- 0
    }
    b1 <- 0.5 - (m1 * HalfLDD)
    LDDHabitat <- m1 * defol + b1
  } else if (defol > HalfLDD && defol <= MaxLDD) {
    if (MaxLDD - HalfLDD > 0) {
      m1 <- 0.5 / (MaxLDD - HalfLDD)
    } else {
      m1 <- 0
    }
    b1 <- 1.0 - (m1 * MaxLDD)
    LDDHabitat <- m1 * defol + b1
  } else {
    LDDHabitat <- 1
  }
  return(LDDHabitat)
}

#' \code{GetLDDFlight}
#'
#' Default values of Max,Half,Min bind the data between 0 and 1 and use the defol value as the lddflight value
#'
#' @param rZY TODO
#' @param MaxLDDProp TODO
#' @param PosRel TODO
#'
#' @return TODO
#' @export
GetLDDFlight <- function(rZY, MaxLDDProp, PosRel) {
  # catch missing values #TODO verify consequences over whole model
  if (is.na(rZY)) {
    LDDFlight <- 0
  } else {
    slope <- (MaxLDDProp - (1 - MaxLDDProp)) / (1 - 0.46)
    intercept <- MaxLDDProp - slope
    if (PosRel) {
      if (rZY < 0.46) {
        LDDFlight <- 0
      } else {
        LDDFlight <- slope * rZY + intercept
      }
    } else {
      slope2 <- -1.0 * slope
      intercept2 <- -1.0 * intercept + 1.0

      if (rZY < 0.46) {
        LDDFlight <- 1.0
      } else {
        LDDFlight <- slope2 * rZY + intercept2
      }
    }
  }

  return(LDDFlight)
}

#' Calculate long-distance dispersal ratio
#'
#' @param PctDefol TODO
#' @param MinLDD TODO
#' @param HalfLDD TODO
#' @param MaxLDD TODO
#' @param MaxLDDProp TODO
#' @param PositiveRelation TODO
#' @param Constant TODO
#'
#' @export
#' @importFrom raster as.matrix extent raster
CalcLDDRatio <- function(PctDefol, MinLDD, HalfLDD, MaxLDD, MaxLDDProp, PositiveRelation, Constant) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test

  # need to create output in metadata
  ### XXXXXXXXXX PctDefol IS INCORRECT Check pseudo code when changing biomass input
  LDDHabitat <- apply(raster::as.matrix((PctDefol / 100)), MARGIN = c(1, 2), FUN = GetLDDHabitat,
                      MinLDD = MinLDD, HalfLDD = HalfLDD, MaxLDD = MaxLDD)

  # LDD flight computation
  # TODO refer to recruit and defoliate for arguments
  rprimeZY <- Constant * (PctDefol / 100) + 1
  LDDFlight <- apply(raster::as.matrix(rprimeZY), MARGIN = c(1, 2), FUN = GetLDDFlight,
                     MaxLDDProp = MaxLDDProp, PosRel = PositiveRelation)

  # Calculate the LDDratio
  LDDRatio <- raster::raster(LDDHabitat * LDDFlight)
  extent(LDDRatio) <- extent(rprimeZY)

  return(LDDRatio)
}

#' \code{DispByDistance}
#'
#' @param r DESCRIPTION NEEDED
#' @param maxDistThreshold DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom data.table data.table uniqueN
#' @importFrom fields fields.rdist.near
#' @importFrom raster rasterFromXYZ rasterToPoints
DispByDistance <- function(r, maxDistThreshold = 10) { # r@nrows * 1000
  ## 1. transform to points
  r2 <- rasterToPoints(r) # can use function argument to subset raster example: dont use 0

  ## 2. calculate distance using a max threshold distance (delta)
  distance <- fields.rdist.near(r2[, 1:2], r2[, 1:2], delta = maxDistThreshold, max.points = nrow(r2) * 1000)

  ## 3. convert distance to probabilities

  # calculate the probability of dispersal to the "to" cell as function of distance (normalized inversed probability)
  dist.invers <- cbind(distance[[1]], 1 / distance[[2]])
  dist.invers[dist.invers[, 3] == Inf, 3] <- 0 # correction for distance = 0

  a <- stats::aggregate(dist.invers[, 3], by = list(dist.invers[, 1]), FUN = sum) # get the sum of distance by "from" cell
  # introduce in the data frame
  b <- merge(dist.invers, a, by = 1)
  c <- merge(b, r2[, 3], by.x = 1, by.y = 0)
  t.data <- cbind(distance[[1]], c[, 3] / c[, 4], c[, 5])

  ## 4. convert to data.table
  t.data <- data.table(t.data)

  ## take group of data (by V1) and sample with weight V3 and size V4 (unique value) return as data.table
  t.data2 <- t.data[, lapply(.SD, sample, size = uniqueN(V4), replace = TRUE, prob = V3), by = V1]

  # count and sort by cell index
  t.data3 <- t.data2[, .N, keyby = V2]

  ## 5. assign to matrix of point and back to raster
  # r.check <- cbind(r2, newValue) # this is to test the output for debug
  r3 <- merge(r2, t.data3, by.x = 0, by.y = 1, all.x = TRUE, sort = TRUE)
  r3[is.na(r3[, 5]), 5] <- 0
  r4 <- r3[, c(2, 3, 5)]
  LDDDispersed <- rasterFromXYZ(r4)

  return(LDDDispersed)
}

#' Short-distance dispersal
#'
#' @param r TODO
#'
#' @return TODO
#'
#' @export
#' @importFrom raster focal
DispSDD <- function(r) {
  ## establish window matrix of cellvalues=1 where the number of rows/column = 2x+1 and
  ## total amount of cells is (2x+1)^2, x=DistanceofFocalWindow

  FocalWindow <- matrix(rep(1, 9), 3, 3)

  ## Focal mean. pad value for length of window matrix to account for edges.
  ## Then take the floor value since fraction of population is impossible
  SDDDispersed <- floor(focal(r, FocalWindow, sum, pad = TRUE, padValue = 0))
  # SDDDispersed <- floor(focal(r, FocalWindow, meanNZ, pad = TRUE, padValue = 0))

  return(SDDDispersed)
}

#' LocalDispersalEnemy
#'
#' Calculate overwintering loss and do a focal mean to simulate local dispersal.
#' Focal mean is done with equal weight of the 8 neighbour cells.
#'
#' @param EnemyCount `RasterLayer`. DESCRIPTION NEEDED
#' @param InsectCount `RasterLayer`. DESCRIPTION NEEDED
#' @param WindowDistanceEnemy DESCRIPTION NEEDED
#' @param Proportion DESCRIPTION NEEDED
#' @param PreyBiasProportion DESCRIPTION NEEDED
#'
#' @return DESCRIPTION NEEDED
#'
#' @export
#' @importFrom raster focal maxValue
LocalDispersalEnemy <- function(EnemyCount, InsectCount, WindowDistanceEnemy, Proportion, PreyBiasProportion) {
  ## establish window matrix of cellvalues=1 where the number of rows/column = 2x+1 and
  ## total amount of cells is (2x+1)^2, x=DistanceofFocalWindow
  FocalWindow <- matrix(rep(1, (2 * WindowDistanceEnemy + 1)^2), nrow = 2 * WindowDistanceEnemy + 1,
                        ncol = 2 * WindowDistanceEnemy + 1)

  ## Disperse Enemy with a spatial filter. A portion of the pop. stays on site depending on user-defined setting
  maxVal <- ifelse(maxValue(InsectCount) == 0, 1, maxValue(InsectCount)) # Avoid division by 0 if insect count is all 0
  TotalPropDisperse <- Proportion + (PreyBiasProportion * (1 - (InsectCount / maxVal)))

  ## Coerce proportion value to 1 if bigger then 1
  TotalPropDisperse[TotalPropDisperse > 1] <- 1
  ## Determine dispersing population
  PoptoDisperse <- EnemyCount * TotalPropDisperse
  PopStaying <- floor(EnemyCount - PoptoDisperse)

  ## Focal mean. pad value for length of window matrix to account for edges.
  ## Then take the ceiling value since fraction of population is impossible and want to avoid complete extinction
  Dispersed <- floor(focal(PoptoDisperse, FocalWindow, meanNZ, pad = TRUE, padValue = 0))

  ## Add back the population count that stayed on site
  EnemyCount <- PopStaying + Dispersed

  return(EnemyCount)
}
