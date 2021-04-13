#' ConvertToEnemyDensity
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
