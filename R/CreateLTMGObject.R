CreateLTMGObject <- function(x, cond) {
  assertthat::assert_that(ncol(x) == ncol(cond))
  return(new("LTMG", Data = x, Normalized.Data = x, Conds = cond, Conds.Label = 2^seq(0, nrow(cond) - 1, by = 1) %*% cond))
}
