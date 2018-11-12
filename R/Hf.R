InverseMillsRatio <- function(q, mean, sd) {
  x <- (q - mean) / sd
  pdf <- dnorm(x, log = TRUE)
  cdf <- pnorm(x, log = TRUE)
  d <- exp(pdf - cdf)
  d[is.na(d)] <- 0
  return(d)
}
