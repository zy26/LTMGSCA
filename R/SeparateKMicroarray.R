SeparateKMicroarray <- function(x, n, k, err = 1e-10) {
  mean <- c()
  for (i in 1:k) {
    mean <- c(mean, sort(x)[floor(i * length(x)/(k + 1))])
  }
  mean[1] <- min(x) - 1  #######################################
  mean[length(mean)] <- max(x) + 1  #######################################
  p <- rep(1/k, k)
  sd <- rep(sqrt(var(x)), k)
  t <- matrix(0, length(x), k)

  for (i in 1:n) {
    p0 <- p
    mean0 <- mean
    sd0 <- sd

    for (row in 1:nrow(t)) {
      all <- p * dnorm(x[row], mean, sd)
      t[row, ] <- all/sum(all)
    }
    p <- colMeans(t)
    for (col in 1:length(sd)) {
      sd[col] <- sqrt(sum((x - mean[col])^2 * t[, col])/(length(x) * p[col]))
    }

    mean <- colSums(crossprod(x, t))/(length(x) * p)

    if ((mean(abs(p - p0)) <= err) && (mean(abs(mean - mean0)) <= err) && (mean(abs(sd - sd0)) <= err)) {
      break
    }
  }
  return(cbind(p, mean, sd))
}
