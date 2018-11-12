SeparateKRpkmNew2 <- function(x, n, q, err = 1e-10) {
  k <- 1
  q <- max(q, min(x))
  c <- sum(x < q)
  x <- x[which(x >= q)]
  if (length(x) <= k) {
    warning(sprintf("The length of x is %i. Sorry, too little conditions\n", length(x)))
    return(cbind(0, 0, 0))
  }
  mean <- c()
  for (i in 1:k) {
    mean <- c(mean, sort(x)[floor(i * length(x)/(k + 1))])
  }
  if(k>1)
  {
    mean[1] <- min(x) - 1
    mean[length(mean)] <- max(x) + 1
  }
  p <- rep(1/k, k)
  sd <- rep(sqrt(var(x)), k)
  pdf.x.portion <- matrix(0, length(x), k)
  for (i in 1:n) {
    p0 <- p
    mean0 <- mean
    sd0 <- sd
    pdf.x.all <- t(p0 * vapply(x, function(x) dnorm(x, mean0, sd0), rep(0, k)))
    pdf.x.portion <- pdf.x.all/rowSums(t(pdf.x.all))
    cdf.q <- pnorm(q, mean0, sd0)
    cdf.q.all <- p0 * cdf.q
    cdf.q.portion <- cdf.q.all/sum(cdf.q.all)
    cdf.q.portion.c <- cdf.q.portion * c
    denom <- colSums(t(pdf.x.portion)) + cdf.q.portion.c
    p <- denom/(nrow(t(pdf.x.portion)) + c)
    im <- dnorm(q, mean0, sd0)/cdf.q * sd0
    im[is.na(im)] <- 0
    mean <- colSums(crossprod(x, t(pdf.x.portion)) + (mean0 - sd0 * im) * cdf.q.portion.c)/denom
    sd <- sqrt((colSums((x - matrix(mean0, ncol = length(mean0), nrow = length(x), byrow = TRUE)) ^ 2
                        * t(pdf.x.portion)) + sd0 ^ 2 * (1 - (q - mean0) / sd0 * im) * cdf.q.portion.c) / denom)
    if (!is.na(match(NaN, sd))) {
      break
    }
    if ((mean(abs(p - p0)) <= err) &&
        (mean(abs(mean - mean0)) <= err) &&
        (mean(abs(sd - sd0)) <= err)) {
      break
    }
  }
  return(cbind(p, mean, sd))
}

#' Calcuate
#'
#' @param x data, example: x<-runif(100,0,1)
#' @param n rounds
#' @param q cutoff
#' @param k k=1..5
#' @param err
#'
#' @return a matrix contains pi, mean and sd
#' @export
#'
#' @examples
SeparateKRpkmNew <- function(x, n, q, k, err = 1e-10) {
  if (k == 1) return(SeparateKRpkmNew2(x, n, q, err = 1e-10))
  q <- max(q, min(x))
  c <- sum(x < q)
  x <- x[which(x >= q)]
  if (length(x) <= k) {
    warning(sprintf("The length of x is %i. Sorry, too little conditions\n", length(x)))
    return(cbind(0, 0, 0))
  }
  mean <- c()
  for (i in 1:k) {
    mean <- c(mean, sort(x)[floor(i * length(x) / (k + 1))])
  }
  mean[1] <- min(x) - 1  # What is those two lines for?
  mean[length(mean)] <- max(x) + 1  # Without them the result of mean[1] is slightly different.
  p <- rep(1 / k, k)
  sd <- rep(sqrt(var(x)), k)
  pdf.x.portion <- matrix(0, length(x), k)

  for (i in 1:n) {
    p0 <- p
    mean0 <- mean
    sd0 <- sd

    pdf.x.all <- t(p0 * vapply(x, function(x) dnorm(x, mean0, sd0), rep(0, k)))
    pdf.x.portion <- pdf.x.all / rowSums(pdf.x.all)
    cdf.q <- pnorm(q, mean0, sd0)
    cdf.q.all <- p0 * cdf.q
    cdf.q.portion <- cdf.q.all / sum(cdf.q.all)
    cdf.q.portion.c <- cdf.q.portion * c
    denom <- colSums(pdf.x.portion) + cdf.q.portion.c
    p <- denom / (nrow(pdf.x.portion) + c)
    im <- dnorm(q, mean0, sd0) / cdf.q * sd0
    im[is.na(im)] <- 0
    mean <- colSums(crossprod(x, pdf.x.portion) + (mean0 - sd0 * im) * cdf.q.portion.c) / denom
    sd <- sqrt((colSums((x - matrix(mean0, ncol = length(mean0), nrow = length(x),
      byrow = TRUE)) ^ 2 * pdf.x.portion) + sd0 ^ 2 * (1 - (q - mean0) / sd0 * im) *
      cdf.q.portion.c) / denom)
    if (!is.na(match(NaN, sd))) {
      break
    }
    if ((mean(abs(p - p0)) <= err) && (mean(abs(mean - mean0)) <= err) &&
      (mean(abs(sd - sd0)) <= err)) {
      break
    }
  }
  return(cbind(p, mean, sd))
}


# x: data, x<-runif(100,0,1), n: 300, (cutoff for stop, diff of total up ABS<1e-6, q=0 for testing data, k=1...5)

SeparateKRpkmNewp <- function(x, n, q, k, err = 1e-10) {
  q <- max(q, min(x))
  c <- sum(x < q)
  x <- x[which(x >= q)]
  mean <- c()
  for (i in 1:k) {
    mean <- c(mean, sort(x)[floor(i * length(x) / (k + 1))])
  }
  mean[1] <- min(x) - 1  # What is those two lines for?
  mean[length(mean)] <- max(x) + 1  # Without them the result of mean[1] is slightly different.
  p <- rep(1 / k, k)
  sd <- rep(sqrt(var(x)), k)
  t <- matrix(0, length(x), k)

  for (i in 1:n) {
    p0 <- p
    mean0 <- mean
    sd0 <- sd

    for (row in 1:nrow(t)) {
      all <- p0 * dnorm(x[row], mean0, sd0)
      t[row, ] <- all / sum(all)
    }

    pZil <- Pi_Zj_Zcut_new(q, mean0, sd0, p0) ################################################################
    denom <- (colSums(t) + pZil * c)
    all <- denom / (nrow(t) + c)
    p <- all / sum(all)
    im <- InverseMillsRatio(q, mean0, sd0)
    mean <- colSums(crossprod(x, t) + (mean0 - sd0 * im) * pZil * c) / denom

    a <- rep(0, length(sd))

    for (col in 1:length(sd)) {
      a[col] <- sum((x - mean0[col])^2 * t[, col])
    }

    sd <- sqrt((a + (sd0)^2 * (1 - (q - mean0) / sd0 * im) * pZil * c) / denom)
    if ((mean(abs(p - p0)) <= err) && (mean(abs(mean - mean0)) <= err) && (mean(abs(sd - sd0)) <= err)) {
      break
    }
  }
  return (cbind(p, mean, sd))
}

LogSeparateKRpkmNew <- function(x, n, q, k, err = 1e-10) {
  return (SeparateKRpkmNew(log(x), n, log(q), k, err))
}
