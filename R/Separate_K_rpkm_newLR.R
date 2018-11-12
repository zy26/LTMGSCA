#' Calcuate the LTMG_2LR for some genes
#'
#' @param x data, a List of NumericVector
#' @param n rounds
#' @param q cutoff of the elements in x
#' @param r maximum value of the standard diversion
#' @param s minimum value of the standard diversion
#' @param k number of peaks, should be 2
#' @param err the upper bound on the absolute error
#'
#' @return a matrix contains pi, mean and sd
#' @export
#'
#' @examples
SeparateKRpkmNewLRPlus <- function(x, n, q, r, s = 0.05, k = 2, err = 1e-10, M = Inf, m = -Inf) {
  c <- sapply(x, function(x) sum(x < q), simplify = "array")

  x_r0 <- lapply(x, function(x) x[which(x < q)])
  x_r <- lapply(x, function(x) x[which(x >= q)])
  x_r0.length <- sapply(x_r0, length)
  x_r.length <- sapply(x_r, length)
  x_r.non.zero <- x_r.length > 0
  x_r.non.allpos<- x_r0.length > 0
  x_r.input<-(x_r.non.allpos*x_r.non.zero)>0

  x_r <- x_r[x_r.input]

  if (sum(x_r.non.zero*x_r.non.allpos) == 0) {
    warning("Completely all 0/all pos conditions\n")
    results_c<-list()
    for(i in 1:length(x))
    {
      ccc<-matrix(0,2,3)
      colnames(ccc)<-c("p","mean","sd")
      if(x_r.non.zero[i]==0)
      {
          ccc[1,1]<-1
          ccc[2,1]<-0
          ccc[1,2]<-m
          ccc[2,2]<-M
          ccc[1,3]<-s
          ccc[2,3]<-s
      }
      if(x_r.non.allpos[i]==0)
      {
        ccc[1,1]<-0
        ccc[2,1]<-1
        ccc[1,2]<-m
        ccc[2,2]<-mean(x[[i]])
        ccc[1,3]<-s
        ccc[2,3]<-sd(x[[i]])
      }
      results_c[[i]]<-ccc
    }
    return(results_c)
  }

  if(is.na(max(x_r.length[x_r.input])<5))
  {
    browser()
  }
  if(max(x_r.length[x_r.input])<5)
  {
    warning("Too little non-zero part, forced ZIG\n")
    results_c<-list()
    for(i in 1:length(x))
    {
      ccc<-matrix(0,2,3)
      colnames(ccc)<-c("p","mean","sd")
      if(x_r.non.zero[i]==0)
      {
        ccc[1,1]<-1
        ccc[2,1]<-0
        ccc[1,2]<-m
        ccc[2,2]<-M
        ccc[1,3]<-s
        ccc[2,3]<-s
      }
      if(x_r.non.allpos[i]==0)
      {
        ccc[1,1]<-0
        ccc[2,1]<-1
        ccc[1,2]<-m
        ccc[2,2]<-mean(x[[i]])
        ccc[1,3]<-s
        ccc[2,3]<-sd(x[[i]])
      }
      if((x_r.non.allpos[i]!=0)&(x_r.non.zero[i]!=0))
      {
        ccc[1,1]<-sum(x[[i]]<q)/length(x[[i]])
        ccc[2,1]<-sum(x[[i]]>=q)/length(x[[i]])
        ccc[1,2]<-m
        ccc[2,2]<-mean(x[[i]][which(x[[i]]>=q)])
        ccc[1,3]<-s
        ccc[2,3]<-sd(x[[i]][which(x[[i]]>=q)])
        if(is.na(ccc[2,3]))
        {
          ccc[2,3]<-s
        }
      }
      results_c[[i]]<-ccc
    }
    return(results_c)
  }

  c<- c[x_r.input]
  c_sum <- sum(c)

  ncol <- length(x_r)

  x_all <- c(x_r, recursive = TRUE)

  p <- matrix(1 / k, nrow = k, ncol = ncol)

  mean <- matrix(nrow = k, ncol = ncol)
  for (col in 1:ncol) {
    ni <- length(x_r[[col]])
    for (row in 1:k) {
      tg_ic <- floor(row * ni / (k + 1))
      cc <- sort(x_r[[col]])[tg_ic]
      if (tg_ic < 1) {
        cc <- sort(x_r[[col]])[1] - 0.5
      }
      if (tg_ic > ni) {
        cc <- sort(x_r[[col]])[ni] + 0.5
      }
      mean[row, col] <- cc
    }
  }

  sd <- matrix(sqrt(vapply(x_r, var, 0)), nrow = k, ncol = ncol, byrow = TRUE)

  if(anyNA(sd)) {
    sd[is.na(sd)] <- 1
  }

  sd[which(sd<s)]<-s

  p0 <- p
  mean0 <- mean
  sd0 <- sd

  t <- lapply(x_r, function(x) matrix(nrow = length(x), ncol = k))
  for (i in 1:n) {
    ccc <- matrix(nrow = k, ncol = ncol)
    wad <- rep(0, ncol + 1)
    mean_all <- rep(0, ncol + 1)
    sd_all <- rep(0, ncol + 1)
    sd_all_1_sum <- 0
    for (col in 1:ncol) {
      t_u <- t(p[, col] * vapply(x_r[[col]], function(x) dnorm(x, mean[, col], sd[, col]), rep(0, k)))
      t[[col]] <- t_u / rowSums(t_u)
      pZil2 <- Pi_Zj_Zcut_new(q, mean[, col], sd[, col], p[, col])
      denom2 <- (colSums(t[[col]]) + pZil2 * c[[col]])
      im <- InverseMillsRatio(q, mean[, col], sd[, col])
      mean0[, col] <- colSums(crossprod(x_r[[col]], t[[col]]) + (mean[, col] - sd[, col] * im) * pZil2 * c[[col]]) / denom2
      if(denom2[1]==0)
      {
        mean0[1, col]<-mean[1, col]
      }
      if(denom2[2]==0)
      {
        mean0[2, col]<-mean[2, col]
      }
      if (anyNA(denom2)) {
        warning("denom2 conttains NA\n")
        print(x_r)
        print(col)
        print(dnorm(x_r[[col]][1], mean[1, col], sd[1, col]))
        print(dnorm(x_r[[col]][2], mean[2, col], sd[2, col]))
        browser()
      }
      if (anyNA(mean0)) {
        warning("mean0 conttains NA\n")
      }
      if (mean0[1, col] > q) {#denom2[1] == 0 ||
        mean0[1, col] <- q
      }
      if (mean0[2, col] < q) {#denom2[2] == 0 ||
        mean0[2, col] <- q
      }
      sd0[, col] <- sqrt((colSums((x_r[[col]] - matrix(mean[, col], ncol = length(mean[, col]), nrow = length(x_r[[col]]),
        byrow = TRUE)) ^ 2 * t[[col]]) + (sd[, col]) ^ 2 * (1 - (q - mean[, col]) / sd[, col] * im) * pZil2 * c[[col]]) / denom2)
      if (denom2[1] == 0 || sd0[1, col] > r) {
        sd0[1, col] <- r
      }
      if (denom2[2] == 0 || sd0[2, col] > r) {
        sd0[2, col] <- r
      }

      wl2 <- denom2 / (nrow(t[[col]]) + c[[col]])

      # Reorder by mean0
      tg_R <- order(mean0[, col])
      ccc[, col] <- wl2[tg_R] * (c[[col]] + length(x_r[[col]])) / (sum(c) + sum(lengths(x_r)))
      mean0[, col] <- mean0[, col][tg_R]
      sd0[, col] <- sd0[, col][tg_R]

      wad[[col + 1]] <- ccc[-1, col]
      wad[1] <- wad[1] + ccc[1, col]
      mean_all[[col + 1]] <- mean0[, col][-1]
      mean_all[1] <- mean_all[1] + mean0[1, col] * c[[col]] / c_sum
      sd_all[[col + 1]] <- sd0[, col][-1]
      sd_all_1_sum <- sd_all_1_sum + sd0[1, col] ^ 2 * c[[col]]
    }
    sd_all[1] <- sqrt(sd_all_1_sum / c_sum)

    t0 <- matrix(nrow = sum(vapply(x_r, length, 0)), ncol = k + ncol - 1)
    for (row in 1:nrow(t0)) {
      t0_u <- wad * dnorm(x_all[row], mean_all, sd_all)
      t0[row, ] <- t0_u / sum(t0_u)
    }

    pZil0 <- Pi_Zj_Zcut_new(q, mean_all, sd_all, wad)
    denom0 <- colSums(t0) + pZil0 * c_sum

    im1 <- InverseMillsRatio(q, mean_all[1], sd_all[1])
    mean0[1, ] <- (sum(x_all * t0[, 1]) + (mean_all[1] - sd_all[1] * im1) * pZil0[1] * c_sum) / denom0[1]
    sd0[1, ] <- sqrt((sum((x_all - mean_all[1]) ^ 2 * t0[, 1]) + sd_all[1] ^ 2 * (1 - (q - mean_all[1]) / sd_all[1] * im1) *
      pZil0[1] * c_sum) / denom0[1])

    for (col in 1:ncol) {
      t_u <- t(p[, col] * vapply(x_r[[col]], function(x) dnorm(x, mean[, col], sd[, col]), c(0, 0)))
      t[[col]] <- t_u / rowSums(t_u)
      pZil <- Pi_Zj_Zcut_new(q, mean[, col], sd[, col], p[, col])
      p_u <- (colSums(t[[col]]) + pZil * c[[col]]) / (nrow(t[[col]]) + c[[col]])
      p0[, col] <- p_u / sum(p_u)
    }

    for (col in 1:ncol) {
      sd0[1, col] <- min(sd0[1, col], r)
      sd0[2, col] <- max(min(sd0[2, col], r), s)
      mean0[1, col] <- min(mean0[1, col], q)
      mean0[2, col] <- max(mean0[2, col], q)
    }

    #print(i)
    #print(p0)
    #print(mean0)
    #print(sd0)

    if (anyNA(sd0)) {
      warning("Found at least one NA in sd")
      break
    }

    if ((mean(abs(p - p0)) <= err) && (mean(abs(mean - mean0)) <= err) && (mean(abs(sd - sd0)) <= err)) {
      break
    }

    p <- p0
    mean <- mean0
    sd <- sd0
  }

  print(i)

  ret <- vector("list", length(x))

  mean_peak1 = mean[1, 1]
  sd_peak1 = sd[1, 1]

  col <- 1
  for (index in which((x_r.non.zero*x_r.non.allpos)==1)){
    ret[[index]] <- cbind(p = p[, col], mean = mean[, col], sd = sd[, col])
    col <- col + 1
  }
  for (index in 1:length(x)){
    if((x_r.non.zero[index]*x_r.non.allpos[index])==0)
    {
    ccc<-matrix(0,2,3)
    colnames(ccc)<-c("p","mean","sd")
    if(x_r.non.zero[index]==0)
    {
      ccc[1,1]<-1
      ccc[2,1]<-0
      ccc[1,2]<-mean_peak1
      ccc[2,2]<-M
      ccc[1,3]<-sd_peak1
      ccc[2,3]<-s
    }
    if(x_r.non.allpos[index]==0)
    {
      ccc[1,1]<-0
      ccc[2,1]<-1
      ccc[1,2]<-mean_peak1
      ccc[2,2]<-mean(x[[index]])
      ccc[1,3]<-sd_peak1
      ccc[2,3]<-sd(x[[index]])
    }
    ret[[index]] <- ccc
    }
  }
  return(list(ret, i))
}

SeparateKRpkmNewLR <- function(x, n, q, r, s = 0.05, k = 2, err = 1e-10, M = Inf, m = -Inf) {
 return (SeparateKRpkmNewLRPlus(x, n, q, r, s, k, err, M, m)[[1]])
}

LogSeparateKRpkmNewLR <- function(x, n, q, r, k = 2) {
  return(SeparateKRpkmNewLR(log(x), n, log(q), r, k))
}
