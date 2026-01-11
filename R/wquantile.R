#' Function from Chirs Sherlock notes on Bootstrap filter for compute
#' weighted quantiles
wquantile <- function(xs, ws, prob) {
  ws <- ws[order(xs)]
  xs <- xs[order(xs)]
  sws <- sum(ws)
  cws <- cumsum(ws)
  np <- length(prob)
  qs <- rep(0,np)
  for (i in 1:np) {
    qs[i] <- xs[which.max(cws / sws >= prob[i])]
  }
  qs
}

