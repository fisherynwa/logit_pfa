pfa = function (Z, 
                Sigma, t, Kmax,
                reg = "L1", e = 0.05,
                gamma, K, inspection = FALSE,
                plot = "-log") {
  
#############################################################################
# There is a slight difference in the PFA functions between simulation study and case study.
# Namely, we have erased one line of code that computes the vector of Z-values -
# i.e. Z <- Z/SD is deleted, where Z are the marginal estimates.
#Since we applied an empirical correction method, 
#and this vector was computed prior to calling the pfa
#############################################################################
  Z <- as.vector(Z)
  Sigma <- as.matrix(Sigma)
  p <- length(Z)
  
#############################################################################
#Estimating correlation matrix based on the empirical covariance matrix ######
#############################################################################
    
  SD <- sqrt(diag(Sigma))
  Sigma <- diag(1/SD) %*% Sigma %*% diag(1/SD)
  
############################################################################# 
#### Running PCA on the correlation matrix based on the test statistics
#############################################################################
  pca <- svd(Sigma, nu = 0, nv = Kmax)
  
  lambda <- pca$d
  eigvec <- pca$v
############################################################################# 
## Automatic way for determing k, for more details:
## The eigenvalue ratio (ER) estimator in Ahn & Horenstein (2013)
############################################################################# 
  if (missing(K)) {
    K = 1
    while (K < Kmax & sqrt(sum(lambda[(K + 1):length(lambda)]^2)) >= 
           e * sum(lambda)) K = K + 1
  }
  
  sqrt_lambda <- as.matrix(sqrt(lambda[1:K]))
  b <- as.matrix(eigvec[, 1:K])
  
  for (i in 1:K) {
    b[, i] <- b[, i] * sqrt_lambda[i]
  }
  
############################################################################# 
#### The \hat{w} estimator for w, then these k values of w are plugged in (15).
# If 'L1' then we estimate w as given in (17); if 'L2' then we get w as in (16).
#############################################################################
  if (reg == "L1") {
    W.hat <- rq(Z ~ b - 1, 0.5)$coef
  }
  else if (reg == "L2") {
    o = order(abs(Z))
    Zperm = Z[o]
    Lperm = as.matrix(b[o, ])
    Z.reduce = Zperm[1:(p * 0.95)]
    L.reduce = as.matrix(Lperm[1:(p * 0.95), ])
    W.hat = lsfit(x = L.reduce, y = Z.reduce, intercept = F)$coef
  }
  rs <- rowSums(b^2)
  inv_a <- sqrt(((1 - rs) + abs(1 - rs))/2)
  bW.est <- b %*% (W.hat)
  P <- 2 * (1 - pnorm(abs(Z)))
  sort <- sort(P, index.return = TRUE)
  index <- sort$ix
  P <- sort$x
  if (missing(gamma)) 
    gamma <- as.numeric(quantile(P, probs = 0.4))
  p0.est <- min(p, sum(P > gamma)/(1 - gamma))
  t.default <- TRUE

############################################################################# 
#### Selection of the plausible cutoffs of (t). Fan et al. have implemented a default way 
#############################################################################
  
  if (!missing(t)) {
    if (t[1] == "pval") {
      t <- P
      t.default <- FALSE
    }
    
    if (is.numeric(t)) 
      t.default = (sum(t >= 0) + sum(t <= 1) < 2 * length(t))
  }
  if (t.default) {
    logt.l <- max(min(log(P)), log(0.00000000000001))
    logt.u <- max(log(P))
    grid <- (logt.u - logt.l) * seq(from = 0.01, to = 1, 
                                    by = 0.025) * 0.5 + 0.85 * logt.l + 0.15 * logt.u
    t <- exp(grid)
  }
############################################################################# 
# Evaluation of FDP(t) over the set of t values.
#############################################################################  
  print(paste0('......just one more moment')) 
  FDPt <- Vt <- Rt <- rep(0, length(t))
  for (l in 1:length(t)) {
    P1 <- 2 * (1 - pnorm(abs(Z)))
    Rt[l] <- sum(P1 <= t[l]) ## Compute R(t), as given in Step 2 in Section 2.5
    a <- rep(0, p)
    for (j in 1:p) {
      qtl <- qnorm(t[l]/2) # 
      if (inv_a[j] > 0) {
        a[j] <- pnorm((qtl + bW.est[j])/inv_a[j]) + 
          pnorm((qtl - bW.est[j])/inv_a[j])
      }
      else {
        a[j] <- as.numeric(abs(bW.est[j]) > abs(qtl))
      }
    }
    Vt[l] <- min(sum(a), Rt[l])
    if (Rt[l] == 0) {
      FDPt[l] <- 0
    }
    else {
      FDPt[l] <- Vt[l]/Rt[l]
    }
  }
############################################################################# 
# The dependency-adjusted p-values, i.e., Eq. (18) in the article.  
############################################################################# 
  adj.P <- as.vector(rep(0, p))
  for (j in 1:p) {
    if (inv_a[j] > 0) {
      adj.P[j] <- 2 * (1 - pnorm(abs(Z[j] - bW.est[j])/inv_a[j]))
    }
    else {
      adj.P[j] <- as.numeric(abs(Z[j] - bW.est[j]) == 0)
    }
  }
  sort <- sort(adj.P, index.return = TRUE)
  adj.index <- sort$ix
  adj.P <- sort$x
  Pvals <- data.frame(p.value = P, Index = index)
  adjPvals <- data.frame(p.value = adj.P, Index = adj.index)
  if (t.default) {
    FDPvals <- data.frame(minus.logt = -log(t), rejects = Rt, 
                          false.rejects = Vt, FDP = FDPt)
  }
  else {
    FDPvals <- data.frame(t = t, rejects = Rt, false.rejects = Vt, 
                          FDP = FDPt)
  }
  results <- list(Pvalue = Pvals, adjPvalue = adjPvals, FDP = FDPvals, 
                  pi0 = p0.est/p, K = K, sigma = NULL)
  
  class(results) <- "FDPresults"
  
  
  return(results)
}
