###############################################################################
# compare_models.R - reusable model-comparison pipeline for binomial data
#   Models: Binomial, Beta-Binomial, Logit-Normal-Binomial, Kumaraswamy-Binomial,
#           Bernstein-Binomial (K grid 1..Kmax)
#   Criteria: logLik, AIC, BIC, L1/MSE vs empirical pmf, Pearson GOF (bootstrap)
#   Usage:
#     source('compare_models.R')
#     res <- compare_one(N, R, label='my data', Kmax=15, B=199)
#     compare_many(list(name = list(N=N, R=R), ...), out_csv='screening.csv', B=199)
###############################################################################

# ---------------- core functions ----------------
a_mat <- function(N, R, K) {
  k <- 0:K
  sapply(k, function(kk) exp(lbeta(kk + 1 + R, K - kk + 1 + N - R) - lbeta(kk + 1, K - kk + 1)))
}
bern_pmf_obs <- function(N, R, lambda) {
  K <- length(lambda) - 1
  A <- a_mat(N, R, K)
  exp(lchoose(N, R)) * as.numeric(A %*% lambda)
}
ll_bern <- function(lambda, N, R) sum(log(bern_pmf_obs(N, R, lambda)))
bern_mm <- function(N, R, K, lambda_init = NULL, ep = 1e-10, maxit = 2000) {
  if (is.null(lambda_init)) lambda <- rep(1 / (K + 1), K + 1)
  A <- a_mat(N, R, K)
  lam <- lambda; ll <- ll_bern(lam, N, R)
  for (it in 1:maxit) {
    prop <- A * matrix(lam, nrow = length(R), ncol = K + 1, byrow = TRUE)
    gamma <- prop / rowSums(prop)
    lam_new <- colMeans(gamma)
    ll_new <- ll_bern(lam_new, N, R)
    if (abs(ll_new - ll) < ep * (abs(ll) + 1e-12)) break
    lam <- lam_new; ll <- ll_new
  }
  list(lambda = lam, ll = ll, iter = it)
}

# collapsed Bern-Bino fitter for constant m (counts over x = 0..m)
bern_mm_tab <- function(m, counts, K, lambda_init = NULL, ep = 1e-10, maxit = 2000) {
  n <- sum(counts)
  k <- 0:K; x <- 0:m
  A <- sapply(k, function(kk) exp(lbeta(kk + 1 + x, K - kk + 1 + m - x) - lbeta(kk + 1, K - kk + 1)))
  lc <- lchoose(m, x)
  if (is.null(lambda_init)) lam <- rep(1 / (K + 1), K + 1) else lam <- lambda_init
  ll <- sum(counts * (lc + log(as.numeric(A %*% lam))))
  for (it in 1:maxit) {
    prop <- A * matrix(lam, nrow = m + 1, ncol = K + 1, byrow = TRUE)
    gamma <- prop / rowSums(prop)
    lam_new <- as.numeric(t(counts) %*% gamma) / n
    ll_new <- sum(counts * (lc + log(as.numeric(A %*% lam_new))))
    if (abs(ll_new - ll) < ep * (abs(ll) + 1e-12)) break
    lam <- lam_new; ll <- ll_new
  }
  list(lambda = lam, ll = ll, iter = it)
}

ll_binom <- function(p, N, R) sum(dbinom(R, N, p, log = TRUE))

ll_bb <- function(theta, N, R) {
  a <- exp(theta[1]); b <- exp(theta[2])
  sum(lchoose(N, R) + lbeta(a + R, b + N - R) - lbeta(a, b))
}

f_ln <- function(p, mu, sigma) {
  out <- rep(0, length(p)); ok <- p > 0 & p < 1; q <- p[ok]
  out[ok] <- exp(-(log(q / (1 - q)) - mu)^2 / (2 * sigma^2)) / (sqrt(2 * pi) * sigma * q * (1 - q))
  out
}
int_ln <- function(Rv, Nv, mu, sigma, ppow = 0) {
  f <- function(u) { q <- plogis(u); q^(Rv + ppow) * (1 - q)^(Nv - Rv) * dnorm(u, mu, sigma) }
  v <- tryCatch(integrate(f, -20, 20, rel.tol = 1e-8, subdivisions = 400,
                          stop.on.error = FALSE)$value, error = function(e) NA)
  if (!is.finite(v) || v <= 0) 1e-300 else v
}
ll_ln <- function(theta, N, R) {
  mu <- theta[1]; sigma <- exp(theta[2])
  sum(log(sapply(seq_along(R), function(i) int_ln(R[i], N[i], mu, sigma, 0)))) + sum(lchoose(N, R))
}

f_km <- function(p, a, b) {
  out <- rep(0, length(p)); ok <- p > 0 & p < 1; q <- p[ok]
  out[ok] <- a * b * q^(a - 1) * (1 - q^a)^(b - 1)
  out
}
int_km <- function(Rv, Nv, a, b) {
  f <- function(u) u^(Rv / a) * (1 - u^(1 / a))^(Nv - Rv) * b * (1 - u)^(b - 1)
  v <- tryCatch(integrate(f, 0, 1, rel.tol = 1e-8, subdivisions = 400,
                          stop.on.error = FALSE)$value, error = function(e) NA)
  if (!is.finite(v) || v <= 0) 1e-300 else v
}
ll_km <- function(theta, N, R) {
  a <- exp(theta[1]); b <- exp(theta[2])
  sum(log(sapply(seq_along(R), function(i) int_km(R[i], N[i], a, b)))) + sum(lchoose(N, R))
}

# fast tabulated log-likelihood when all N are equal to m
ll_ln_tab <- function(theta, N, R, m, counts) {
  mu <- theta[1]; sigma <- exp(theta[2])
  px <- sapply(0:m, function(xx) int_ln(xx, m, mu, sigma, 0))
  sum(counts * (lchoose(m, 0:m) + log(px)))
}
ll_km_tab <- function(theta, N, R, m, counts) {
  a <- exp(theta[1]); b <- exp(theta[2])
  px <- sapply(0:m, function(xx) int_km(xx, m, a, b))
  sum(counts * (lchoose(m, 0:m) + log(px)))
}

pearson_stat <- function(N, R, mu1, mu2) {
  Var <- N * mu1 - N * mu2 + N * (N - 1) * (mu2 - mu1^2)
  sum((R - N * mu1)^2 / Var)
}

# ---------------- fit one model and return estimates ----------------
fit_model <- function(model, N, R, K = NULL) {
  n <- length(N)
  switch(model,
    bin = { p <- sum(R) / sum(N); list(ll = ll_binom(p, N, R), npar = 1,
                                       extra = list(p = p)) },
    bb  = { opt <- optim(c(0, 0), function(th) -ll_bb(th, N, R), method = 'L-BFGS-B',
                         lower = c(-8, -8), upper = c(8, 8))
            a <- exp(opt$par[1]); b <- exp(opt$par[2])
            list(ll = -opt$value, npar = 2, extra = list(a = a, b = b)) },
    ln  = { if (length(unique(N)) == 1) {
              m <- N[1]; counts <- tabulate(R + 1, nbins = m + 1)
              opt <- optim(c(0, log(0.5)), function(th) -ll_ln_tab(th, N, R, m, counts),
                           method = 'L-BFGS-B', lower = c(-5, log(0.05)), upper = c(5, log(3)))
            } else {
              opt <- optim(c(0, log(0.5)), function(th) -ll_ln(th, N, R), method = 'L-BFGS-B',
                           lower = c(-5, log(0.05)), upper = c(5, log(3)))
            }
            mu <- opt$par[1]; sig <- exp(opt$par[2])
            list(ll = -opt$value, npar = 2, extra = list(mu = mu, sigma = sig)) },
    km  = { if (length(unique(N)) == 1) {
              m <- N[1]; counts <- tabulate(R + 1, nbins = m + 1)
              opt <- optim(c(0, 0), function(th) -ll_km_tab(th, N, R, m, counts),
                           method = 'L-BFGS-B', lower = c(-8, -8), upper = c(8, 8))
            } else {
              opt <- optim(c(0, 0), function(th) -ll_km(th, N, R), method = 'L-BFGS-B',
                           lower = c(-8, -8), upper = c(8, 8))
            }
            a <- exp(opt$par[1]); b <- exp(opt$par[2])
            list(ll = -opt$value, npar = 2, extra = list(a = a, b = b)) },
    bern = { fit <- bern_mm(N, R, K)
             list(ll = fit$ll, npar = K + 1, extra = list(lambda = fit$lambda, K = K)) }
  )
}

moments_of <- function(model, fit, N, R, K = NULL) {
  ex <- fit$extra
  switch(model,
    bin = c(mu1 = ex$p, mu2 = ex$p^2),
    bb  = c(mu1 = ex$a / (ex$a + ex$b),
            mu2 = ex$a * (ex$a + 1) / ((ex$a + ex$b) * (ex$a + ex$b + 1))),
    ln  = c(mu1 = int_ln(0, 0, ex$mu, ex$sigma, 1),
            mu2 = int_ln(0, 0, ex$mu, ex$sigma, 2)),
    km  = c(mu1 = ex$b * beta(1 + 1 / ex$a, ex$b),
            mu2 = ex$b * beta(1 + 2 / ex$a, ex$b)),
    bern = { k <- 0:K
             c(mu1 = sum(ex$lambda * (k + 1) / (K + 2)),
               mu2 = sum(ex$lambda * (k + 1) * (k + 2) / ((K + 2) * (K + 3)))) }
  )
}

pmf_of <- function(model, fit, Nv, x) {
  ex <- fit$extra
  switch(model,
    bin = dbinom(x, Nv, ex$p),
    bb  = exp(lchoose(Nv, x) + lbeta(ex$a + x, ex$b + Nv - x) - lbeta(ex$a, ex$b)),
    ln  = sapply(x, function(xx) exp(lchoose(Nv, xx)) * int_ln(xx, Nv, ex$mu, ex$sigma, 0)),
    km  = sapply(x, function(xx) exp(lchoose(Nv, xx)) * int_km(xx, Nv, ex$a, ex$b)),
    bern = bern_pmf_obs(rep(Nv, length(x)), x, ex$lambda)
  )
}

sim_of <- function(model, fit, N, B = 199, seed = 12345) {
  ex <- fit$extra
  set.seed(seed)
  switch(model,
    bin = replicate(B, rbinom(length(N), N, ex$p)),
    bb  = replicate(B, { pp <- rbeta(length(N), ex$a, ex$b); rbinom(length(N), N, pp) }),
    ln  = replicate(B, { z <- rnorm(length(N), ex$mu, ex$sigma); rbinom(length(N), N, plogis(z)) }),
    km  = replicate(B, { u <- runif(length(N)); pp <- (1 - (1 - u)^(1 / ex$b))^(1 / ex$a); rbinom(length(N), N, pp) }),
    bern = replicate(B, { K <- ex$K; kk <- sample(0:K, length(N), replace = TRUE, prob = ex$lambda)
                          pp <- rbeta(length(N), kk + 1, K - kk + 1); rbinom(length(N), N, pp) })
  )
}

# empirical vs fitted discrepancy (L1 = half TV, MSE), on litter sizes with >= min_cnt units
disc_of <- function(model, fit, N, R, min_cnt = 5) {
  tab <- table(N)
  Ns <- as.integer(names(tab)[tab >= min_cnt])
  if (length(Ns) == 0) Ns <- as.integer(names(sort(tab, decreasing = TRUE))[1])
  L1 <- 0; MSE <- 0; wsum <- 0
  for (Nv in Ns) {
    idx <- which(N == Nv)
    obs <- as.numeric(table(factor(R[idx], levels = 0:Nv))) / length(idx)
    pred <- pmf_of(model, fit, Nv, 0:Nv)
    pred <- pred / sum(pred)
    w <- length(idx)
    L1  <- L1 + w * 0.5 * sum(abs(obs - pred))
    MSE <- MSE + w * mean((obs - pred)^2)
    wsum <- wsum + w
  }
  c(L1 = L1 / wsum, MSE = MSE / wsum)
}

# ---------------- main comparison ----------------
compare_one <- function(N, R, label = 'data', Kmax = 15, B = 199, seed = 12345) {
  n <- length(N)
  if (any(R > N) || any(R < 0)) stop('invalid binomial counts')
  if (length(unique(N)) == 1) return(compare_one_fastm(N, R, label, Kmax, B, seed))
  
  mods <- c('bin', 'bb', 'ln', 'km')
  fits <- lapply(mods, function(m) fit_model(m, N, R))
  names(fits) <- mods

  K_grid <- 1:Kmax
  bern_ll <- sapply(K_grid, function(k) fit_model('bern', N, R, k)$ll)
  AIC_bern <- 2 * (K_grid + 1) - 2 * bern_ll
  BIC_bern <- (K_grid + 1) * log(n) - 2 * bern_ll
  K_bic  <- K_grid[which.min(BIC_bern)]
  K_ll   <- K_grid[which.max(bern_ll)]
  fit_bern_bic <- fit_model('bern', N, R, K_bic)
  fit_bern_ll  <- fit_model('bern', N, R, K_ll)

  rows <- list()
  for (m in mods) {
    fit <- fits[[m]]
    mom <- moments_of(m, fit, N, R)
    X2 <- pearson_stat(N, R, mom['mu1'], mom['mu2'])
    sim <- sim_of(m, fit, N, B, seed)
    X2s <- apply(sim, 2, function(RR) pearson_stat(N, RR, mom['mu1'], mom['mu2']))
    pval <- (1 + sum(X2s >= X2)) / (B + 1)
    d <- disc_of(m, fit, N, R)
    rows[[m]] <- data.frame(Model = m, K = NA, npar = fit$npar, logLik = fit$ll,
                            AIC = 2 * fit$npar - 2 * fit$ll, BIC = fit$npar * log(n) - 2 * fit$ll,
                            L1 = d['L1'], MSE = d['MSE'], X2 = X2, pval = pval)
  }
  for (tag in c('bic', 'll')) {
    fit <- if (tag == 'bic') fit_bern_bic else fit_bern_ll
    K <- fit$extra$K
    mom <- moments_of('bern', fit, N, R, K)
    X2 <- pearson_stat(N, R, mom['mu1'], mom['mu2'])
    sim <- sim_of('bern', fit, N, B, seed)
    X2s <- apply(sim, 2, function(RR) pearson_stat(N, RR, mom['mu1'], mom['mu2']))
    pval <- (1 + sum(X2s >= X2)) / (B + 1)
    d <- disc_of('bern', fit, N, R)
    rows[[paste0('bern_', tag)]] <- data.frame(Model = paste0('bern_', tag), K = K,
      npar = fit$npar, logLik = fit$ll, AIC = 2 * fit$npar - 2 * fit$ll,
      BIC = fit$npar * log(n) - 2 * fit$ll, L1 = d['L1'], MSE = d['MSE'], X2 = X2, pval = pval)
  }
  tab <- do.call(rbind, rows)
  tab$label <- label
  tab
}

compare_many <- function(datasets, out_csv = NULL, Kmax = 15, B = 199, seed = 12345) {
  res <- lapply(names(datasets), function(nm) {
    d <- datasets[[nm]]
    tryCatch(compare_one(d$N, d$R, label = nm, Kmax = Kmax, B = B, seed = seed),
             error = function(e) data.frame(Model = 'ERROR', label = nm, err = conditionMessage(e)))
  })
  out <- do.call(rbind, res)
  if (!is.null(out_csv)) write.csv(out, out_csv, row.names = FALSE)
  out
}

# fast comparison for constant m (collapsed counts)
compare_one_fastm <- function(N, R, label = 'data', Kmax = 15, B = 199, seed = 12345) {
  m <- N[1]; n <- length(N)
  counts <- tabulate(R + 1, nbins = m + 1)
  x <- 0:m
  # fitted pmfs at x=0..m
  pmf_bin <- dbinom(x, m, sum(R) / (n * m))
  fit_bb <- optim(c(0, 0), function(th) {
    a <- exp(th[1]); b <- exp(th[2])
    -sum(counts * (lchoose(m, x) + lbeta(a + x, b + m - x) - lbeta(a, b)))
  }, method = 'L-BFGS-B', lower = c(-8, -8), upper = c(8, 8))
  a_bb <- exp(fit_bb$par[1]); b_bb <- exp(fit_bb$par[2])
  pmf_bb <- exp(lchoose(m, x) + lbeta(a_bb + x, b_bb + m - x) - lbeta(a_bb, b_bb))
  fit_ln <- optim(c(0, log(0.5)), function(th) -ll_ln_tab(th, N, R, m, counts),
                  method = 'L-BFGS-B', lower = c(-5, log(0.05)), upper = c(5, log(3)))
  mu_ln <- fit_ln$par[1]; sig_ln <- exp(fit_ln$par[2])
  pmf_ln <- sapply(x, function(xx) exp(lchoose(m, xx)) * int_ln(xx, m, mu_ln, sig_ln, 0))
  fit_km <- optim(c(0, 0), function(th) -ll_km_tab(th, N, R, m, counts),
                  method = 'L-BFGS-B', lower = c(-8, -8), upper = c(8, 8))
  a_km <- exp(fit_km$par[1]); b_km <- exp(fit_km$par[2])
  pmf_km <- sapply(x, function(xx) exp(lchoose(m, xx)) * int_km(xx, m, a_km, b_km))

  K_grid <- 1:Kmax
  bern_ll <- sapply(K_grid, function(k) bern_mm_tab(m, counts, k)$ll)
  AIC_bern <- 2 * (K_grid + 1) - 2 * bern_ll
  BIC_bern <- (K_grid + 1) * log(n) - 2 * bern_ll
  K_bic <- K_grid[which.min(BIC_bern)]; K_ll <- K_grid[which.max(bern_ll)]
  lam_bic <- bern_mm_tab(m, counts, K_bic)$lambda
  lam_ll  <- bern_mm_tab(m, counts, K_ll)$lambda
  pmf_bern <- function(lam, K) {
    k <- 0:K
    A <- sapply(k, function(kk) exp(lbeta(kk + 1 + x, K - kk + 1 + m - x) - lbeta(kk + 1, K - kk + 1)))
    exp(lchoose(m, x)) * as.numeric(A %*% lam)
  }
  emp <- counts / n
  # Pearson from counts: Var(X)=m*mu1 - m*mu2 + m*(m-1)*(mu2-mu1^2)
  pearson_counts <- function(pmf) {
    mu1 <- sum(x * pmf) / m
    mu2 <- sum(x^2 * pmf) / m^2
    Var <- m * mu1 - m * mu2 + m * (m - 1) * (mu2 - mu1^2)
    sum(counts * (x - m * mu1)^2 / Var)
  }
  gof_p <- function(pmf) {
    set.seed(seed)
    X2 <- pearson_counts(pmf)
    X2s <- replicate(B, {
      cc <- as.numeric(rmultinom(1, n, pmf))
      mu1 <- sum(x * cc) / (n * m); mu2 <- sum(x^2 * cc) / (n * m^2)
      Var <- m * mu1 - m * mu2 + m * (m - 1) * (mu2 - mu1^2)
      sum(cc * (x - m * mu1)^2 / Var)
    })
    (1 + sum(X2s >= X2)) / (B + 1)
  }
  L1 <- function(pmf) 0.5 * sum(abs(emp - pmf))
  MSE <- function(pmf) mean((emp - pmf)^2)

  models <- c('bin', 'bb', 'ln', 'km')
  pmfs <- list(bin = pmf_bin, bb = pmf_bb, ln = pmf_ln, km = pmf_km)
  lls <- c(bin = sum(counts * log(pmf_bin)), bb = -fit_bb$value, ln = -fit_ln$value, km = -fit_km$value)
  rows <- lapply(models, function(mm) {
    data.frame(Model = mm, K = NA, npar = ifelse(mm == 'bin', 1, 2), logLik = lls[[mm]],
               AIC = 2 * (ifelse(mm == 'bin', 1, 2)) - 2 * lls[[mm]],
               BIC = (ifelse(mm == 'bin', 1, 2)) * log(n) - 2 * lls[[mm]],
               L1 = L1(pmfs[[mm]]), MSE = MSE(pmfs[[mm]]), X2 = pearson_counts(pmfs[[mm]]),
               pval = gof_p(pmfs[[mm]]))
  })
  rows[['bern_bic']] <- data.frame(Model = 'bern_bic', K = K_bic, npar = K_bic + 1, logLik = bern_ll[K_bic],
    AIC = AIC_bern[K_bic], BIC = BIC_bern[K_bic], L1 = L1(pmf_bern(lam_bic, K_bic)),
    MSE = MSE(pmf_bern(lam_bic, K_bic)), X2 = pearson_counts(pmf_bern(lam_bic, K_bic)),
    pval = gof_p(pmf_bern(lam_bic, K_bic)))
  rows[['bern_ll']] <- data.frame(Model = 'bern_ll', K = K_ll, npar = K_ll + 1, logLik = bern_ll[K_ll],
    AIC = AIC_bern[K_ll], BIC = BIC_bern[K_ll], L1 = L1(pmf_bern(lam_ll, K_ll)),
    MSE = MSE(pmf_bern(lam_ll, K_ll)), X2 = pearson_counts(pmf_bern(lam_ll, K_ll)),
    pval = gof_p(pmf_bern(lam_ll, K_ll)))
  tab <- do.call(rbind, rows)
  tab$label <- label
  tab
}

# winner summary: best model by logLik and by L1/MSE
winner_summary <- function(tab) {
  sub <- tab[!grepl('^bern_', tab$Model), ]
  bern <- tab[tab$Model == 'bern_ll', ]
  cat(sprintf('best logLik (parametric): %s (%.2f) | Bern-Bino(K=%d) logLik %.2f | Bern wins logLik: %s\n',
              sub$Model[which.max(sub$logLik)], max(sub$logLik), bern$K, bern$logLik,
              ifelse(bern$logLik > max(sub$logLik), 'YES', 'no')))
  cat(sprintf('best L1: %s (%.4f) | Bern L1 %.4f\n', sub$Model[which.min(sub$L1)], min(sub$L1), bern$L1))
  cat(sprintf('best MSE: %s (%.5f) | Bern MSE %.5f\n', sub$Model[which.min(sub$MSE)], min(sub$MSE), bern$MSE))
  cat(sprintf('best GOF p: %s (%.3f) | Bern p %.3f\n', sub$Model[which.max(sub$pval)], max(sub$pval), bern$pval))
}
