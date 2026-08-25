###############################################################################
# Real-data analysis: low-iron rat teratology data (58 litters)
# Section 5 of the paper. Bern-Bino vs Binomial / Beta-Binomial /
# Logit-Normal-Binomial / Kumaraswamy-Binomial, plus Bayesian DA inference.
#
# Run from the folder containing lirat.csv:
#   Rscript teratology_analysis.R
###############################################################################

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(coda)
  library(gridExtra)
})

# ---------------- paths (relative to working directory) ----------------
base_dir <- normalizePath(getwd())
repo_root <- normalizePath(file.path(base_dir, '..', '..'))
data_path <- file.path(base_dir, 'lirat.csv')
fig_dir   <- file.path(repo_root, 'Thesis', 'figures')
out_dir   <- file.path(base_dir, 'output')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(12345)

save_png <- function(p, file, w = 7, h = 5) {
  ggsave(file.path(fig_dir, file), p, width = w, height = h, dpi = 300, bg = 'white')
}
save_grob <- function(plots, file, w = 10, h = 4.5, ncol = 2) {
  png(file.path(fig_dir, file), width = w, height = h, units = 'in', res = 300, bg = 'white')
  do.call(gridExtra::grid.arrange, c(plots, list(ncol = ncol)))
  dev.off()
}

# ============================================================================
# 1. Data description and overdispersion check
# ============================================================================
d <- read.csv(data_path)
d$grp <- factor(d$grp)
n <- nrow(d)
N <- d$N; R <- d$R
p_hat <- sum(R) / sum(N)

grp_tab <- d %>% group_by(grp) %>%
  summarise(litters = n(),
            mean_N = mean(N),
            sum_N = sum(N),
            sum_R = sum(R),
            prop = sum(R) / sum(N),
            mean_hb = mean(hb), .groups = 'drop')
print(grp_tab)

# extra-binomial dispersion statistic (Binomial model)
phi_hat <- sum((R - N * p_hat)^2 / (N * p_hat * (1 - p_hat))) / (n - 1)
disp_p  <- pchisq((n - 1) * phi_hat, df = n - 1, lower.tail = FALSE)

cat(sprintf('n = %d,  sum N = %d,  sum R = %d\n', n, sum(N), sum(R)))
cat(sprintf('p_hat = %.4f,  phi_hat = %.3f,  chi2(%d) p-value = %.4g\n',
            p_hat, phi_hat, n - 1, disp_p))

write.csv(as.data.frame(grp_tab), file.path(out_dir, 'group_summary.csv'), row.names = FALSE)

# ---- plot16: per-litter death proportion by group + histogram ----
p16a <- ggplot(d, aes(x = grp, y = R / N, fill = grp)) +
  geom_boxplot(alpha = 0.35, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.75) +
  scale_fill_brewer(palette = 'Set2') +
  labs(x = 'Group', y = 'Proportion of dead fetuses per litter',
       title = '(a) Per-litter death proportion by group') +
  theme_bw() + theme(legend.position = 'none')
p16b <- ggplot(d, aes(x = R / N)) +
  geom_histogram(bins = 12, fill = 'steelblue', color = 'white') +
  geom_vline(xintercept = p_hat, linetype = 2, color = 'firebrick') +
  labs(x = 'Proportion of dead fetuses', y = 'Number of litters',
       title = '(b) Overall distribution of death proportions') +
  theme_bw()
save_grob(list(p16a, p16b), 'plot16.png')

# ============================================================================
# 2. Model log-likelihoods and MLEs
# ============================================================================
# --- Binomial ---
ll_binom <- function(p, N, R) sum(dbinom(R, N, p, log = TRUE))
ll_bin <- ll_binom(p_hat, N, R)

# --- Beta-Binomial ---
ll_bb <- function(theta, N, R) {
  a <- exp(theta[1]); b <- exp(theta[2])
  sum(lchoose(N, R) + lbeta(a + R, b + N - R) - lbeta(a, b))
}
opt_bb <- optim(c(0, 0), function(th) -ll_bb(th, N, R), method = 'L-BFGS-B',
                lower = c(-8, -8), upper = c(8, 8))
alpha_bb <- exp(opt_bb$par[1]); beta_bb <- exp(opt_bb$par[2]); ll_bb_v <- -opt_bb$value

# --- Logit-Normal-Binomial ---
f_ln <- function(p, mu, sigma) {
  out <- rep(0, length(p)); ok <- p > 0 & p < 1; q <- p[ok]
  out[ok] <- exp(-(log(q / (1 - q)) - mu)^2 / (2 * sigma^2)) / (sqrt(2 * pi) * sigma * q * (1 - q))
  out
}
int_ln <- function(Rv, Nv, mu, sigma, ppow = 0) {
  f <- function(u) {
    q <- plogis(u)
    q^(Rv + ppow) * (1 - q)^(Nv - Rv) * dnorm(u, mu, sigma)
  }
  v <- tryCatch(integrate(f, -20, 20, rel.tol = 1e-8, subdivisions = 400,
                          stop.on.error = FALSE)$value, error = function(e) NA)
  if (!is.finite(v) || v <= 0) 1e-300 else v
}
ll_ln <- function(theta, N, R) {
  mu <- theta[1]; sigma <- exp(theta[2])
  s <- 0
  for (i in seq_along(R)) {
    s <- s + log(int_ln(R[i], N[i], mu, sigma, 0))
  }
  s + sum(lchoose(N, R))
}
opt_ln <- optim(c(0, log(0.5)), function(th) -ll_ln(th, N, R), method = 'L-BFGS-B',
                lower = c(-5, log(0.05)), upper = c(5, log(3)))
mu_ln <- opt_ln$par[1]; sigma_ln <- exp(opt_ln$par[2]); ll_ln_v <- -opt_ln$value

# --- Kumaraswamy-Binomial ---
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
  s <- 0
  for (i in seq_along(R)) {
    s <- s + log(int_km(R[i], N[i], a, b))
  }
  s + sum(lchoose(N, R))
}
opt_km <- optim(c(0, 0), function(th) -ll_km(th, N, R), method = 'L-BFGS-B',
                lower = c(-8, -8), upper = c(8, 8))
a_km <- exp(opt_km$par[1]); b_km <- exp(opt_km$par[2]); ll_km_v <- -opt_km$value

# --- Bernstein-Binomial (MM algorithm, eq. 2.8 with unit-specific N_i) ---
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
  lam <- lambda
  ll <- ll_bern(lam, N, R)
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

# ---- K selection (BIC primary, AIC secondary) ----
K_grid <- 1:15
res_K <- data.frame(K = K_grid, ll = NA, AIC = NA, BIC = NA)
bern_fits <- vector('list', length(K_grid))
for (j in seq_along(K_grid)) {
  k <- K_grid[j]
  fit <- bern_mm(N, R, k)
  bern_fits[[j]] <- fit
  p <- k + 1
  res_K$ll[j]  <- fit$ll
  res_K$AIC[j] <- 2 * p - 2 * fit$ll
  res_K$BIC[j] <- p * log(n) - 2 * fit$ll
}
K_bic <- res_K$K[which.min(res_K$BIC)]
K_aic <- res_K$K[which.min(res_K$AIC)]
cat(sprintf('K by BIC = %d,  K by AIC = %d\n', K_bic, K_aic))
write.csv(res_K, file.path(out_dir, 'k_selection.csv'), row.names = FALSE)

p18a <- ggplot(res_K, aes(x = K, y = BIC)) + geom_line() + geom_point() +
  geom_point(data = res_K[which.min(res_K$BIC), ], color = 'firebrick', size = 3) +
  labs(x = 'K', y = 'BIC', title = '(a) BIC vs K') + theme_bw()
p18b <- ggplot(res_K, aes(x = K, y = AIC)) + geom_line() + geom_point() +
  geom_point(data = res_K[which.min(res_K$AIC), ], color = 'firebrick', size = 3) +
  labs(x = 'K', y = 'AIC', title = '(b) AIC vs K') + theme_bw()
save_grob(list(p18a, p18b), 'plot18.png')

# main Bern-Bino fit at K_bic
fit_bern <- bern_mm(N, R, K_bic)
lambda_hat <- fit_bern$lambda
ll_bern_v <- fit_bern$ll
cat('Bern-Bino lambda_hat:', round(lambda_hat, 4), '\n')

# ============================================================================
# 3. Model comparison: moments, Pearson GOF (parametric bootstrap), AIC/BIC, MAE
# ============================================================================
moments <- function(model) {
  switch(model,
    bin  = c(mu1 = p_hat, mu2 = p_hat^2),
    bb   = c(mu1 = alpha_bb / (alpha_bb + beta_bb),
             mu2 = alpha_bb * (alpha_bb + 1) / ((alpha_bb + beta_bb) * (alpha_bb + beta_bb + 1))),
    ln   = c(mu1 = int_ln(0, 0, mu_ln, sigma_ln, 1),
             mu2 = int_ln(0, 0, mu_ln, sigma_ln, 2)),
    km   = c(mu1 = b_km * beta(1 + 1 / a_km, b_km),
             mu2 = b_km * beta(1 + 2 / a_km, b_km)),
    bern = {
      K <- K_bic; k <- 0:K
      c(mu1 = sum(lambda_hat * (k + 1) / (K + 2)),
        mu2 = sum(lambda_hat * (k + 1) * (k + 2) / ((K + 2) * (K + 3))))
    }
  )
}

pearson_stat <- function(N, R, mu1, mu2) {
  Var <- N * mu1 - N * mu2 + N * (N - 1) * (mu2 - mu1^2)
  sum((R - N * mu1)^2 / Var)
}

sim_fun <- function(model) {
  switch(model,
    bin  = function(N) rbinom(length(N), N, p_hat),
    bb   = function(N) { pp <- rbeta(length(N), alpha_bb, beta_bb); rbinom(length(N), N, pp) },
    ln   = function(N) { z <- rnorm(length(N), mu_ln, sigma_ln); pp <- plogis(z); rbinom(length(N), N, pp) },
    km   = function(N) { u <- runif(length(N)); pp <- (1 - (1 - u)^(1 / b_km))^(1 / a_km); rbinom(length(N), N, pp) },
    bern = function(N) {
      K <- K_bic
      kk <- sample(0:K, length(N), replace = TRUE, prob = lambda_hat)
      pp <- rbeta(length(N), kk + 1, K - kk + 1)
      rbinom(length(N), N, pp)
    }
  )
}

gof_boot <- function(model, B = 999, seed = 12345) {
  set.seed(seed)
  m <- moments(model); sf <- sim_fun(model)
  X2 <- pearson_stat(N, R, m['mu1'], m['mu2'])
  X2s <- replicate(B, pearson_stat(N, sf(N), m['mu1'], m['mu2']))
  c(X2 = X2, pval = (1 + sum(X2s >= X2)) / (B + 1))
}

pmf_fun <- function(model) {
  switch(model,
    bin  = function(Nv, x) dbinom(x, Nv, p_hat),
    bb   = function(Nv, x) exp(lchoose(Nv, x) + lbeta(alpha_bb + x, beta_bb + Nv - x) - lbeta(alpha_bb, beta_bb)),
    ln   = function(Nv, x) sapply(x, function(xx) {
      exp(lchoose(Nv, xx)) * int_ln(xx, Nv, mu_ln, sigma_ln, 0) }),
    km   = function(Nv, x) sapply(x, function(xx) {
      exp(lchoose(Nv, xx)) * int_km(xx, Nv, a_km, b_km) }),
    bern = function(Nv, x) bern_pmf_obs(rep(Nv, length(x)), x, lambda_hat)
  )
}

mae_dist <- function(model, min_cnt = 5) {
  tab <- table(N)
  Ns <- as.integer(names(tab)[tab >= min_cnt])
  if (length(Ns) == 0) return(NA)
  pf <- pmf_fun(model)
  mae <- 0; wsum <- 0
  for (Nv in Ns) {
    idx <- which(N == Nv)
    obs <- as.numeric(table(factor(R[idx], levels = 0:Nv))) / length(idx)
    pred <- pf(Nv, 0:Nv)
    mae <- mae + length(idx) * mean(abs(obs - pred))
    wsum <- wsum + length(idx)
  }
  mae / wsum
}

models <- c('bin', 'bb', 'ln', 'km', 'bern')
lls  <- c(ll_bin, ll_bb_v, ll_ln_v, ll_km_v, ll_bern_v)
npar <- c(1, 2, 2, 2, K_bic + 1)
AICs <- 2 * npar - 2 * lls
BICs <- npar * log(n) - 2 * lls
gofs <- sapply(models, gof_boot)
maes <- sapply(models, mae_dist)
rankA <- rank(AICs)
rankp <- rank(-gofs['pval', ])

comp <- data.frame(
  Model = c('Binomial', 'Beta-Binomial', 'Logit-Normal', 'Kumaraswamy',
            sprintf('Bern-Bino (K=%d)', K_bic)),
  npar = npar, logLik = lls, AIC = AICs, BIC = BICs,
  X2 = gofs['X2', ], pval = gofs['pval', ],
  MAE = maes, Rank_A = rankA, Rank_p = rankp)
comp[, 2:9] <- round(comp[, 2:9], 4)
print(comp)
write.csv(comp, file.path(out_dir, 'model_comparison.csv'), row.names = FALSE)

# ---- plot17: fitted pmf (modal litter size) + Pearson residuals ----
modal_N <- as.integer(names(sort(table(N), decreasing = TRUE))[1])
xgrid <- 0:modal_N
obs_pmf <- as.numeric(table(factor(R[N == modal_N], levels = xgrid))) / sum(N == modal_N)
pmf_df <- data.frame(x = xgrid, obs = obs_pmf)
for (m in models) {
  pmf_df[[m]] <- pmf_fun(m)(modal_N, xgrid)
}
pmf_long <- pmf_df %>% pivot_longer(cols = all_of(c('obs', models)), names_to = 'model', values_to = 'prob')
pmf_long$model <- factor(pmf_long$model,
                         levels = c('obs', 'bin', 'bb', 'ln', 'km', 'bern'),
                         labels = c('Observed', 'Binomial', 'Beta-Binomial',
                                    'Logit-Normal', 'Kumaraswamy', 'Bern-Bino'))
p17a <- ggplot(pmf_long, aes(x = x, y = prob, color = model, shape = model)) +
  geom_line() + geom_point(size = 2) +
  labs(x = paste0('Number of dead fetuses (litter size N = ', modal_N, ')'),
       y = 'Probability', title = paste0('(a) Fitted pmf at the modal litter size N = ', modal_N)) +
  theme_bw() + theme(legend.position = 'bottom')

mu1_bern <- moments('bern')['mu1']; mu2_bern <- moments('bern')['mu2']
Var_bern <- N * mu1_bern - N * mu2_bern + N * (N - 1) * (mu2_bern - mu1_bern^2)
pres <- (R - N * mu1_bern) / sqrt(Var_bern)
res_df <- data.frame(N = N, grp = d$grp, pres = pres)
p17b <- ggplot(res_df, aes(x = N, y = pres, color = grp)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 2) +
  labs(x = 'Litter size N', y = 'Pearson residual (Bern-Bino)',
       title = '(b) Pearson residuals vs litter size') +
  theme_bw() + theme(legend.position = 'bottom')
save_grob(list(p17a, p17b), 'plot17.png', w = 11, h = 4.8)

# ============================================================================
# 4. Estimated prior densities of p (Bernstein vs Beta vs LN vs KM)
# ============================================================================
pgrid <- seq(0.005, 0.995, length.out = 200)
dens_bern <- sapply(pgrid, function(p) {
  k <- 0:K_bic
  sum(lambda_hat * dbeta(p, k + 1, K_bic - k + 1))
})
dens_bb  <- dbeta(pgrid, alpha_bb, beta_bb)
dens_ln  <- f_ln(pgrid, mu_ln, sigma_ln)
dens_km  <- f_km(pgrid, a_km, b_km)
dens_df <- data.frame(p = pgrid, BernBino = dens_bern, Beta = dens_bb,
                      LogitNormal = dens_ln, Kumaraswamy = dens_km)
dens_long <- dens_df %>% pivot_longer(-p, names_to = 'model', values_to = 'density')
p19 <- ggplot(dens_long, aes(x = p, y = density, color = model)) +
  geom_line(linewidth = 1) +
  labs(x = 'p', y = 'Density',
       title = 'Estimated prior density of the fetal death probability p') +
  theme_bw() + theme(legend.position = 'bottom')
save_png(p19, 'plot19.png', w = 7, h = 5)

# ============================================================================
# 5. Bayesian inference via the DA algorithm (K = K_bic, Dirichlet(1,...,1))
# ============================================================================
rdirichlet <- function(alpha) { x <- rgamma(length(alpha), alpha); x / sum(x) }

da_one_chain <- function(N, R, K, alpha_prior, G, burnin, lambda_init, A) {
  lam <- lambda_init
  keep <- G - burnin
  out <- matrix(NA, keep, K + 1)
  for (g in 1:G) {
    # I-step: W_i | (R_i, N_i, lambda)
    logq <- log(lam) + log(A)              # n x (K+1)
    logq <- logq - apply(logq, 1, max)
    q <- exp(logq); q <- q / rowSums(q)
    W <- apply(q, 1, function(pr) sample(0:K, 1, prob = pr))
    Nk <- tabulate(W + 1, nbins = K + 1)
    # P-step: lambda | W
    lam <- rdirichlet(alpha_prior + Nk)
    if (g > burnin) out[g - burnin, ] <- lam
  }
  out
}

G <- 10000; burnin <- 2500; M <- 4
K <- K_bic
alpha_prior <- rep(1, K + 1)
A <- a_mat(N, R, K)
chains <- lapply(seq_len(M), function(m) {
  set.seed(1000 + m)
  da_one_chain(N, R, K, alpha_prior, G, burnin, rdirichlet(alpha_prior), A)
})
post <- do.call(rbind, chains)             # 30000 x (K+1)

post_mean <- colMeans(post)
post_sd   <- apply(post, 2, sd)
post_q <- t(apply(post, 2, quantile, probs = c(0.025, 0.975)))
ess <- effectiveSize(mcmc(post))
mcmc_list <- mcmc.list(lapply(chains, function(ch) mcmc(ch)))
rhat <- gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]

bay_tab <- data.frame(
  k = 0:K,
  post_mean = round(post_mean, 4), post_sd = round(post_sd, 4),
  q2.5 = round(post_q[, 1], 4), q97.5 = round(post_q[, 2], 4),
  ESS = round(ess, 1), Rhat = round(rhat, 4))
print(bay_tab)
write.csv(bay_tab, file.path(out_dir, 'bayesian_results.csv'), row.names = FALSE)

# posterior summaries of E[p] and Var(p)
Ep_post <- rowSums(post * matrix(0:K, nrow = nrow(post), ncol = K + 1, byrow = TRUE) + post) / (K + 2)
cat(sprintf('Posterior E[p]: mean %.4f, 95%% CI [%.4f, %.4f]\n',
            mean(Ep_post), quantile(Ep_post, 0.025), quantile(Ep_post, 0.975)))

# ---- plot20: posterior mean prior density with 95% pointwise band ----
sub_idx <- seq(1, nrow(post), by = 10)     # 3000 draws
psub <- post[sub_idx, , drop = FALSE]
pgrid2 <- seq(0.005, 0.995, length.out = 120)
dens_post <- t(apply(psub, 1, function(lam) {
  k <- 0:K
  sapply(pgrid2, function(p) sum(lam * dbeta(p, k + 1, K - k + 1)))
}))                                        # 3000 x 120
band <- data.frame(
  p = pgrid2,
  mean = colMeans(dens_post),
  lo = apply(dens_post, 2, quantile, 0.025),
  hi = apply(dens_post, 2, quantile, 0.975))
p20 <- ggplot(band, aes(x = p)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, fill = 'steelblue') +
  geom_line(aes(y = mean), color = 'steelblue', linewidth = 1) +
  geom_line(data = data.frame(p = pgrid, y = dens_bern), aes(x = p, y = y),
            color = 'firebrick', linetype = 2, linewidth = 0.8) +
  labs(x = 'p', y = 'Density',
       title = 'Posterior mean prior density of p with 95% pointwise credible band',
       subtitle = 'solid: posterior mean (DA), dashed: MLE (Bern-Bino)') +
  theme_bw()
save_png(p20, 'plot20.png', w = 7, h = 5)

# ---- plot21: posterior predictive check ----
pred <- t(apply(psub, 1, function(lam) {
  kk <- sample(0:K, n, replace = TRUE, prob = lam)
  pp <- rbeta(n, kk + 1, K - kk + 1)
  rbinom(n, N, pp)
}))                                        # 3000 x 58
prop_obs <- R / N
prop_pred <- pred / matrix(N, nrow = nrow(pred), ncol = n, byrow = TRUE)
brks <- seq(0, 1, length.out = 16)
h_obs <- hist(prop_obs, breaks = brks, plot = FALSE)
mid <- h_obs$mids
pred_cnt <- t(apply(prop_pred, 1, function(v) hist(v, breaks = brks, plot = FALSE)$counts))
ppc_df <- data.frame(
  mid = mid,
  obs = h_obs$counts,
  mean = colMeans(pred_cnt),
  lo = apply(pred_cnt, 2, quantile, 0.025),
  hi = apply(pred_cnt, 2, quantile, 0.975))
p21 <- ggplot(ppc_df, aes(x = mid)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = 'gray40') +
  geom_line(aes(y = mean), color = 'steelblue', linewidth = 1) +
  geom_point(aes(y = obs), color = 'firebrick', size = 2) +
  labs(x = 'Proportion of dead fetuses per litter', y = 'Number of litters',
       title = 'Posterior predictive check',
       subtitle = 'points: observed, line: posterior predictive mean, band: 95% interval') +
  theme_bw()
save_png(p21, 'plot21.png', w = 7, h = 5)

# Bayesian p-value based on Pearson discrepancy
T_rep <- sapply(seq_len(nrow(psub)), function(j) {
  lam <- psub[j, ]
  k <- 0:K
  mu1 <- sum(lam * (k + 1) / (K + 2)); mu2 <- sum(lam * (k + 1) * (k + 2) / ((K + 2) * (K + 3)))
  pearson_stat(N, pred[j, ], mu1, mu2)
})
T_obs <- pearson_stat(N, R, moments('bern')['mu1'], moments('bern')['mu2'])
bay_p <- mean(T_rep >= T_obs)
cat(sprintf('Bayesian p-value (Pearson discrepancy) = %.4f\n', bay_p))

# ============================================================================
# 6. Summary output for the paper
# ============================================================================
cat('\n========== SUMMARY FOR SECTION 5 ==========\n')
cat(sprintf('modal litter size N = %d\n', modal_N))
cat(sprintf('Beta-Binomial: alpha = %.4f, beta = %.4f\n', alpha_bb, beta_bb))
cat(sprintf('Logit-Normal: mu = %.4f, sigma = %.4f\n', mu_ln, sigma_ln))
cat(sprintf('Kumaraswamy: a = %.4f, b = %.4f\n', a_km, b_km))
cat(sprintf('Bern-Bino: K = %d, lambda = (%s)\n', K_bic,
            paste(round(lambda_hat, 4), collapse = ', ')))
print(comp)
print(bay_tab)
cat(sprintf('Bayesian p-value = %.4f\n', bay_p))
cat('All outputs written to:', out_dir, 'and', fig_dir, '\n')
