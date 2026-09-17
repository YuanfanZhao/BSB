# ============================================================================
# 03_final_results.R --- Women's Mobility: training/test results (final tables)
# ----------------------------------------------------------------------------
# Outputs required for the thesis:
#   (A) training and test samples (score frequency distributions)
#   (B) training-set fits: parameter estimates, logLik, AIC, BIC and the
#       chi-square goodness-of-fit statistic with its simulated p-value
#       (identical to the chi-square test used in the simulation study, Sec. 4)
#   (C) test-set validation: test logLik, test AIC, test BIC, test MSE, test L1
#
# Models: Binomial / Beta-Binomial / Logit-Normal-Binomial /
#         Kumaraswamy-Binomial / Bern-Bino (K = 25, chosen inside the training
#         set by 5-fold x 3 repeated CV with the validation likelihood).
# Note: the MM algorithm needs many iterations for K = 25; the default cap of
#       2000 was too small, so maxit = 2e5 with eps = 1e-12 is used here.
# Run: Rscript 03_final_results.R
# ============================================================================
source("00_functions.R")

scores <- read.csv("mobility_data.csv")$score
m <- max(scores); n <- length(scores)

# ---- (0) 70/30 split (same seed as script 02) ----------------------------
set.seed(12345)
test_idx  <- sample(1:n, size = floor(0.3 * n))
train_idx <- setdiff(1:n, test_idx)
counts_tr <- tabulate(scores[train_idx] + 1, nbins = m + 1)
counts_te <- tabulate(scores[test_idx]  + 1, nbins = m + 1)
n1 <- length(train_idx); n2 <- length(test_idx)
cat(sprintf("n = %d  ->  training n1 = %d (%.1f%%), test n2 = %d (%.1f%%)\n",
            n, n1, 100 * n1 / n, n2, 100 * n2 / n))

samp <- data.frame(x = 0:m, train_n = counts_tr, train_pct = 100 * counts_tr / n1,
                   test_n = counts_te, test_pct = 100 * counts_te / n2)
write.csv(samp, "results/samples_train_test.csv", row.names = FALSE)
write.csv(data.frame(score = scores[train_idx]), "results/scores_train.csv", row.names = FALSE)
write.csv(data.frame(score = scores[test_idx]),  "results/scores_test.csv",  row.names = FALSE)

# ---- (1) fit the five models on the training set -------------------------
K_final <- 25
bern_f  <- bern_mm(m, counts_tr, K_final, maxit = 200000, eps = 1e-12)

fits <- list(
  bin  = fit_model("bin", m, counts_tr),
  bb   = fit_model("bb",  m, counts_tr),
  ln   = fit_model("ln",  m, counts_tr),
  km   = fit_model("km",  m, counts_tr),
  bern = list(ll = bern_f$ll, npar = K_final + 1,
              pmf = bern_pmf(m, counts_tr, bern_f$lambda),
              extra = list(lambda = bern_f$lambda, K = K_final, iter = bern_f$iter)))
cat(sprintf("Bern-Bino K=%d: MM converged in %d iterations, logLik = %.4f\n",
            K_final, bern_f$iter, bern_f$ll))

# ---- (2) goodness of fit: classical Pearson chi-square + simulated p -----
gof_classic <- function(counts, pmf, B = 999, seed = 12345) {
  obs <- as.numeric(counts); p <- pmf / sum(pmf)
  E   <- sum(obs) * p
  X2  <- sum((obs - E)^2 / E)
  set.seed(seed)
  pv  <- suppressWarnings(chisq.test(x = obs, p = p, simulate.p.value = TRUE, B = B)$p.value)
  c(X2 = X2, p = pv)
}

# ---- (3) training-set results table --------------------------------------
mnames <- c(bin = "Binomial", bb = "Beta-Binomial", ln = "Logit-Normal-Binomial",
            km = "Kumaraswamy-Binomial", bern = "Bern-Bino")
pnames <- c(bin = "p", bb = "a, b", ln = "mu, sigma", km = "a, b",
            bern = paste0("lambda_0,...,lambda_", K_final, " (K = ", K_final, ")"))
pv <- function(md) {
  ex <- fits[[md]]$extra
  switch(md,
    bin  = sprintf("p = %.4f", ex$p),
    bb   = sprintf("a = %.4f, b = %.4f", ex$a, ex$b),
    ln   = sprintf("mu = %.4f, sigma = %.4f", ex$mu, ex$sigma),
    km   = sprintf("a = %.4f, b = %.4f", ex$a, ex$b),
    bern = paste(sprintf("%.5f", ex$lambda), collapse = ", "))
}

train_tab <- do.call(rbind, lapply(names(fits), function(md) {
  g  <- gof_classic(counts_tr, fits[[md]]$pmf)
  gm <- gof_p(m, counts_tr, fits[[md]]$pmf, B = 999, seed = 12345)   # moment-based
  k  <- fits[[md]]$npar
  data.frame(Model = mnames[md], Parameters = pnames[md], npar = k,
             Estimate = pv(md), logLik = fits[[md]]$ll,
             AIC = 2 * k - 2 * fits[[md]]$ll,
             BIC = k * log(n1) - 2 * fits[[md]]$ll,
             Pearson_X2 = g["X2"], GOF_p_sim = g["p"],
             Moment_X2 = pearson_stat(m, counts_tr, fits[[md]]$pmf),
             Moment_p = gm, stringsAsFactors = FALSE)
}))
write.csv(train_tab, "results/train_fit_results.csv", row.names = FALSE)

# ---- (4) test-set validation table ---------------------------------------
test_eval <- function(pmf) {
  p  <- pmf / sum(pmf)
  ll <- sum(ifelse(counts_te > 0, counts_te * log(p), 0))
  d  <- dist_emp(m, counts_te, p)
  c(test_logLik = ll, test_L1 = unname(d["L1"]), test_MSE = unname(d["MSE"]))
}
test_tab <- do.call(rbind, lapply(names(fits), function(md) {
  k <- fits[[md]]$npar; e <- test_eval(fits[[md]]$pmf)
  data.frame(Model = mnames[md], npar = k, test_logLik = e["test_logLik"],
             test_AIC = 2 * k - 2 * e["test_logLik"],
             test_BIC = k * log(n2) - 2 * e["test_logLik"],
             test_MSE = e["test_MSE"], test_L1 = e["test_L1"],
             stringsAsFactors = FALSE)
}))
write.csv(test_tab, "results/test_validation_results.csv", row.names = FALSE)

# ---- (5) print -----------------------------------------------------------
cat("\n===== (A) Training and test samples =====\n"); print(samp, row.names = FALSE)
cat("\n===== (B) Training-set fits =====\n")
print(train_tab[, c("Model","npar","logLik","AIC","BIC","Pearson_X2","GOF_p_sim",
                    "Moment_X2","Moment_p")], row.names = FALSE)
cat("\nParameter estimates (training set):\n")
for (md in names(fits)) cat(sprintf("  %-22s %s\n", mnames[md], pv(md)))
cat("\n===== (C) Test-set validation =====\n"); print(test_tab, row.names = FALSE)
cat("\nSaved: results/samples_train_test.csv, train_fit_results.csv, test_validation_results.csv\n")
