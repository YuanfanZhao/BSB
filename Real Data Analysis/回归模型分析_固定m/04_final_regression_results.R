# NOTE: superseded by 06_corrected_regression_analysis.R, which uses the
# corrected Bern-Bino specification (intercept carried by phi_k, not by
# beta) and a multi-start search adequate for this multimodal likelihood.
# ============================================================================
# 04_final_regression_results.R --- MathExam14W regression: final tables
# ----------------------------------------------------------------------------
# Produces, for the thesis, the following outputs for the MathExam14W data:
#   (1) training/test samples and the 70/30 split
#   (2) parameter estimates on the TRAINING set for
#         - logistic regression
#         - beta-binomial regression
#         - Bern-Bino regression (Section 2.5)
#   (3) training-set fit criteria: logLik, AIC, BIC
#   (4) test-set validation: test logLik, test AIC, test BIC, test MSE, test MAE
#   (5) the K-selection curves inside the training set, for BOTH cross-validation
#       criteria (validation log-likelihood and validation MSE) and for BIC
#
# Design matrix (p = 5): intercept, tests (standardised), genderM, study571,
#                        semester (standardised)
# Run: Rscript 04_final_regression_results.R
# ============================================================================
source("00_regression_lib.R")

std <- function(v) as.numeric(scale(v))
make_design <- function(d) cbind(1,
    tests    = std(d$tests),
    genderM  = as.numeric(d$gender == "male"),
    study571 = as.numeric(d$study == "571"),
    semester = std(d$semester))

dat <- read.csv("data/mathexam_data.csv", stringsAsFactors = TRUE)
X <- make_design(dat); m <- dat$m[1]; x <- dat$x; n <- length(x)
p <- ncol(X)

set.seed(12345)
test_idx  <- sample(1:n, size = floor(0.3 * n))
train_idx <- setdiff(1:n, test_idx)
X_tr <- X[train_idx, , drop = FALSE]; x_tr <- x[train_idx]
X_te <- X[test_idx,  , drop = FALSE]; x_te <- x[test_idx]
n1 <- length(x_tr); n2 <- length(x_te)
cat(sprintf("n = %d -> train n1 = %d (%.1f%%), test n2 = %d (%.1f%%); m = %d, p = %d\n",
            n, n1, 100*n1/n, n2, 100*n2/n, m, p))

# ---- samples -------------------------------------------------------------
samp <- data.frame(x = 0:m,
                   train_n = tabulate(x_tr + 1, nbins = m + 1),
                   test_n  = tabulate(x_te + 1, nbins = m + 1))
samp$train_pct <- 100 * samp$train_n / n1
samp$test_pct  <- 100 * samp$test_n / n2
write.csv(samp, "results/me_samples.csv", row.names = FALSE)

# ---- (1) logistic regression --------------------------------------------
glm_fit <- glm(cbind(x_tr, m - x_tr) ~ X_tr[, -1, drop = FALSE], family = binomial)
fit_log_eval <- fit_logistic(X_tr, m, x_tr)   # same model, with a predict closure
sm <- summary(glm_fit)$coefficients
log_tab <- data.frame(term = rownames(sm), estimate = sm[, 1], se = sm[, 2],
                      z = sm[, 3], p_value = sm[, 4], row.names = NULL)
write.csv(log_tab, "results/me_params_logistic.csv", row.names = FALSE)
ll_log <- as.numeric(logLik(glm_fit))

# ---- (2) beta-binomial regression ---------------------------------------
fit_bb <- fit_bbreg(X_tr, m, x_tr)
th_bb  <- fit_bb$theta
H_bb   <- optimHess(th_bb, function(th) -bbreg_ll(th, X_tr, m, x_tr))
se_bb  <- sqrt(diag(solve(H_bb)))
bb_tab <- data.frame(term = c("beta0", "beta_tests", "beta_genderM",
                              "beta_study571", "beta_semester", "logit_rho"),
                     estimate = th_bb, se = se_bb, row.names = NULL)
bb_tab <- rbind(bb_tab, data.frame(term = "rho", estimate = plogis(th_bb[p+1]),
                                   se = NA, row.names = NULL))
write.csv(bb_tab, "results/me_params_bbreg.csv", row.names = FALSE)

# ---- (3) Bern-Bino: K choice inside the training set --------------------
Kmax <- 10
# BIC curve (training set only)
bic_curve <- bernreg_fit_bic(X_tr, m, x_tr, Kmax, nstart = 12)$curve
K_bic <- as.integer(bic_curve$K[which.min(bic_curve$BIC)])

# 5-fold x 2 repeated CV inside the training set, with BOTH criteria
folds <- 5; reps <- 2; K_grid <- 1:Kmax
set.seed(6789)
cv_ll <- cv_mse <- matrix(NA, Kmax, folds * reps)
for (K in K_grid) {
  col <- 0; start <- NULL
  for (r in 1:reps) {
    foldid <- sample(rep(1:folds, length.out = n1))
    for (f in 1:folds) {
      col <- col + 1
      va <- which(foldid == f); tr <- which(foldid != f)
      fitK <- fit_bernreg(X_tr[tr, , drop = FALSE], m, x_tr[tr], K,
                          start = start, nstart = 5)
      start <- fitK$theta
      cv_ll[K, col]  <- test_loglik("bernreg", fitK, X_tr[va, , drop = FALSE], m, x_tr[va], K)
      cv_mse[K, col] <- test_mse("bernreg", fitK, X_tr[va, , drop = FALSE], m, x_tr[va], K)
    }
  }
}
cv_tab <- data.frame(K = K_grid, train_ll = bic_curve$train_ll,
                     BIC = bic_curve$BIC, CV_mean_logLik = rowMeans(cv_ll),
                     CV_mean_MSE = rowMeans(cv_mse))
write.csv(cv_tab, "results/me_k_selection.csv", row.names = FALSE)
K_cv_ll  <- cv_tab$K[which.max(cv_tab$CV_mean_logLik)]
K_cv_mse <- cv_tab$K[which.min(cv_tab$CV_mean_MSE)]
cat(sprintf("K selected: BIC = %d ; CV(logLik) = %d ; CV(MSE) = %d\n",
            K_bic, K_cv_ll, K_cv_mse))

# final Bern-Bino fit on the whole training set (K from CV by log-likelihood)
K_final <- K_cv_ll
fit_bern <- fit_bernreg(X_tr, m, x_tr, K_final, start = NULL, nstart = 12)
th_b <- fit_bern$theta
H_b  <- tryCatch(optimHess(th_b, function(th)
                  -bernreg_ll_grad(th, X_tr, m, x_tr, bern_coef_table(m, K_final))$ll),
                 error = function(e) NULL)
se_b <- if (!is.null(H_b)) tryCatch(sqrt(diag(solve(H_b))), error = function(e) rep(NA, length(th_b))) else rep(NA, length(th_b))

phi  <- c(0, th_b[1:K_final])
gam  <- c(1, th_b[(K_final + 1):(2 * K_final)])
beta <- th_b[(2 * K_final + 1):(2 * K_final + p)]
se_beta <- se_b[(2 * K_final + 1):(2 * K_final + p)]

bern_tab <- data.frame(
  term = c(paste0("phi_", 0:K_final), paste0("gamma_", 0:K_final),
           c("beta0","beta_tests","beta_genderM","beta_study571","beta_semester")),
  estimate = c(phi, gam, beta),
  se = c(0, se_b[1:K_final], 0, se_b[(K_final+1):(2*K_final)], se_beta),
  group = c(rep("phi", K_final+1), rep("gamma", K_final+1), rep("beta", p)))
write.csv(bern_tab, "results/me_params_bernreg.csv", row.names = FALSE)
write.csv(data.frame(K = K_final, ll = fit_bern$ll, npar = fit_bern$npar),
          "results/me_bern_final_K.csv", row.names = FALSE)

# ---- (4) training and test criteria -------------------------------------
train_tab <- data.frame(
  Model = c("Logistic", "Beta-Binomial regression", "Bern-Bino regression"),
  npar  = c(length(coef(glm_fit)), fit_bb$npar, fit_bern$npar),
  K     = c(NA, NA, K_final),
  train_logLik = c(ll_log, fit_bb$ll, fit_bern$ll))
train_tab$train_AIC <- 2 * train_tab$npar - 2 * train_tab$train_logLik
train_tab$train_BIC <- train_tab$npar * log(n1) - 2 * train_tab$train_logLik

ev <- function(md, fit, K = NULL) c(
  test_logLik = test_loglik(md, fit, X_te, m, x_te, K),
  test_MAE    = test_mae(md, fit, X_te, m, x_te, K),
  test_MSE    = test_mse(md, fit, X_te, m, x_te, K))
test_tab <- data.frame(
  Model = c("Logistic", "Beta-Binomial regression", "Bern-Bino regression"),
  npar  = c(length(coef(glm_fit)), fit_bb$npar, fit_bern$npar),
  K     = c(NA, NA, K_final),
  rbind(ev("logistic", fit_log_eval), ev("bbreg", fit_bb),
        ev("bernreg", fit_bern)),
  row.names = NULL)
test_tab$test_AIC <- 2 * test_tab$npar - 2 * test_tab$test_logLik
test_tab$test_BIC <- test_tab$npar * log(n2) - 2 * test_tab$test_logLik
write.csv(train_tab, "results/me_train_metrics.csv", row.names = FALSE)
write.csv(test_tab,  "results/me_test_metrics.csv",  row.names = FALSE)

# ---- print --------------------------------------------------------------
cat("\n===== parameters: logistic =====\n"); print(log_tab, row.names = FALSE)
cat("\n===== parameters: beta-binomial =====\n"); print(bb_tab, row.names = FALSE)
cat("\n===== parameters: Bern-Bino (K =", K_final, ") =====\n")
print(bern_tab, row.names = FALSE)
cat("\n===== training criteria =====\n"); print(train_tab, row.names = FALSE)
cat("\n===== test criteria =====\n"); print(test_tab, row.names = FALSE)
cat("\nSaved CSVs in results/.\n")
