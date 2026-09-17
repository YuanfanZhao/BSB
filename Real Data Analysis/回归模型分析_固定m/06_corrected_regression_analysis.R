# ============================================================================
# 06_corrected_regression_analysis.R
#   MathExam14W: corrected Bern-Bino regression (Section 2.5)
# ----------------------------------------------------------------------------
# MODEL (corrected):
#   eta_ik = phi_k + gamma_k * ( z_i^T beta ),   k = 0,1,...,K,
#   with phi_0 = 0 and gamma_0 = 1 as the baseline, and
#     * z_i  : covariates WITHOUT an intercept (the intercept lives in phi_k),
#     * beta : ONE regression vector shared by all Bernstein components,
#     * gamma_k : component-specific scale applied to the shared index,
#     * phi_k   : component-specific intercept.
#   Number of parameters = 2K + p', with p' = 4 covariates.
#
# The earlier version wrongly put an intercept column into z, so that
# phi_k and gamma_k*beta_0 were confounded; this script uses the corrected
# specification throughout.
#
# Protocol: 70/30 split (seed 12345); logistic, beta-binomial and Bern-Bino
# regression are fitted on the training set only; K is chosen inside the
# training set by 5-fold repeated CV with the validation log-likelihood;
# BIC and CV-MSE are recorded for comparison; all models are then evaluated
# on the held-out test set.
# Run: Rscript 06_corrected_regression_analysis.R
# ============================================================================
source("00_regression_lib.R")

std <- function(v) as.numeric(scale(v))
dat <- read.csv("data/mathexam_data.csv", stringsAsFactors = TRUE)
m <- dat$m[1]; x <- dat$x; n <- length(x)
Z <- cbind(tests    = std(dat$tests),
           genderM  = as.numeric(dat$gender == "male"),
           study571 = as.numeric(dat$study == "571"),
           semester = std(dat$semester))          # covariates only (NO intercept)
X <- cbind(1, Z)                                  # design matrix for logistic / BB
pz <- ncol(Z)

set.seed(12345)
te <- sample(1:n, size = floor(0.3 * n)); tr <- setdiff(1:n, te)
Z_tr <- Z[tr, , drop = FALSE]; Z_te <- Z[te, , drop = FALSE]
X_tr <- X[tr, , drop = FALSE]; X_te <- X[te, , drop = FALSE]
x_tr <- x[tr]; x_te <- x[te]
n1 <- length(x_tr); n2 <- length(x_te)
cat(sprintf("n = %d -> train n1 = %d, test n2 = %d; m = %d; covariates p' = %d\n",
            n, n1, n2, m, pz))

# ---- (1) logistic and beta-binomial regression on the training set --------
glm_fit  <- glm(cbind(x_tr, m - x_tr) ~ Z_tr, family = binomial)
fit_log  <- fit_logistic(X_tr, m, x_tr)          # same model, with pred closure
fit_bb   <- fit_bbreg(X_tr, m, x_tr)
ll_log   <- as.numeric(logLik(glm_fit))
H_bb     <- optimHess(fit_bb$theta, function(th) -bbreg_ll(th, X_tr, m, x_tr))
se_bb    <- sqrt(diag(solve(H_bb)))

# ---- (2) choice of K inside the training set -----------------------------
Kmax <- 15
bic_curve <- bernreg_fit_bic(Z_tr, m, x_tr, Kmax, nstart = 8, nindep = 8)$curve
K_bic <- as.integer(bic_curve$K[which.min(bic_curve$BIC)])

folds <- 5; reps <- 2
set.seed(6789)
cv_ll <- cv_mse <- matrix(NA, Kmax, folds * reps)
for (K in 1:Kmax) {
  col <- 0; start <- NULL
  for (r in 1:reps) {
    foldid <- sample(rep(1:folds, length.out = n1))
    for (f in 1:folds) {
      col <- col + 1
      va <- which(foldid == f); tr2 <- which(foldid != f)
      fK <- fit_bernreg(Z_tr[tr2, , drop = FALSE], m, x_tr[tr2], K,
                        start = start, nstart = 8, nindep = 8)
      start <- fK$theta
      cv_ll[K, col]  <- test_loglik("bernreg", fK, Z_tr[va, , drop = FALSE], m, x_tr[va], K)
      cv_mse[K, col] <- test_mse("bernreg", fK, Z_tr[va, , drop = FALSE], m, x_tr[va], K)
    }
  }
}
ksel <- data.frame(K = 1:Kmax, train_ll = bic_curve$train_ll, BIC = bic_curve$BIC,
                   CV_mean_logLik = rowMeans(cv_ll), CV_mean_MSE = rowMeans(cv_mse))
write.csv(ksel, "results/me2_k_selection.csv", row.names = FALSE)
K_cvll  <- ksel$K[which.max(ksel$CV_mean_logLik)]
K_cvmse <- ksel$K[which.min(ksel$CV_mean_MSE)]
cat(sprintf("K selected: BIC = %d ; CV(logLik) = %d ; CV(MSE) = %d\n", K_bic, K_cvll, K_cvmse))

# ---- (3) final Bern-Bino fit on the training set -------------------------
K_final <- K_cvll
fit_bern <- fit_bernreg(Z_tr, m, x_tr, K_final, start = NULL, nstart = 12, nindep = 60)
pr_tr <- fit_bern$pred(Z_tr)
phi <- pr_tr$phi; gam <- pr_tr$gam; beta <- pr_tr$beta
cat(sprintf("Bern-Bino K=%d: npar=%d, train logLik=%.4f\n", K_final, fit_bern$npar, fit_bern$ll))

# standard errors (may be unavailable: the Hessian is singular when some
# Bernstein components receive zero weight)
A_tab <- bern_coef_table(m, K_final)
H_b <- tryCatch(optimHess(fit_bern$theta,
        function(th) -bernreg_ll_grad(th, Z_tr, m, x_tr, A_tab)$ll), error = function(e) NULL)
se_b <- if (!is.null(H_b)) tryCatch(sqrt(diag(solve(H_b))), error = function(e) rep(NA, length(fit_bern$theta))) else rep(NA, length(fit_bern$theta))
ev_b <- if (!is.null(H_b)) eigen(H_b, symmetric = TRUE, only.values = TRUE)$values else NULL

# ---- (4) parameter tables -------------------------------------------------
log_tab <- data.frame(term = rownames(summary(glm_fit)$coefficients),
                      estimate = coef(glm_fit), se = summary(glm_fit)$coefficients[, 2],
                      z = summary(glm_fit)$coefficients[, 3],
                      p_value = summary(glm_fit)$coefficients[, 4], row.names = NULL)
write.csv(log_tab, "results/me2_params_logistic.csv", row.names = FALSE)

bb_tab <- data.frame(term = c(names(coef(glm_fit)), "logit_rho"),
                     estimate = fit_bb$theta, se = se_bb, row.names = NULL)
bb_tab <- rbind(bb_tab, data.frame(term = "rho", estimate = plogis(fit_bb$theta[length(fit_bb$theta)]),
                                   se = NA, row.names = NULL))
write.csv(bb_tab, "results/me2_params_bbreg.csv", row.names = FALSE)

bern_tab <- data.frame(
  term = c(paste0("phi_", 0:K_final), paste0("gamma_", 0:K_final),
           paste0("beta_", colnames(Z))),
  estimate = c(phi, gam, beta),
  se = c(0, se_b[1:K_final], 0, se_b[(K_final+1):(2*K_final)],
         se_b[(2*K_final+1):(2*K_final+pz)]))
bern_tab$se[1] <- NA; bern_tab$se[K_final + 2] <- NA   # phi_0, gamma_0 fixed
write.csv(bern_tab, "results/me2_params_bernreg.csv", row.names = FALSE)

# ---- (5) identification diagnostics --------------------------------------
lam_tr <- pr_tr$lambda
ident <- data.frame(
  component = 0:K_final,
  mean_weight = colMeans(lam_tr),
  max_weight = apply(lam_tr, 2, max))
write.csv(ident, "results/me2_identification.csv", row.names = FALSE)
scale_chk <- data.frame(c = c(0.25, 0.5, 1, 2, 5, 10),
                        logLik = sapply(c(0.25, 0.5, 1, 2, 5, 10), function(cc) {
                          if (cc == 1) return(fit_bern$ll)
                          bernreg_ll_grad(c(phi[-1], gam[-1]/cc, cc*beta), Z_tr, m, x_tr, A_tab)$ll
                        }))
write.csv(scale_chk, "results/me2_scale_check.csv", row.names = FALSE)

# ---- (6) training and test criteria --------------------------------------
npar_log <- length(coef(glm_fit)); npar_bb <- fit_bb$npar; npar_b <- fit_bern$npar
train_tab <- data.frame(
  Model = c("Logistic", "Beta-Binomial regression", "Bern-Bino regression"),
  npar = c(npar_log, npar_bb, npar_b), K = c(NA, NA, K_final),
  train_logLik = c(ll_log, fit_bb$ll, fit_bern$ll))
train_tab$train_AIC <- 2*train_tab$npar - 2*train_tab$train_logLik
train_tab$train_BIC <- train_tab$npar*log(n1) - 2*train_tab$train_logLik
write.csv(train_tab, "results/me2_train_metrics.csv", row.names = FALSE)

ev <- function(md, fit, D, K = NULL) c(
  test_logLik = test_loglik(md, fit, D, m, x_te, K),
  test_MAE    = test_mae(md, fit, D, m, x_te, K),
  test_MSE    = test_mse(md, fit, D, m, x_te, K))
test_tab <- data.frame(
  Model = c("Logistic", "Beta-Binomial regression", "Bern-Bino regression"),
  npar = c(npar_log, npar_bb, npar_b), K = c(NA, NA, K_final),
  rbind(ev("logistic", fit_log, X_te), ev("bbreg", fit_bb, X_te),
        ev("bernreg", fit_bern, Z_te)))
test_tab$test_AIC <- 2*test_tab$npar - 2*test_tab$test_logLik
test_tab$test_BIC <- test_tab$npar*log(n2) - 2*test_tab$test_logLik
write.csv(test_tab, "results/me2_test_metrics.csv", row.names = FALSE)


# small helper for printing
round_df <- function(d, dig = 4) {
  num <- sapply(d, is.numeric); d[num] <- lapply(d[num], function(z) round(z, dig)); d
}

# ---- print ---------------------------------------------------------------
cat("\n===== parameter estimates (training set) =====\n")
cat("\n-- logistic --\n");           print(round_df(log_tab), row.names = FALSE)
cat("\n-- beta-binomial --\n");      print(round_df(bb_tab), row.names = FALSE)
cat("\n-- Bern-Bino (K =", K_final, ") --\n"); print(round_df(bern_tab), row.names = FALSE)
cat("\n===== identification diagnostics =====\n")
print(ident, row.names = FALSE)
cat("\nscale check (beta -> c*beta, gamma_k -> gamma_k/c):\n"); print(scale_chk, row.names = FALSE)
if (!is.null(ev_b)) cat("\nHessian: min eigenvalue =", format(min(ev_b), digits=4),
                        " max =", format(max(ev_b), digits=4), "\n")
cat("\n===== training criteria =====\n"); print(round_df(train_tab), row.names = FALSE)
cat("\n===== test criteria =====\n");     print(round_df(test_tab), row.names = FALSE)
cat("\nSaved results/me2_*.csv\n")
