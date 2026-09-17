# ============================================================================
# 05_report_figures.R --- export data for the figures of the regression report
#   (corrected Bern-Bino regression: eta_ik = phi_k + gamma_k * z_i^T beta,
#    z WITHOUT an intercept; K chosen by 5-fold repeated CV, see 06 script)
# Run: Rscript 05_report_figures.R
# ============================================================================
source("00_regression_lib.R")
std <- function(v) as.numeric(scale(v))
dat <- read.csv("data/mathexam_data.csv", stringsAsFactors = TRUE)
m <- dat$m[1]; x <- dat$x; n <- length(x)
Z <- cbind(tests = std(dat$tests), genderM = as.numeric(dat$gender == "male"),
           study571 = as.numeric(dat$study == "571"), semester = std(dat$semester))
set.seed(12345); te <- sample(1:n, floor(0.3*n)); tr <- setdiff(1:n, te)
Z_tr <- Z[tr,,drop=FALSE]; Z_te <- Z[te,,drop=FALSE]; x_tr <- x[tr]; x_te <- x[te]

ksel <- read.csv("results/me2_k_selection.csv")
K <- ksel$K[which.max(ksel$CV_mean_logLik)]
cat("K chosen by CV (log-likelihood) =", K, "\n")
fit <- fit_bernreg(Z_tr, m, x_tr, K, start = NULL, nstart = 12, nindep = 60)
pr <- fit$pred(Z_tr)
phi <- pr$phi; gam <- pr$gam; beta <- pr$beta

# fitted Bernstein prior density pi(p | z) at representative covariate profiles
pg <- seq(0.001, 0.999, length.out = 241)
density_at <- function(u) {
  eta <- phi + gam * u; eta <- eta - max(eta); lam <- exp(eta); lam <- lam/sum(lam)
  as.numeric(lam %*% t(sapply(0:K, function(k) dbeta(pg, k+1, K-k+1))))
}
u_tr <- as.numeric(Z_tr %*% beta); qs <- quantile(u_tr, c(0.10, 0.50, 0.90))
prior <- data.frame(p = pg, q10 = density_at(qs[1]), q50 = density_at(qs[2]),
                    q90 = density_at(qs[3]))
lam_te <- fit$pred(Z_te)$lambda
prior$mean_test <- colMeans(lam_te %*% t(sapply(0:K, function(k) dbeta(pg, k+1, K-k+1))))
write.csv(prior, "results/me2_prior_profiles.csv", row.names = FALSE)

# calibration on the test set
fl <- fit_logistic(cbind(1, Z_tr), m, x_tr); bb <- fit_bbreg(cbind(1, Z_tr), m, x_tr)
pl <- fl$pred(cbind(1, Z_te)); pb <- bb$pred(cbind(1, Z_te))$mu
pbern <- as.numeric(fit$pred(Z_te)$lambda %*% ((0:K+1)/(K+2)))
grp <- cut(rank(pbern, ties.method="first"), breaks=10, labels=FALSE)
agg <- data.frame(bin = 1:10, obs = tapply(x_te/m, grp, mean),
                  logistic = tapply(pl, grp, mean), bbreg = tapply(pb, grp, mean),
                  bern = tapply(pbern, grp, mean))
write.csv(agg, "results/me2_calibration.csv", row.names = FALSE)
cat("Exported results/me2_prior_profiles.csv and results/me2_calibration.csv\n")
cat("train logLik =", round(fit$ll,4), " test logLik =",
    round(test_loglik("bernreg", fit, Z_te, m, x_te),4), "\n")
