# ============================================================================
# 02_train_test_regression.R ?? ?????????/????
# ----------------------------------------------------------------------------
# ?????"??/??"????????????
#   ? ??? 70/30 ????? / ??????????
#   ? ?????????????
#        Logistic??????glm?
#        Beta-Binomial ????? logit ?? + ???? rho?
#        Bern-Bino ?????? 2.5 ??K ?????????????
#           ?? A?BIC??? B?????? 5 ? x 3 ? CV?????????
#      ???Logit-Normal ???????????????
#   ? ???????????test logLik ? test MAE?
#      ?? BIC ??? CV ?? K ??????????????????? K?
# ??? data/ ??????????
# ============================================================================
source('00_regression_lib.R')

make_design <- function(dat) {
  if ('TypeF' %in% names(dat)) {
    X <- cbind(1, TypeF_W = as.numeric(dat$TypeF == 'W'),
                  TypeM_W = as.numeric(dat$TypeM == 'W'))
  } else if ('trisk' %in% names(dat)) {
    X <- cbind(1, group_TREAT = as.numeric(dat$group == 'TREAT'), trisk = dat$trisk)
  } else if ('seed' %in% names(dat)) {
    X <- cbind(1, seed_O73 = as.numeric(dat$seed == 'O73'),
                  root_CUK = as.numeric(dat$root == 'CUCUMBER'))
  } else if ('period' %in% names(dat)) {
    X <- cbind(1, period2 = as.numeric(dat$period == 2),
                  period3 = as.numeric(dat$period == 3),
                  period4 = as.numeric(dat$period == 4))
  } else if ('hb' %in% names(dat)) {
    X <- cbind(1, grp2 = as.numeric(dat$grp == 2), grp3 = as.numeric(dat$grp == 3),
                  grp4 = as.numeric(dat$grp == 4), hb = dat$hb)
  } else stop('unknown dataset')
  X
}

run_dataset <- function(file, label, Kmax = 8) {
  dat <- read.csv(file)
  X <- make_design(dat)
  m <- dat$m; x <- dat$x
  n <- length(m)

  set.seed(12345)                       # ? 70/30 ????????
  test_idx  <- sample(1:n, size = floor(0.3 * n))
  train_idx <- setdiff(1:n, test_idx)
  X_tr <- X[train_idx, , drop = FALSE]; m_tr <- m[train_idx]; x_tr <- x[train_idx]
  X_te <- X[test_idx,  , drop = FALSE]; m_te <- m[test_idx];  x_te <- x[test_idx]

  # ? ?????
  fit_log  <- fit_logistic(X_tr, m_tr, x_tr)
  fit_bb   <- fit_bbreg(X_tr, m_tr, x_tr)
  bern_bic <- bernreg_fit_bic(X_tr, m_tr, x_tr, Kmax)
  bern_cv  <- bernreg_fit_cv(X_tr, m_tr, x_tr, Kmax)
  cat('  [', label, '] Bern-Bino ?? K?BIC =', bern_bic$K,
      '???CV =', bern_cv$K, '\n')

  # ? ??????logistic / BB ?? / Bern-Bino(BIC) / Bern-Bino(CV)?
  ev <- function(md, fit) c(
    test_logLik = test_loglik(md, fit, X_te, m_te, x_te),
    test_MAE    = test_mae(md, fit, X_te, m_te, x_te))
  tab <- rbind(
    c(Model = 'Logistic', npar = fit_log$npar, K = NA, ev('logistic', fit_log)),
    c(Model = 'Beta-Binomial-reg', npar = fit_bb$npar, K = NA, ev('bbreg', fit_bb)),
    c(Model = 'Bern-Bino-reg(BIC)', npar = bern_bic$fit$npar, K = bern_bic$K, ev('bernreg', bern_bic$fit)),
    c(Model = 'Bern-Bino-reg(CV)', npar = bern_cv$fit$npar, K = bern_cv$K, ev('bernreg', bern_cv$fit)))
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  tab$npar <- as.integer(tab$npar); tab$K <- as.integer(tab$K)
  tab$test_logLik <- as.numeric(tab$test_logLik); tab$test_MAE <- as.numeric(tab$test_MAE)
  tab$label <- label
  tab
}

cat('========== ?????????/???? ==========\n\n')
res <- list()
for (f in c('data/salamander_data.csv', 'data/dja_data.csv', 'data/orob2_data.csv',
            'data/cbpp_data.csv', 'data/lirat_data.csv')) {
  lab <- sub('_data.csv', '', basename(f))
  cat('=====', lab, '=====\n')
  t <- tryCatch(run_dataset(f, lab), error = function(e) { cat('  ERROR:', conditionMessage(e), '\n'); NULL })
  if (!is.null(t)) { print(t, row.names = FALSE); res[[lab]] <- t; cat('\n') }
}
out <- do.call(rbind, res)
write.csv(out, 'results/regression_comparison.csv', row.names = FALSE)
cat('?????? results/regression_comparison.csv\n')
