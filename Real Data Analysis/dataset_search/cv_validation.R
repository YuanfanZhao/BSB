###############################################################################
# cv_validation.R - out-of-sample (train/test) model comparison
#   Fit all models on the training fold, evaluate test log-likelihood (and,
#   for constant m, test L1/MSE) on the held-out fold. Bern-Bino's K is
#   selected on the training fold by BIC (no test information leakage).
###############################################################################
suppressMessages(library(ltm))
source('compare_models.R')

# fit all models on training data; return list with fitted objects + bern K_bic
fit_models_train <- function(N_tr, R_tr, Kmax = 30) {
  fits <- lapply(c('bin','bb','ln','km'), function(m) fit_model(m, N_tr, R_tr))
  names(fits) <- c('bin','bb','ln','km')
  # Bern-Bino: select K by BIC on the training data
  if (length(unique(N_tr)) == 1) {
    m <- N_tr[1]; counts <- tabulate(R_tr + 1, nbins = m + 1); n_tr <- length(R_tr)
    llK <- sapply(1:Kmax, function(k) bern_mm_tab(m, counts, k)$ll)
    BICk <- (1:Kmax + 1) * log(n_tr) - 2 * llK
    Kb <- which.min(BICk)
    lam <- bern_mm_tab(m, counts, Kb)$lambda
    fits[['bern']] <- list(ll = llK[Kb], npar = Kb + 1, extra = list(lambda = lam, K = Kb))
  } else {
    llK <- sapply(1:Kmax, function(k) bern_mm(N_tr, R_tr, k)$ll)
    BICk <- (1:Kmax + 1) * log(length(R_tr)) - 2 * llK
    Kb <- which.min(BICk)
    f <- bern_mm(N_tr, R_tr, Kb)
    fits[['bern']] <- list(ll = f$ll, npar = Kb + 1, extra = list(lambda = f$lambda, K = Kb))
  }
  list(fits = fits, Kb = Kb)
}

# test log-likelihood of a fitted model on held-out obs
test_ll <- function(model, fit, N_te, R_te) {
  if (length(unique(N_te)) == 1 && length(unique(N_te)) == 1 && length(unique(N_te)) == 1) {
    m <- N_te[1]
    counts <- tabulate(R_te + 1, nbins = m + 1)
    pmf <- pmf_of(model, fit, m, 0:m); pmf <- pmf / sum(pmf)
    sum(counts * log(pmf))
  } else {
    sum(log(pmf_of(model, fit, N_te, R_te)))
  }
}

# test L1 / MSE vs empirical test distribution (constant m only)
test_disc <- function(model, fit, N_te, R_te) {
  m <- N_te[1]
  counts <- tabulate(R_te + 1, nbins = m + 1); n_te <- length(R_te)
  emp <- counts / n_te
  pmf <- pmf_of(model, fit, m, 0:m); pmf <- pmf / sum(pmf)
  c(L1 = 0.5 * sum(abs(emp - pmf)), MSE = mean((emp - pmf)^2))
}

# main CV comparison
cv_compare <- function(N, R, label, Kmax = 30, folds = 5, reps = 3, seed = 12345) {
  n <- length(N)
  models <- c('bin','bb','ln','km','bern')
  set.seed(seed)
  per_fold <- list()
  for (r in 1:reps) {
    foldid <- sample(rep(1:folds, length.out = n))
    for (f in 1:folds) {
      te <- which(foldid == f); tr <- which(foldid != f)
      ft <- fit_models_train(N[tr], R[tr], Kmax)
      fits <- ft$fits
      ll <- sapply(models, function(m) test_ll(m, fits[[m]], N[te], R[te]))
      disc <- sapply(models, function(m) {
        if (length(unique(N)) == 1) test_disc(m, fits[[m]], N[te], R[te]) else c(L1 = NA, MSE = NA)
      })
      per_fold[[length(per_fold) + 1]] <- data.frame(rep = r, fold = f,
        t(as.data.frame(ll)), t(as.data.frame(disc['L1', ])), t(as.data.frame(disc['MSE', ])),
        Kb = ft$Kb)
    }
  }
  res <- do.call(rbind, per_fold)
  names(res)[3:7] <- paste0('ll_', models)
  names(res)[8:12] <- paste0('L1_', models)
  names(res)[13:17] <- paste0('MSE_', models)
  # totals
  tot <- colSums(res[, paste0('ll_', models)])
  out <- data.frame(Model = models, CV_logLik = tot, mean_fold_ll = colMeans(res[, paste0('ll_', models)]))
  if (length(unique(N)) == 1) {
    out$CV_L1 <- colMeans(res[, paste0('L1_', models)])
    out$CV_MSE <- colMeans(res[, paste0('MSE_', models)])
  }
  out$label <- label
  attr(out, 'K_bic_per_split') <- res$Kb
  out
}

# ---- run on key datasets ----
run_cv <- function(N, R, label, Kmax = 30, folds = 5, reps = 3) {
  cat('\n=====', label, '=====\n')
  tab <- cv_compare(N, R, label, Kmax = Kmax, folds = folds, reps = reps, seed = 12345)
  print(tab, row.names = FALSE)
  cat('K (BIC on train) per split:', paste(attr(tab, 'K_bic_per_split'), collapse = ','), '\n')
  sub <- tab[tab$Model != 'bern', ]
  b <- tab[tab$Model == 'bern', ]
  cat(sprintf('Bern-Bino wins out-of-sample logLik: %s (delta %.3f)\n',
      ifelse(b$CV_logLik > max(sub$CV_logLik), 'YES', 'no'), b$CV_logLik - max(sub$CV_logLik)))
  invisible(tab)
}

d <- ltm::Mobility; s <- rowSums(d)
mob <- run_cv(rep(8, length(s)), s, 'Mobility (n=8445, m=8)', Kmax = 30, folds = 5, reps = 3)

suppressMessages(library(mirt))
dd <- mirt::deAyala; ds <- rowSums(dd[, 1:5]); rr <- rep(ds, dd[, 'Frequency'])
dea <- run_cv(rep(5, length(rr)), rr, 'deAyala (n=19601, m=5)', Kmax = 30, folds = 5, reps = 3)

suppressMessages(library(aod))
d2 <- aod::orob1
run_cv(d2$n, d2$y, 'orob1 (n=16, small biomedical)', Kmax = 20, folds = 4, reps = 2)

suppressMessages(library(lme4))
d3 <- lme4::cbpp
run_cv(d3$size, d3$incidence, 'cbpp (n=56)', Kmax = 20, folds = 5, reps = 2)

suppressMessages(library(VGAM))
d4 <- VGAM::toxop
run_cv(d4$ssize, d4$positive, 'toxop (n=34)', Kmax = 20, folds = 5, reps = 2)

saveRDS(list(Mobility = mob, deAyala = dea), 'cv_results.rds')
cat('\nSaved cv_results.rds\n')
