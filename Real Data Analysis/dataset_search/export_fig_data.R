# ============================================================================
# export_fig_data.R ?? ???? 5 ?????????????/?????
# ----------------------------------------------------------------------------
# ? 02_train_test_outsample.R ????? 70/30 ???seed=12345??
#   - ???????????Bern-Bino ?? 02 ??????????????? K
#     ?Mobility: K=25?deAyala: K=30???????????
#   - fig_*.csv ? observed ????????????? pmf ?????????
#     ????????????
#   - prior_*.csv ?????????????
#   - kgrid_*.csv ?????? K ?????BIC??? CV ???????/MSE??
# ??? make_section5_figures.py ??????
# ============================================================================
suppressMessages({library(ltm); library(mirt)})
source('compare_models.R')

export <- function(name, m, counts_all, Kb) {
  n <- sum(counts_all)
  # ---- ? 02 ????? 70/30 ?? ----
  set.seed(12345)
  scores <- rep(0:m, counts_all)                 # ???????
  test_idx  <- sample(1:n, size = floor(0.3 * n))
  counts_tr <- tabulate(scores[setdiff(1:n, test_idx)] + 1, nbins = m + 1)
  counts_te <- tabulate(scores[test_idx] + 1, nbins = m + 1)
  n1 <- sum(counts_tr); n2 <- sum(counts_te)

  # ---- ??????fit_model ?????? N?R ???????----
  N_tr <- rep(m, n1); R_tr <- rep(0:m, counts_tr)
  fits <- lapply(c('bin','bb','ln','km'), function(md) fit_model(md, N_tr, R_tr))
  names(fits) <- c('bin','bb','ln','km')
  fb <- bern_mm_tab(m, counts_tr, Kb)            # ?? K ? Bern-Bino

  # ---- fig???????? + ????? pmf ----
  x <- 0:m
  emp <- counts_te / n2
  pmf <- sapply(c('bin','bb','ln','km'), function(md) pmf_of(md, fits[[md]], m, x))
  pmf <- cbind(pmf, bern = bern_pmf_obs(rep(m, m + 1), x, fb$lambda))
  colnames(pmf) <- c('bin','bb','ln','km','bern')
  write.csv(data.frame(x = x, observed = emp, pmf),
            sprintf('fig_%s.csv', name), row.names = FALSE)

  # ---- prior??????????? ----
  pg <- seq(0.002, 0.998, length.out = 300)
  pr <- data.frame(
    p    = pg,
    bern = sapply(pg, function(pp) sum(fb$lambda * dbeta(pp, (0:Kb) + 1, Kb - (0:Kb) + 1))),
    bb   = dbeta(pg, fits$bb$extra$a, fits$bb$extra$b),
    ln   = f_ln(pg, fits$ln$extra$mu, fits$ln$extra$sigma),
    km   = f_km(pg, fits$km$extra$a, fits$km$extra$b))
  write.csv(pr, sprintf('prior_%s.csv', name), row.names = FALSE)

  cat(sprintf('%s: n1=%d n2=%d, ?? K=%d\n', name, n1, n2, Kb))
}

d1 <- ltm::Mobility; s1 <- rowSums(d1)
export('Mobility', 8, tabulate(s1 + 1, nbins = 9), Kb = 25)
d2 <- mirt::deAyala; ds <- rowSums(d2[, 1:5]); rr <- rep(ds, d2[, 'Frequency'])
export('deAyala', 5, tabulate(rr + 1, nbins = 6), Kb = 30)
cat('Exported fig_*.csv, prior_*.csv (train/test protocol)\n')
