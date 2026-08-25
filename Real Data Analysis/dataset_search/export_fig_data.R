suppressMessages({library(ltm); library(mirt)})
source('compare_models.R')

export <- function(name, m, counts, Kmax = 30) {
  n <- sum(counts); x <- 0:m
  emp <- counts / n
  # fit each model via the pipeline helpers (constant m)
  N <- rep(m, n); R <- rep(x, counts)
  fits <- lapply(c('bin','bb','ln','km'), function(mm) fit_model(mm, N, R))
  names(fits) <- c('bin','bb','ln','km')
  # Bern-Bino: BIC K
  llK <- sapply(1:Kmax, function(k) bern_mm_tab(m, counts, k)$ll)
  BICk <- (1:Kmax + 1) * log(n) - 2 * llK
  AICk <- 2 * (1:Kmax + 1) - 2 * llK
  Kb <- which.min(BICk)
  fb <- bern_mm_tab(m, counts, Kb)
  pmf <- sapply(c('bin','bb','ln','km'), function(mm) pmf_of(mm, fits[[mm]], m, x))
  pmf <- cbind(pmf, bern = bern_pmf_obs(rep(m, length(x)), x, fb$lambda))
  colnames(pmf) <- c('bin','bb','ln','km','bern')
  out <- data.frame(x = x, observed = emp, pmf)
  write.csv(out, sprintf('fig_%s.csv', name), row.names = FALSE)
  # prior densities
  pg <- seq(0.002, 0.998, length.out = 300)
  d_bb <- dbeta(pg, fits$bb$extra$a, fits$bb$extra$b)
  d_ln <- f_ln(pg, fits$ln$extra$mu, fits$ln$extra$sigma)
  d_km <- f_km(pg, fits$km$extra$a, fits$km$extra$b)
  d_bern <- sapply(pg, function(pp) sum(fb$lambda * dbeta(pp, (0:Kb) + 1, Kb - (0:Kb) + 1)))
  pr <- data.frame(p = pg, bern = d_bern, bb = d_bb, ln = d_ln, km = d_km)
  write.csv(pr, sprintf('prior_%s.csv', name), row.names = FALSE)
  # K grid
  kg <- data.frame(K = 1:Kmax, ll = llK, AIC = AICk, BIC = BICk)
  write.csv(kg, sprintf('kgrid_%s.csv', name), row.names = FALSE)
  cat(sprintf('%s: K_bic=%d, logLik: bin=%.2f bb=%.2f ln=%.2f km=%.2f bern=%.2f\n',
      name, Kb, fits$bin$ll, fits$bb$ll, fits$ln$ll, fits$km$ll, fb$ll))
  invisible(list(fits = fits, Kb = Kb, fb = fb))
}

d1 <- ltm::Mobility; s1 <- rowSums(d1)
c1 <- tabulate(s1 + 1, nbins = 9)
r1 <- export('mobility', 8, c1)
d2 <- mirt::deAyala; s2 <- rowSums(d2[, 1:5]); rr <- rep(s2, d2[, 'Frequency'])
c2 <- tabulate(rr + 1, nbins = 6)
r2 <- export('deayala', 5, c2)
cat('Exported fig_*.csv, prior_*.csv, kgrid_*.csv\n')
