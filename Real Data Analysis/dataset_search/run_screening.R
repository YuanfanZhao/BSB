source('compare_models.R')
datasets <- readRDS('datasets_cache.rds')
cat('Screening', length(datasets), 'datasets with B=99, Kmax=15 ...\n\n')
res_all <- list()
for (nm in names(datasets)) {
  d <- datasets[[nm]]
  cat('=====', nm, '=====\n')
  tab <- tryCatch(compare_one(d$N, d$R, label = nm, Kmax = 15, B = 99, seed = 12345),
                  error = function(e) { cat('ERROR:', conditionMessage(e), '\n'); NULL })
  if (is.null(tab)) next
  res_all[[nm]] <- tab
  sub <- tab[!grepl('^bern_', tab$Model), ]
  bb <- tab[tab$Model == 'bern_ll', ]
  bc <- tab[tab$Model == 'bern_bic', ]
  cat(sprintf('  best parametric: %s logLik=%.3f (AIC=%.1f BIC=%.1f)\n',
      sub$Model[which.max(sub$logLik)], max(sub$logLik),
      sub$AIC[which.max(sub$logLik)], sub$BIC[which.max(sub$logLik)]))
  cat(sprintf('  Bern-Bino: bestK=%d logLik=%.3f (AIC=%.1f BIC=%.1f) | BIC K=%d logLik=%.3f\n',
      bb$K, bb$logLik, bb$AIC, bb$BIC, bc$K, bc$logLik))
  cat(sprintf('  Bern wins logLik: %s | wins L1: %s | wins MSE: %s | best GOF p: %s (%.3f), Bern p: %.3f\n',
      ifelse(bb$logLik > max(sub$logLik), 'YES', 'no'),
      ifelse(bb$L1 <= min(sub$L1) + 1e-9, 'YES', 'no'),
      ifelse(bb$MSE <= min(sub$MSE) + 1e-12, 'YES', 'no'),
      sub$Model[which.max(sub$pval)], max(sub$pval), bb$pval))
  cat('\n')
}
saveRDS(res_all, 'screening_results.rds')
cat('Saved screening_results.rds\n')
