suppressMessages({library(ltm); library(mirt)})
source('compare_models.R')
run <- function(nm, N, R, B = 999) {
  cat('\n=====', nm, '=====\n')
  tab <- compare_one(N, R, label=nm, Kmax=30, B=B, seed=12345)
  print(tab[, c('Model','K','npar','logLik','AIC','BIC','L1','MSE','pval')], row.names=FALSE)
  sub <- tab[!grepl('^bern_', tab$Model), ]
  bb <- tab[tab$Model=='bern_ll', ]; bc <- tab[tab$Model=='bern_bic', ]
  cat(sprintf('Bern(bestK=%d) wins logLik vs best parametric: %.3f vs %.3f (delta %.3f); BIC-selected K=%d wins BIC: %.1f vs %.1f\n',
      bb$K, bb$logLik, max(sub$logLik), bb$logLik - max(sub$logLik), bc$K, bc$BIC, min(sub$BIC)))
}
d <- ltm::Mobility; s <- rowSums(d)
run('Mobility', rep(8, length(s)), s)
dd <- mirt::deAyala; ds <- rowSums(dd[,1:5]); rr <- rep(ds, dd[,'Frequency'])
run('deAyala', rep(5, length(rr)), rr)
da <- mirt::ASVAB; items <- da[,1:4]; sc <- rowSums(items); rr2 <- NULL
for (g in 1:4) rr2 <- c(rr2, rep(sc, da[,4+g]))
run('ASVAB', rep(4, length(rr2)), rr2)
