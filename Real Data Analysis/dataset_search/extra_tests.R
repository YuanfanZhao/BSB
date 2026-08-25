suppressMessages({library(mirt); library(aod)})
source('compare_models.R')
# ASVAB: 4 binary items x 4 subgroups with counts
d <- mirt::ASVAB
items <- d[, 1:4]
cnt <- d[, 5:8]
sc <- rowSums(items)
rr <- NULL
for (g in 1:4) rr <- c(rr, rep(sc, cnt[, g]))
cat('ASVAB: n =', length(rr), ' m = 4\n')
tab <- compare_one(rep(4, length(rr)), rr, label='ASVAB', Kmax=20, B=999, seed=12345)
print(tab[, c('Model','K','npar','logLik','AIC','BIC','L1','MSE','pval')], row.names=FALSE)
cat('\n')

# orob1 definitive
d2 <- aod::orob1
tab2 <- compare_one(d2$n, d2$y, label='orob1', Kmax=20, B=999, seed=12345)
cat('--- orob1 (n=16, m varies) ---\n')
print(tab2[, c('Model','K','npar','logLik','AIC','BIC','L1','MSE','pval')], row.names=FALSE)
cat('\n')

# dja definitive
d3 <- aod::dja
tab3 <- compare_one(d3$n, d3$y, label='dja', Kmax=20, B=999, seed=12345)
cat('--- dja (n=75, m varies) ---\n')
print(tab3[, c('Model','K','npar','logLik','AIC','BIC','L1','MSE','pval')], row.names=FALSE)
