suppressMessages({library(ltm); library(mirt); library(ggplot2); library(tidyr); library(dplyr)})
source('compare_models.R')
d <- ltm::Mobility; s <- rowSums(d); m <- 8
counts <- tabulate(s + 1, nbins = m + 1); n <- length(s)

# fitted pmfs (reuse fast-m internals by calling compare_one and refitting)
tab <- compare_one(rep(m, length(s)), s, label='Mobility', Kmax = 30, B = 999, seed = 12345)
K <- tab$K[tab$Model == 'bern_bic']
lam <- bern_mm_tab(m, counts, K)$lambda
x <- 0:m
A <- sapply(0:K, function(kk) exp(lbeta(kk + 1 + x, K - kk + 1 + m - x) - lbeta(kk + 1, K - kk + 1)))
pmf_bern <- exp(lchoose(m, x)) * as.numeric(A %*% lam)
p_hat <- sum(s)/(n*m)
pmf_bin <- dbinom(x, m, p_hat)
fb <- optim(c(0,0), function(th){a<-exp(th[1]);b<-exp(th[2]); -sum(counts*(lchoose(m,x)+lbeta(a+x,b+m-x)-lbeta(a,b)))}, method='L-BFGS-B', lower=c(-8,-8), upper=c(8,8))
a_bb <- exp(fb$par[1]); b_bb <- exp(fb$par[2])
pmf_bb <- exp(lchoose(m,x) + lbeta(a_bb+x, b_bb+m-x) - lbeta(a_bb,b_bb))
emp <- counts/n
pdf2 <- data.frame(x=x, Observed=emp, Binomial=pmf_bin, `Beta-Binomial`=pmf_bb, `Bern-Bino`=pmf_bern, check.names=FALSE)
long <- pdf2 %>% pivot_longer(-x, names_to='model', values_to='prob')
long$model <- factor(long$model, levels=c('Observed','Binomial','Beta-Binomial','Bern-Bino'))
p <- ggplot(long, aes(x=x, y=prob, color=model, shape=model)) + geom_line(linewidth=1) + geom_point(size=2.5) +
  labs(x='Number of permitted mobility activities (score)', y='Probability',
       title='Mobility data: observed vs fitted score distributions (n=8445, m=8)') +
  theme_bw() + theme(legend.position='bottom')
ggsave('mobility_fit.png', p, width=8, height=5, dpi=300, bg='white')
# Bernstein prior density pi(p)
pg <- seq(0.005, 0.995, length.out=300)
dens <- sapply(pg, function(pp) sum(lam * dbeta(pp, (0:K)+1, K-(0:K)+1)))
pdf('mobility_prior_density.pdf')
plot(pg, dens, type='l', col='steelblue', lwd=2, xlab='p (mobility level)', ylab='Density',
     main=paste0('Estimated Bernstein prior density of p (K=', K, ')'))
dev.off()
cat('Bern-Bino K=', K, ' lambda=', paste(round(lam,4), collapse=', '), '\n')
cat('Beta-Binomial a=', round(a_bb,3), ' b=', round(b_bb,3), '\n')
cat('score distribution:', paste(counts, collapse=','), '\n')

# deAyala secondary
dd <- mirt::deAyala
ds <- rowSums(dd[,1:5]); rr <- rep(ds, dd[,'Frequency'])
tab2 <- compare_one(rep(5, length(rr)), rr, label='deAyala', Kmax=30, B=999, seed=12345)
cat('\n--- deAyala (n=', length(rr), ', m=5) ---\n')
print(tab2[, c('Model','K','npar','logLik','AIC','BIC','L1','MSE','pval')], row.names=FALSE)
