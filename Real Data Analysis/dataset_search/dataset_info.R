suppressMessages({library(ltm); library(mirt); library(ggplot2); library(tidyr); library(dplyr)})
source('compare_models.R')

# ---------- Mobility ----------
d <- ltm::Mobility; s <- rowSums(d); m <- ncol(d); n <- length(s)
cat('=== Mobility ===\n')
cat('n =', n, ' m =', m, '\n')
cat('score distribution:\n'); print(table(s))
cat('sum of counts:', sum(table(s)), '\n')
cat('proportion at score 8:', mean(s == 8), '\n')

# ---------- deAyala ----------
dd <- mirt::deAyala
cat('\n=== deAyala ===\n')
cat('dim:', dim(dd), ' names:', paste(colnames(dd), collapse=', '), '\n')
ds <- rowSums(dd[, 1:5]); rr <- rep(ds, dd[, 'Frequency'])
cat('total examinees (sum Freq):', sum(dd[, 'Frequency']), ' m =', 5, '\n')
cat('score distribution:\n'); print(table(rr))
cat('sum of counts:', length(rr), '\n')
cat('proportion at score 5:', mean(rr == 5), '\n')

# ---------- deAyala fit figure ----------
counts <- tabulate(rr + 1, nbins = 6); x <- 0:5
tab <- compare_one(rep(5, length(rr)), rr, label='deAyala', Kmax=30, B=999, seed=12345)
K <- tab$K[tab$Model == 'bern_bic']
lam <- bern_mm_tab(5, counts, K)$lambda
A <- sapply(0:K, function(kk) exp(lbeta(kk+1+x, K-kk+1+5-x) - lbeta(kk+1, K-kk+1)))
pmf_bern <- exp(lchoose(5, x)) * as.numeric(A %*% lam)
p_hat <- sum(rr)/(length(rr)*5)
pmf_bin <- dbinom(x, 5, p_hat)
fb <- optim(c(0,0), function(th){a<-exp(th[1]);b<-exp(th[2]); -sum(counts*(lchoose(5,x)+lbeta(a+x,b+5-x)-lbeta(a,b)))}, method='L-BFGS-B', lower=c(-8,-8), upper=c(8,8))
a_bb <- exp(fb$par[1]); b_bb <- exp(fb$par[2])
pmf_bb <- exp(lchoose(5,x) + lbeta(a_bb+x, b_bb+5-x) - lbeta(a_bb,b_bb))
emp <- counts/length(rr)
pdf2 <- data.frame(x=x, Observed=emp, Binomial=pmf_bin, `Beta-Binomial`=pmf_bb, `Bern-Bino`=pmf_bern, check.names=FALSE)
long <- pdf2 %>% pivot_longer(-x, names_to='model', values_to='prob')
long$model <- factor(long$model, levels=c('Observed','Binomial','Beta-Binomial','Bern-Bino'))
p <- ggplot(long, aes(x=x, y=prob, color=model, shape=model)) + geom_line(linewidth=1) + geom_point(size=2.5) +
  labs(x='Test score (number correct out of 5)', y='Probability',
       title='deAyala data: observed vs fitted score distributions (n=19601, m=5)') +
  theme_bw() + theme(legend.position='bottom')
ggsave('deayala_fit.png', p, width=8, height=5, dpi=300, bg='white')
cat('\nBern-Bino K=', K, ' lambda=', paste(round(lam,4), collapse=', '), '\n')
cat('Beta-Binomial a=', round(a_bb,3), ' b=', round(b_bb,3), '\n')

# deAyala prior density
pg <- seq(0.005, 0.995, length.out=300)
dens <- sapply(pg, function(pp) sum(lam * dbeta(pp, (0:K)+1, K-(0:K)+1)))
pdf('deayala_prior_density.pdf')
plot(pg, dens, type='l', col='steelblue', lwd=2, xlab='p (ability)', ylab='Density',
     main=paste0('deAyala: estimated Bernstein prior density of p (K=', K, ')'))
dev.off()

# Mobility prior density values (for report)
d2 <- ltm::Mobility; s2 <- rowSums(d2); counts2 <- tabulate(s2+1, nbins=9); K2 <- 15
lam2 <- bern_mm_tab(8, counts2, K2)$lambda
cat('\nMobility Bern-Bino K=', K2, ' lambda=', paste(round(lam2,4), collapse=', '), '\n')
