suppressMessages(library(ltm))
source('compare_models.R')
d <- ltm::Mobility; s <- rowSums(d); counts <- tabulate(s+1, nbins=9); n <- length(s); m <- 8
kgrid <- data.frame(K = 1:30)
kgrid$ll  <- sapply(kgrid$K, function(k) bern_mm_tab(m, counts, k)$ll)
kgrid$AIC <- 2*(kgrid$K+1) - 2*kgrid$ll
kgrid$BIC <- (kgrid$K+1)*log(n) - 2*kgrid$ll
write.csv(kgrid, 'mobility_K_grid.csv', row.names = FALSE)
cat('Mobility K-grid (K, ll, AIC, BIC):\n'); print(round(kgrid, 2))
cat('\nmin BIC at K =', kgrid$K[which.min(kgrid$BIC)], '\n')
