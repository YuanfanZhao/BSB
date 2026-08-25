suppressMessages(library(ltm))
source('compare_models.R')
d <- ltm::Mobility
s <- rowSums(d); m <- ncol(d)
cat('Mobility: n =', length(s), ' m =', m, '\n')
cat('Score distribution:\n'); print(table(s))
# save dataset
df <- data.frame(woman_id = seq_len(length(s)), score = s)
write.csv(df, 'mobility_scores.csv', row.names = FALSE)
write.csv(d, 'mobility_items.csv', row.names = FALSE)

tab <- compare_one(rep(m, length(s)), s, label = 'Mobility', Kmax = 30, B = 999, seed = 12345)
print(tab[, c('Model','K','npar','logLik','AIC','BIC','L1','MSE','pval')], row.names = FALSE)
write.csv(tab, 'mobility_full_comparison.csv', row.names = FALSE)
cat('\nBIC-selected K:', tab$K[tab$Model=='bern_bic'], ' best-logLik K:', tab$K[tab$Model=='bern_ll'], '\n')
