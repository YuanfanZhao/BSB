suppressMessages({
  library(lme4); library(ltm); library(mirt); library(aod)
  library(dispmod); library(hglm.data); library(VGAM); library(MLGdata)
})
source('compare_models.R')

load_ds <- function(pkg, name) { e <- new.env(); data(list = name, package = pkg, envir = e); get(name, envir = e) }

get_ds <- function(name) {
  switch(name,
    cbpp   = { d <- load_ds('lme4','cbpp'); list(N = d$size, R = d$incidence, note = 'CBPP 56 herd-period binomial counts') },
    lirat  = { d <- read.csv('../\u5927\u9f20\u81f4\u7578\u6570\u636e/lirat.csv'); list(N = d$N, R = d$R, note = 'low-iron rat teratology 58 (baseline)') },
    lsat   = { d <- load_ds('ltm','LSAT'); s <- rowSums(d); list(N = rep(5, length(s)), R = s, note = 'LSAT 1000 examinees m=5 (IRT)') },
    lsat6  = { d <- load_ds('mirt','LSAT6'); s <- rowSums(d[, 1:5]); rr <- rep(s, d$Freq); list(N = rep(5, length(rr)), R = rr, note = 'LSAT6 1000 examinees m=5') },
    lsat7  = { d <- load_ds('mirt','LSAT7'); s <- rowSums(d[, 1:5]); rr <- rep(s, d$freq); list(N = rep(5, length(rr)), R = rr, note = 'LSAT7 1000 examinees m=5') },
    orob1  = { d <- load_ds('aod','orob1'); list(N = d$n, R = d$y, note = 'Orobanche cernua 16 plates') },
    orob2  = { d <- load_ds('aod','orob2'); list(N = d$n, R = d$y, note = 'Crowder Orobanche 21 plates') },
    seeds  = { d <- load_ds('hglm.data','seeds'); list(N = d$n, R = d$r, note = 'Crowder seeds 21 plates') },
    toxop  = { d <- load_ds('VGAM','toxop'); list(N = d$ssize, R = d$positive, note = 'toxoplasmosis 34 cities (NPMLE classic)') },
    pneumo = { d <- load_ds('VGAM','pneumo'); tot <- d$normal + d$mild + d$severe; list(N = tot, R = d$mild + d$severe, note = 'pneumoconiosis 8 groups') },
    prats  = { d <- load_ds('VGAM','prats'); list(N = d$litter.size, R = d$alive, note = 'Weil rats 32 litters survival') },
    mice   = { d <- load_ds('aod','mice'); list(N = d$n, R = d$y, note = 'aod::mice 20 litters') },
    sal    = { d <- load_ds('hglm.data','salamander');
               tab <- tapply(d$Mate, d$Female, function(z) c(R = sum(z), N = length(z)))
               mat <- do.call(rbind, tab); list(N = mat[, 'N'], R = mat[, 'R'], note = 'salamander grouped by female') },
    infant = { d <- load_ds('MLGdata','Infant'); R0 <- ifelse(d$survival == 'No', d$Freq, 0);
               list(N = d$Freq, R = R0, note = 'Infant survival 16 cells (huge m)') },
    beetles= { d <- load_ds('MLGdata','Beetles'); list(N = d$num, R = d$uccisi, note = 'Beetles mortality 8 doses') },
    bioassay={ d <- load_ds('MLGdata','Bioassay'); list(N = d$den, R = d$y, note = 'Bioassay 10 doses') },
    fbeetle= { d <- load_ds('VGAM','fbeetle'); list(N = d$n, R = d$r, note = 'flour beetle mortality (VGAM)') },
    NULL
  )
}

names_all <- c('cbpp','lirat','lsat','lsat6','lsat7','orob1','orob2','seeds','toxop','pneumo',
               'prats','mice','sal','infant','beetles','bioassay','fbeetle')

phi <- function(N, R) {
  p <- sum(R) / sum(N); n <- length(N)
  sum((R - N * p)^2 / (N * p * (1 - p))) / (n - 1)
}

cat(sprintf('%-9s %5s %8s %7s  note\n', 'name', 'n', 'm-range', 'phi'))
datasets <- list()
for (nm in names_all) {
  d <- tryCatch(get_ds(nm), error = function(e) NULL)
  if (is.null(d)) { cat(sprintf('%-9s ERROR\n', nm)); next }
  ph <- tryCatch(phi(d$N, d$R), error = function(e) NA)
  cat(sprintf('%-9s %5d %8s %7.2f  %s\n', nm, length(d$R),
              paste(range(d$N), collapse = '-'), ph, d$note))
  datasets[[nm]] <- d
}
saveRDS(datasets, 'datasets_cache.rds')
cat('\nSaved', length(datasets), 'datasets\n')
