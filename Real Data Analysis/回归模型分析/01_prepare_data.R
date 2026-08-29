# ============================================================================
# 01_prepare_data.R ?? ??????????????
# ----------------------------------------------------------------------------
# ?? data/ ?? CSV????????? (x, m) ?????
# ???????????????? X = cbind(1, ????)??????
# ============================================================================
dir.create('data', showWarnings = FALSE)

# ---- 1. salamander???????????????????----
# ????? 3 ? R ???? 3 ? W ??????? (??, ????) ???
#   x = ???????0..3??m = 3????????? TypeF????? TypeM?R/W??
e <- new.env(); data('salamander', package = 'hglm.data', envir = e); d <- e$salamander
g <- aggregate(Mate ~ Female + TypeM, data = d, FUN = function(z) c(x = sum(z), m = length(z)))
g$mx <- g$Mate[, 1]; g$mm <- g$Mate[, 2]     # aggregate 返回矩阵列，取出 x 与 m
finfo <- unique(d[, c('Female', 'TypeF')])
g <- merge(g, finfo, by = 'Female')
g <- data.frame(x = g$mx, m = g$mm, TypeF = g$TypeF, TypeM = g$TypeM)
write.csv(g, 'data/salamander_data.csv', row.names = FALSE)
cat('salamander: n =', nrow(g), ' m =', unique(g$m), '\n')

# ---- 2. dja???????75 ??----
e2 <- new.env(); data('dja', package = 'aod', envir = e2); d2 <- e2$dja
d2 <- data.frame(x = d2$y, m = d2$n, group = d2$group, trisk = d2$trisk)
write.csv(d2, 'data/dja_data.csv', row.names = FALSE)
cat('dja: n =', nrow(d2), ' m ?? =', paste(range(d2$m), collapse='-'), '\n')

# ---- 3. orob2?Crowder ???????21 ??----
e3 <- new.env(); data('orob2', package = 'aod', envir = e3); d3 <- e3$orob2
d3 <- data.frame(x = d3$y, m = d3$n, seed = d3$seed, root = d3$root)
write.csv(d3, 'data/orob2_data.csv', row.names = FALSE)
cat('orob2: n =', nrow(d3), '\n')

# ---- 4. cbpp??????????56 ? herd-period?----
e4 <- new.env(); data('cbpp', package = 'lme4', envir = e4); d4 <- e4$cbpp
d4 <- data.frame(x = d4$incidence, m = d4$size, period = d4$period)
write.csv(d4, 'data/cbpp_data.csv', row.names = FALSE)
cat('cbpp: n =', nrow(d4), ' m ?? =', paste(range(d4$m), collapse='-'), '\n')


# ---- 5. lirat（大鼠低铁致畸，58 窝；协变量 grp、hb）----
# 注意：lirat.csv 存放在 ../大鼠致畸数据/（中文路径，R 无法直接读取），
# 已用 Python 复制为 data/lirat_data.csv（列：x=R, m=N, grp, hb）。此处直接使用该副本。
if (file.exists('data/lirat_data.csv')) {
  d5 <- read.csv('data/lirat_data.csv')
  cat('lirat: n =', nrow(d5), ' m 范围 =', paste(range(d5$m), collapse='-'), '\n')
} else {
  cat('lirat: 请先用 Python 从 ..\u5927\u9f20\u81f4\u7578\u6570\u636e/lirat.csv 复制为 data/lirat_data.csv\n')
}
