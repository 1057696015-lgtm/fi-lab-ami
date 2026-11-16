# Minimum observed-component thresholds 12-18 (Supplementary Tables S6-S7).
# Usage: Rscript threshold.R DATA_DIR RESULTS_DIR   (DATA_DIR holds M_all.csv / E_all.csv from threshold_export.py)
# At each threshold both models are refitted in MIMIC-IV and applied unchanged to eICU, with exported and
# creatinine-based SOFA; 1,000 paired patient bootstrap samples; paired DeLong test.
args <- commandArgs(TRUE); dd <- args[1]; od <- args[2]; set.seed(20261101)
M <- read.csv(file.path(dd, 'M_all.csv'), check.names = FALSE); E <- read.csv(file.path(dd, 'E_all.csv'), check.names = FALSE)
cv <- setdiff(names(M), c('sofa', 'sofa_cr', 'fi_any', 'ncomp', 'death30', 'fol'))
auc <- function(y, p) { n <- sum(y); (sum(rank(p)[y == 1]) - n * (n + 1) / 2) / (n * (length(y) - n)) }
fit <- function(x, y) { f <- glm.fit(x, y, family = binomial()); stopifnot(f$converged, all(is.finite(f$coefficients))); f }
met <- function(y, l0, l1) { p0 <- plogis(l0); p1 <- plogis(l1)
  c(AUC0 = auc(y, p0), AUC1 = auc(y, p1), dAUC = auc(y, p1) - auc(y, p0), Brier0 = mean((p0 - y)^2), Brier1 = mean((p1 - y)^2),
    OE1 = sum(y) / sum(p1), slope1 = coef(glm(y ~ l1, family = binomial()))[2]) }
res <- list(); k <- 0
for (t in 12:18) for (s in c('sofa', 'sofa_cr')) {
  m <- M[M$ncomp >= t & M$fol > 1 & !is.na(M[[s]]), ]; e <- E[E$ncomp >= t & E$fol > 1 & !is.na(E[[s]]), ]
  x0 <- cbind(1, as.matrix(m[, c(cv, s)])); x1 <- cbind(x0, fi01 = m$fi_any * 10)
  z0 <- cbind(1, as.matrix(e[, c(cv, s)])); z1 <- cbind(z0, fi01 = e$fi_any * 10)
  f0 <- fit(x0, m$death30); f1 <- fit(x1, m$death30); l0 <- drop(z0 %*% f0$coefficients); l1 <- drop(z1 %*% f1$coefficients); y <- e$death30
  point <- met(y, l0, l1); boot <- t(replicate(1000, { ix <- sample.int(length(y), replace = TRUE); met(y[ix], l0[ix], l1[ix]) }))
  test <- pROC::roc.test(pROC::roc(y, l0, quiet = TRUE, direction = '<'), pROC::roc(y, l1, quiet = TRUE, direction = '<'), paired = TRUE, method = 'delong')
  k <- k + 1
  res[[k]] <- data.frame(threshold = t, sofa = s, nM = nrow(m), eventsM = sum(m$death30), nE = nrow(e), eventsE = sum(y), metric = names(point),
                         estimate = point, lower = apply(boot, 2, quantile, .025), upper = apply(boot, 2, quantile, .975), DeLongP = test$p.value)
  cat('Done', t, s, nrow(m), nrow(e), 'delta', point['dAUC'], '\n')
}
write.csv(do.call(rbind, res), file.path(od, 'threshold_performance.csv'), row.names = FALSE)
r <- list()
for (co in c('M', 'E')) {
  d <- if (co == 'M') M else E; d <- d[d$fol > 1 & !is.na(d$sofa), ]
  g <- ifelse(d$ncomp >= 18, 'ge18', ifelse(d$ncomp >= 12, '12to17', 'lt12'))
  for (v in c('age', 'sexM', 'charlson', 'sofa', 'fi_any', 'ncomp', 'death30')) for (gg in c('ge18', '12to17', 'lt12')) {
    a <- d[g == gg, v]; b <- d[g == 'ge18', v]
    r[[length(r) + 1]] <- data.frame(cohort = co, group = gg, variable = v, n = sum(!is.na(a)), mean = mean(a, na.rm = TRUE), sd = sd(a, na.rm = TRUE),
                                     SMD_vs_ge18 = (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2))
  }
}
write.csv(do.call(rbind, r), file.path(od, 'included_excluded.csv'), row.names = FALSE)
