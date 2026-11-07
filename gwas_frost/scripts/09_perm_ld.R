# Этап 9.6. Пермутационный порог и межхромосомное LD.
#  Порог: смешанная модель EMMAX-типа (компоненты дисперсии при нулевой модели, rrBLUP::mixed.solve; P3D),
#   ковариаты схемы A (происхождение + PC по BIC). 1000 перестановок по Фридману–Лейну: остатки «выбеленной» модели
#   с ковариатами переставляются и добавляются к подогнанным значениям; для каждой перестановки — максимум −log10 p
#   по всем SNP; порог = 95-й перцентиль. S1: FT_all_ang; S2: FTcol_rint (основной вариант S2, п. 9.2) и FTcol_ang.
#  LD: ведущие SNP 1B, 4A, 5B со всеми SNP пика 5A (Q5A.1, Langdon 466–546 Мб): r2 и r2_V (с учётом структуры:
#   генотипы выбелены матрицей родства V = K, Mangin et al. 2012) — в S1 и S2.
suppressPackageStartupMessages({library(data.table); library(rrBLUP)})
OUT <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
set.seed(20261025)
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8")
B <- fread(file.path(OUT, "tables/01_blues.csv")); S2T <- fread(file.path(OUT, "tables/09_s2_traits.csv"))
bic <- function(s, t) { f <- file.path(OUT, sprintf("tables/09_bic_%s_%s.csv", s, t)); if (!file.exists(f)) return(0); b <- fread(f); as.integer(b[[1]][which.max(b[[2]])]) }
whiten <- function(K, delta) { e <- eigen(K + delta * diag(nrow(K)), symmetric = TRUE); e$vectors %*% (t(e$vectors) / sqrt(e$values)) }
perm_thr <- function(stream, trait, yv, ids, npc, nperm = 1000) {
  raw <- fread(file.path(WORK, "qc", stream, "final.raw")); X <- as.matrix(raw[, -(1:6)]); rownames(X) <- raw$IID; rm(raw)
  X <- X[ids, ]; X <- X[, apply(X, 2, sd) > 0]
  K <- as.matrix(fread(file.path(WORK, "kin", sprintf("K_vanraden_%s.csv", stream)))[, -1]); rownames(K) <- colnames(K); K <- K[ids, ids]
  o <- factor(OR$origin[match(ids, OR$TP_MO)]); C <- if (nlevels(droplevels(o)) > 1) model.matrix(~ droplevels(o)) else matrix(1, length(ids), 1)
  if (npc > 0) C <- cbind(C, prcomp(X, rank. = npc)$x)
  ms <- mixed.solve(yv, X = C, K = K); delta <- ms$Ve / ms$Vu
  Wm <- whiten(K, delta); ys <- Wm %*% yv; Cs <- Wm %*% C; Xs <- Wm %*% X
  Q <- qr.Q(qr(Cs)); res <- function(M) M - Q %*% crossprod(Q, M)
  Xr <- res(Xs); xx <- colSums(Xr^2); df <- length(yv) - ncol(C) - 1
  stat <- function(Yr) { b <- crossprod(Xr, Yr) / xx
    rss <- matrix(colSums(Yr^2), nrow(b), ncol(b), byrow = TRUE) - b^2 * xx   # остаточная SS для каждого SNP (строка) и перестановки (столбец)
    t2 <- b^2 * xx / (rss / df); -log10(pf(t2, 1, df, lower.tail = FALSE)) }
  yr <- res(ys); obs <- stat(matrix(yr))
  fitted <- ys - yr; mx <- numeric(nperm)
  for (b0 in seq(1, nperm, by = 50)) { nb <- min(50, nperm - b0 + 1)
    Yp <- sapply(seq_len(nb), function(i) fitted + yr[sample(length(yr))]); Ypr <- res(Yp)
    mx[b0:(b0 + nb - 1)] <- apply(stat(Ypr), 2, max) }
  data.table(stream = stream, trait = trait, n = length(yv), M = ncol(X), nPC = npc, delta = round(delta, 3),
             perm_thr_mlog10p_95 = round(quantile(mx, .95), 2), perm_thr_mlog10p_99 = round(quantile(mx, .99), 2),
             bonferroni_mlog10p = round(-log10(0.05 / ncol(X)), 2), observed_max_mlog10p = round(max(obs), 2),
             observed_top_snp = colnames(X)[which.max(obs)], n_snp_above_perm95 = sum(obs > quantile(mx, .95)))
}
s1 <- readLines(file.path(OUT, "tables/02_samples_stream1.txt")); s2 <- readLines(file.path(OUT, "tables/02_samples_stream2.txt"))
d1 <- B[trait == "FT_all" & scale == "ang" & TP_MO %in% s1]
PT <- list(perm_thr("stream1", "FT_all_ang", d1$BLUE, d1$TP_MO, bic("S1", "FT_all_ang")))
for (v in c("rint", "ang")) { d2 <- S2T[trait == "FTcol" & TP_MO %in% s2]; yv <- if (v == "rint") d2$rint else d2$y
  PT[[length(PT) + 1]] <- perm_thr("stream2", paste0("FTcol_", v), yv, d2$TP_MO, bic("S2", paste0("FTcol_", v))) }
PT <- rbindlist(PT); print(PT); fwrite(PT, file.path(OUT, "tables/09_permutation_thresholds.csv"))

## ---------- межхромосомное LD с пиком 5A ----------
LEAD <- data.table(QTL = c("Q1B.6", "Q4A.5", "Q5B.2", "Q5B.2", "Q5B.1"), chr = c("1B", "4A", "5B", "5B", "5B"),
                   pos = c(653803203, 642497807, 532482598, 535618213, 521615359))
LD <- list()
for (S in list(list("S1", "stream1", s1), list("S2", "stream2", s2))) {
  raw <- fread(file.path(WORK, "qc", S[[2]], "final.raw")); X <- as.matrix(raw[, -(1:6)]); rownames(X) <- raw$IID; rm(raw)
  X <- X[S[[3]], ]; key <- sub(":[^:]*:[^:]*_[^_]*$", "", colnames(X)); kc <- sub(":.*", "", key); kp <- as.numeric(sub(".*:", "", key))
  K <- as.matrix(fread(file.path(WORK, "kin", sprintf("K_vanraden_%s.csv", S[[2]])))[, -1]); rownames(K) <- colnames(K); K <- K[S[[3]], S[[3]]]
  e <- eigen(K, symmetric = TRUE); keep <- e$values > 1e-6; Vm <- e$vectors[, keep] %*% (t(e$vectors[, keep]) / sqrt(e$values[keep]))   # K^(-1/2) (псевдообратная)
  pk <- which(kc == "5A" & kp >= 466e6 & kp <= 546e6); Xp <- X[, pk, drop = FALSE]; Xpv <- Vm %*% scale(Xp, scale = FALSE)
  for (k in seq_len(nrow(LEAD))) { j <- match(paste(LEAD$chr[k], LEAD$pos[k], sep = ":"), key); if (is.na(j)) next
    r2 <- suppressWarnings(cor(X[, j], Xp)^2); r2v <- suppressWarnings(cor(Vm %*% (X[, j] - mean(X[, j])), Xpv)^2)
    LD[[length(LD) + 1]] <- data.table(set = S[[1]], QTL = LEAD$QTL[k], lead = paste(LEAD$chr[k], LEAD$pos[k], sep = ":"), n_snp_5A_peak = length(pk),
      r2_max = round(max(r2, na.rm = TRUE), 3), r2_median = round(median(r2, na.rm = TRUE), 3), snp_max_r2 = key[pk][which.max(r2)],
      r2V_max = round(max(r2v, na.rm = TRUE), 3), r2V_median = round(median(r2v, na.rm = TRUE), 3)) }
}
LD <- rbindlist(LD); print(LD); fwrite(LD, file.path(OUT, "tables/09_interchrom_ld_5A.csv"))
