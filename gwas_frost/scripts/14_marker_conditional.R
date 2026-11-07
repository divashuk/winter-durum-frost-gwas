# Этап 14.2 (дополнение). Какой маркер статьи снимает локус: региональный EMMAX (как 10_special_checks.R) — схема A + один маркер за раз и схема B целиком.
#  RU (S1, FT_all, ang и rint): 1B 650.2–660.8 Мб; 4B 28.2–29.2.  IT (ITFT_m12, ang и rint): 7B 145.8–349.7, 2B 767.1–769.5, 2B 775.4–776.7, 1A 596.7–598.1, 5B 502.3–504.8.
suppressPackageStartupMessages({library(data.table); library(rrBLUP)})
OUT <- "/work/project/gwas_frost"; WORK <- "/work/scratch"
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8"); B <- fread(file.path(OUT, "tables/01_blues.csv"))
rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
whiten <- function(K, delta) { e <- eigen(K + delta * diag(nrow(K)), symmetric = TRUE); e$vectors %*% (t(e$vectors) / sqrt(e$values)) }
G <- list(); for (s in c("stream1", "stream2")) { raw <- fread(file.path(WORK, "qc", s, "final.raw")); X <- as.matrix(raw[, -(1:6)]); rownames(X) <- raw$IID
  K <- as.matrix(fread(file.path(WORK, "kin", sprintf("K_vanraden_%s.csv", s)))[, -1]); rownames(K) <- colnames(K)
  key <- sub(":[^:]*:[^:]*_[^_]*$", "", colnames(X)); G[[s]] <- list(X = X, K = K, chr = sub(":.*", "", key), pos = as.numeric(sub(".*:", "", key)), key = key); rm(raw) }
scan <- function(st, ids, y, origin, extra, ch, lo, hi) { g <- G[[st]]; X <- g$X[ids, ]; K <- g$K[ids, ids]
  C <- matrix(1, length(ids), 1)
  if (origin) { o <- droplevels(factor(OR$origin[match(ids, OR$TP_MO)])); if (nlevels(o) > 1) C <- model.matrix(~ o) }
  if (!is.null(extra)) { E <- as.matrix(extra); E <- E[, apply(E, 2, function(v) sd(v, na.rm = TRUE)) > 0, drop = FALSE]; C <- cbind(C, E) }
  C <- C[, qr(C)$pivot[seq_len(qr(C)$rank)], drop = FALSE]
  ms <- mixed.solve(y, X = C, K = K); W <- whiten(K, ms$Ve / ms$Vu)
  j <- which(g$chr == ch & g$pos >= lo & g$pos <= hi); Xj <- X[, j, drop = FALSE]; v <- apply(Xj, 2, sd) > 0; j <- j[v]; Xj <- Xj[, v, drop = FALSE]
  Q <- qr.Q(qr(W %*% C)); res <- function(M) M - Q %*% crossprod(Q, M); yr <- res(W %*% y); Xr <- res(W %*% Xj); xx <- colSums(Xr^2)
  b <- drop(crossprod(Xr, yr)) / xx; df <- length(ids) - ncol(C) - 1; p <- pf(b^2 * xx / ((sum(yr^2) - b^2 * xx) / df), 1, df, lower.tail = FALSE)
  data.table(best = g$key[j][which.min(p)], mlp = round(-log10(min(p)), 2)) }
RUN <- list(
  list(set = "RU (S1)", st = "stream1", cov = "09_covariates_S1.csv", origin = TRUE, ids = readLines(file.path(OUT, "tables/02_samples_stream1.txt")), tr = "FT_all",
       loci = list(c("1B", 650196761, 660809546), c("4B", 28212352, 29212352))),
  list(set = "IT", st = "stream2", cov = "09_covariates_S2.csv", origin = FALSE, ids = OR[collection == "italy", TP_MO], tr = "ITFT_m12",
       loci = list(c("7B", 145807057, 349683836), c("2B", 767142000, 769510000), c("2B", 775380000, 776720000), c("1A", 596700000, 598100000), c("5B", 502300000, 504800000))))
out <- list()
for (r in RUN) { COV <- fread(file.path(OUT, "tables", r$cov)); d <- B[trait == r$tr & scale == "ang" & TP_MO %in% r$ids & TP_MO %in% rownames(G[[r$st]]$X), .(TP_MO, y = BLUE)]
  mk <- setdiff(names(COV), "TP_MO"); Mm <- as.matrix(COV[match(d$TP_MO, TP_MO), ..mk])
  for (sc in c("ang", "rint")) { y <- if (sc == "ang") d$y else rint(d$y)
    for (l in r$loci) { ch <- l[1]; lo <- as.numeric(l[2]); hi <- as.numeric(l[3]); lab <- sprintf("%s %.1f–%.1f", ch, lo / 1e6, hi / 1e6)
      base <- scan(r$st, d$TP_MO, y, r$origin, NULL, ch, lo, hi)
      rows <- list(data.table(set = r$set, scale = sc, locus = lab, covariate = "схема A", base[, .(best, mlp)]))
      for (m in mk) if (sd(Mm[, m], na.rm = TRUE) > 0) rows[[length(rows) + 1]] <- data.table(set = r$set, scale = sc, locus = lab, covariate = paste("A +", m), scan(r$st, d$TP_MO, y, r$origin, Mm[, m, drop = FALSE], ch, lo, hi))
      rows[[length(rows) + 1]] <- data.table(set = r$set, scale = sc, locus = lab, covariate = "схема B (все маркеры)", scan(r$st, d$TP_MO, y, r$origin, Mm, ch, lo, hi))
      rows[[length(rows) + 1]] <- data.table(set = r$set, scale = sc, locus = lab, covariate = "схема Bn (без Mor)", scan(r$st, d$TP_MO, y, r$origin, Mm[, !grepl("^Mor", mk), drop = FALSE], ch, lo, hi))
      x <- rbindlist(rows); x[, drop_vs_A := round(x$mlp[1] - mlp, 2)]; out[[length(out) + 1]] <- x } } }
R <- rbindlist(out); fwrite(R, file.path(OUT, "tables/14_marker_conditional.csv"))
print(R[covariate %in% c("схема A", "схема B (все маркеры)", "схема Bn (без Mor)") | drop_vs_A >= 0.8], nrows = 300)
