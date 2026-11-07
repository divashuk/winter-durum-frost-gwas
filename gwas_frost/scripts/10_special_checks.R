# Этап 10.6. Специальные проверки.
#  (а) Fr-B2 (гены CBF) и Mor534 против Q5B.1/Q5B.2: положения, LD (r2) ведущих SNP с SNP у генов CBF Fr-B2 и с Mor534,
#      условный региональный скан 5B 515–556 Мб (EMMAX, как в 09_perm_ld.R: происхождение + PC по BIC + добавочная ковариата).
#  (б) Fr-A2 / VRN-A1 против пика 5A: LD ведущего SNP 5A:514681016 с маркерами Fr-A2, VRN-A1 и SNP у генов CBF;
#      условный скан 5A 466–600 Мб с каждым маркером отдельно и со всеми маркерами Fr-A2 вместе.
#  (в) Q1B.6: гены с доменами LRR / NB-ARC / киназы — только состав доменов, без вывода о NLR (делается в 10_lrr_Q1B6.py).
suppressPackageStartupMessages({library(data.table); library(rrBLUP)})
OUT <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8")
B <- fread(file.path(OUT, "tables/01_blues.csv")); S2T <- fread(file.path(OUT, "tables/09_s2_traits.csv"))
s1 <- readLines(file.path(OUT, "tables/02_samples_stream1.txt")); s2 <- readLines(file.path(OUT, "tables/02_samples_stream2.txt"))
bic <- function(s, t) { f <- file.path(OUT, sprintf("tables/09_bic_%s_%s.csv", s, t)); if (!file.exists(f)) return(0); b <- fread(f); as.integer(b[[1]][which.max(b[[2]])]) }
whiten <- function(K, delta) { e <- eigen(K + delta * diag(nrow(K)), symmetric = TRUE); e$vectors %*% (t(e$vectors) / sqrt(e$values)) }
G <- list(); for (s in c("stream1", "stream2")) { raw <- fread(file.path(WORK, "qc", s, "final.raw")); X <- as.matrix(raw[, -(1:6)]); rownames(X) <- raw$IID; rm(raw)
  K <- as.matrix(fread(file.path(WORK, "kin", sprintf("K_vanraden_%s.csv", s)))[, -1]); rownames(K) <- colnames(K)
  key <- sub(":[^:]*:[^:]*_[^_]*$", "", colnames(X)); G[[s]] <- list(X = X, K = K, chr = sub(":.*", "", key), pos = as.numeric(sub(".*:", "", key)), key = key) }
COV <- list(stream1 = fread(file.path(OUT, "tables/09_covariates_S1.csv")), stream2 = fread(file.path(OUT, "tables/09_covariates_S2.csv")))
SETS <- list(
  S1 = list(stream = "stream1", d = B[trait == "FT_all" & scale == "ang" & TP_MO %in% s1, .(TP_MO, y = BLUE)], npc = bic("S1", "FT_all_ang"), trait = "FT_all_ang"),
  RU = list(stream = "stream1", d = B[trait == "FT_m12" & scale == "ang" & TP_MO %in% s1, .(TP_MO, y = BLUE)], npc = bic("RU", "FT_m12_ang"), trait = "FT_m12_ang"),
  S2 = list(stream = "stream2", d = S2T[trait == "FT12j" & TP_MO %in% s2, .(TP_MO, y = rint)], npc = bic("S2", "FT12j_rint"), trait = "FT12j_rint"),
  IT = list(stream = "stream2", d = S2T[trait == "FT12j" & coll == "IT" & TP_MO %in% s2, .(TP_MO, y)], npc = bic("IT", "ITFT_m12_ang"), trait = "ITFT_m12_ang"))
# региональный EMMAX-скан с добавочными ковариатами extra (data.table TP_MO + столбцы)
scan <- function(S, chr, lo, hi, extra = NULL, label = "базовая") {
  g <- G[[S$stream]]; d <- S$d[TP_MO %in% rownames(g$X)]
  if (!is.null(extra)) { d <- merge(d, extra, by = "TP_MO"); d <- d[complete.cases(d)] }
  ids <- d$TP_MO; X <- g$X[ids, , drop = FALSE]; K <- g$K[ids, ids]
  o <- droplevels(factor(OR$origin[match(ids, OR$TP_MO)])); C <- if (nlevels(o) > 1) model.matrix(~ o) else matrix(1, length(ids), 1)
  if (S$npc > 0) C <- cbind(C, prcomp(X[, apply(X, 2, sd) > 0], rank. = S$npc)$x)
  if (!is.null(extra)) C <- cbind(C, as.matrix(d[, setdiff(names(extra), "TP_MO"), with = FALSE]))
  C <- C[, qr(C)$pivot[seq_len(qr(C)$rank)], drop = FALSE]
  ms <- mixed.solve(d$y, X = C, K = K); Wm <- whiten(K, ms$Ve / ms$Vu)
  j <- which(g$chr == chr & g$pos >= lo & g$pos <= hi); Xj <- X[, j, drop = FALSE]; keep <- apply(Xj, 2, sd) > 0; j <- j[keep]; Xj <- Xj[, keep, drop = FALSE]
  ys <- Wm %*% d$y; Cs <- Wm %*% C; Xs <- Wm %*% Xj; Q <- qr.Q(qr(Cs)); res <- function(M) M - Q %*% crossprod(Q, M)
  yr <- res(ys); Xr <- res(Xs); xx <- colSums(Xr^2); b <- drop(crossprod(Xr, yr)) / xx; df <- length(ids) - ncol(C) - 1
  rss <- sum(yr^2) - b^2 * xx; p <- pf(b^2 * xx / (rss / df), 1, df, lower.tail = FALSE)
  data.table(model = label, n = length(ids), chr = chr, pos = g$pos[j], snp = g$key[j], mlp = -log10(p)) }
win_min <- function(R, wins) rbindlist(lapply(names(wins), function(w) { x <- R[pos >= wins[[w]][1] & pos <= wins[[w]][2]]
  if (!nrow(x)) return(NULL); x[, .(window = w, min_p_mlog10 = round(max(mlp), 2), best_snp = snp[which.max(mlp)], n_snp = .N), by = .(model, n)] }))
snpcov <- function(stream, chr, pos, nm) { g <- G[[stream]]; j <- which(g$chr == chr)[which.min(abs(g$pos[g$chr == chr] - pos))]
  data.table(TP_MO = rownames(g$X), v = g$X[, j])[, setnames(.SD, "v", nm)] }
markcov <- function(stream, m) { x <- COV[[stream]][, c("TP_MO", m), with = FALSE]; setnames(x, gsub("[^A-Za-z0-9]", "_", names(x))); x }
nearest <- function(stream, chr, pos) { g <- G[[stream]]; k <- which(g$chr == chr); k[which.min(abs(g$pos[k] - pos))] }
r2tab <- function(stream, ids, leads, targets) { g <- G[[stream]]; X <- g$X[intersect(ids, rownames(g$X)), ]
  rbindlist(lapply(names(leads), function(l) { jl <- nearest(stream, leads[[l]][1], as.numeric(leads[[l]][2]))
    rbindlist(lapply(names(targets), function(t) { tg <- targets[[t]]
      v <- if (is.character(tg) && length(tg) == 1) { mm <- COV[[stream]][match(rownames(X), TP_MO), get(tg)]; mm } else X[, nearest(stream, tg[1], as.numeric(tg[2]))]
      tsnp <- if (is.character(tg) && length(tg) == 1) paste("маркер", tg) else g$key[nearest(stream, tg[1], as.numeric(tg[2]))]
      data.table(lead = l, lead_snp = g$key[jl], target = t, target_snp = tsnp, n = sum(complete.cases(X[, jl], v)),
                 r2 = round(suppressWarnings(cor(X[, jl], v, use = "complete.obs"))^2, 3)) })) })) }

## ---------- (а) 5B: Q5B.1, Fr-B2, Q5B.2, Mor534 ----------
WIN5B <- list("Q5B.1 520.97–522.74" = c(520973288, 522742783), "Fr-B2 CBF 523.0–524.2" = c(523000000, 524200000),
              "Q5B.2 525.93–539.50" = c(525925602, 539497626), "Mor534 542.7–551.3" = c(542670002, 551332110))
P5B <- data.table(object = c("Q5B.1 ведущий SNP (S2)", "гены CBF Fr-B2 (Svevo → Langdon): TRITD5Bv1G167560", "TRITD5Bv1G167610", "TRITD5Bv1G167690",
                             "CBF15 → Fr-B2 (T. monococcum, 05)", "CBF12 → Fr-B2 (T. monococcum, 05)", "Q5B.2 ведущий SNP (S2)", "Q5B.2 ведущий SNP (IT)", "Q5B.2 второй SNP", "Mor534 (праймеры, 9.4)", "Mor534 LD-SNP (S1, r2 = 0.92)"),
                  pos = c(521642333, 523133648, 523174845, 523428524, 524074424, 524095658, 535618213, 535723924, 532482598, 550225979, 551172879))
P5B[, `:=`(dist_to_Q5B1_Mb = round((pos - 521642333) / 1e6, 2), dist_to_Q5B2_Mb = round((pos - 535618213) / 1e6, 2), dist_to_Mor534_Mb = round((pos - 550225979) / 1e6, 2))]
fwrite(P5B, file.path(OUT, "tables/10_check_FrB2_Mor534_positions.csv"))
L5B <- rbindlist(lapply(c("S2", "IT", "S1"), function(s) { S <- SETS[[s]]
  r2tab(S$stream, S$d$TP_MO, list("Q5B.1 521.64" = c("5B", 521642333), "Q5B.2 535.62" = c("5B", 535618213), "Q5B.2 532.48" = c("5B", 532482598)),
        list("CBF Fr-B2 523.13" = c("5B", 523133648), "CBF Fr-B2 523.43" = c("5B", 523428524), "CBF15→Fr-B2 524.07" = c("5B", 524074424),
             "Mor534 LD-SNP 551.17" = c("5B", 551172879), "Mor534 (импутированный маркер)" = "Mor534"))[, set := s] }))
fwrite(L5B, file.path(OUT, "tables/10_check_FrB2_Mor534_LD.csv")); print(L5B)
C5B <- rbindlist(lapply(c("S2", "IT"), function(s) { S <- SETS[[s]]; st <- S$stream
  rr <- list(scan(S, "5B", 515e6, 556e6),
             scan(S, "5B", 515e6, 556e6, snpcov(st, "5B", 535618213, "Q5B2_lead"), "+ Q5B.2 5B:535618213"),
             scan(S, "5B", 515e6, 556e6, snpcov(st, "5B", 521642333, "Q5B1_lead"), "+ Q5B.1 5B:521642333"),
             scan(S, "5B", 515e6, 556e6, snpcov(st, "5B", 523133648, "FrB2_snp"), "+ SNP у CBF Fr-B2 (523.13)"),
             scan(S, "5B", 515e6, 556e6, markcov(st, "Mor534"), "+ Mor534"))
  win_min(rbindlist(rr), WIN5B)[, `:=`(set = s, trait = S$trait)] }))
fwrite(C5B, file.path(OUT, "tables/10_check_FrB2_Mor534_conditional.csv")); print(C5B)

## ---------- (б) 5A: пик против Fr-A2 и VRN-A1 ----------
WIN5A <- list("ведущий SNP 5A:514681016" = c(514681016, 514681016), "Q5A.1 узкий 513.1–532.6" = c(513143615, 532591502),
              "Fr-A2 CBF 533.0–535.2" = c(533000000, 535200000), "VRN-A1 ± 1 Мб" = c(597268462, 599273527), "Q5A.4 609.2–622.4" = c(609173502, 622369126))
P5A <- data.table(object = c("ведущий SNP S1/RU", "KASP S1298957", "KASP S2269949", "KASP S1862541", "гены CBF Fr-A2 (Svevo → Langdon), первый", "последний", "KASP S1077313", "S1051014 (9.4)", "VRN-A1 (ген)", "Q5A.4 ведущий SNP"),
                  pos = c(514681016, 513240733, 524288011, 530965192, 533109223, 534643215, 537231599, 549361380, 598268462, 615862573))
P5A[, dist_to_lead_Mb := round((pos - 514681016) / 1e6, 2)]; fwrite(P5A, file.path(OUT, "tables/10_check_FrA2_VRNA1_positions.csv"))
L5A <- rbindlist(lapply(c("S1", "RU", "S2"), function(s) { S <- SETS[[s]]
  r2tab(S$stream, S$d$TP_MO, list("5A 514.68" = c("5A", 514681016)),
        list("S1298957" = "S1298957", "S2269949" = "S2269949", "S1862541" = "S1862541", "S1077313" = "S1077313", "S1051014" = "S1051014", "Vrn1_new" = "Vrn1_new",
             "SNP у CBF 533.11" = c("5A", 533109223), "SNP у CBF 534.64" = c("5A", 534643215), "SNP у VRN-A1 598.27" = c("5A", 598268462)))[, set := s] }))
fwrite(L5A, file.path(OUT, "tables/10_check_FrA2_VRNA1_LD.csv")); print(L5A)
FR <- c("S2269949", "S1077313", "S1862541", "S1298957")
C5A <- rbindlist(lapply(c("S1", "RU"), function(s) { S <- SETS[[s]]; st <- S$stream
  rr <- c(list(scan(S, "5A", 466e6, 625e6)),
          lapply(c(FR, "S1051014", "Vrn1_new"), function(m) scan(S, "5A", 466e6, 625e6, markcov(st, m), paste("+", m))),
          list(scan(S, "5A", 466e6, 625e6, markcov(st, FR), "+ все 4 маркера Fr-A2"),
               scan(S, "5A", 466e6, 625e6, markcov(st, c(FR, "Vrn1_new")), "+ 4 маркера Fr-A2 + Vrn1_new")))
  win_min(rbindlist(rr), WIN5A)[, `:=`(set = s, trait = S$trait)] }))
fwrite(C5A, file.path(OUT, "tables/10_check_FrA2_VRNA1_conditional.csv")); print(C5A[window %in% c("ведущий SNP 5A:514681016", "Fr-A2 CBF 533.0–535.2", "VRN-A1 ± 1 Мб")])
