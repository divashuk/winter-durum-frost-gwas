# Этап 14.1–14.2. Сводка локусов по наборам S2, IT, RU (S1), F2 и мета-анализам RU + F2, IT + RU + F2.
#  Источники (GEMMA, схема A; B и Bn — для «снимается ли маркерами»):
#   S2  — FT12j_ang, FT12j_rint (all_results_09_S2: A, B, Bn);
#   IT  — ITFT_m12_ang (09: A), ITFT_m12_rint (12.1 «в»: A), B и Bn ang/rint (раздел 14, заново);
#   RU  — S1 FT_all_ang (09: A, B; 14: Bn), FT_all_rint (14: A, B, Bn — заново);
#   F2  — скан Хейли–Нотта раздела 8 (логистическая модель, семья — фиксированный эффект, доза по HMM): «без Mor» = A, «с Mor» (Mor-прокси) = B.
#  Мета-анализ: взвешенный Стауффер по z (вес √n), z ориентирован на аллель ALT (идентификатор chr:pos:REF:ALT); только общие SNP; раздельно для шкал ang и rint
#   (F2 — один и тот же бинарный скан в обеих). METAL в окружении нет.
#  Пики: −log10 p >= 4 (максимум из ang и rint); локус = SNP той же хромосомы в пределах ±5 Мб от ведущего или с r² > 0.2 с ведущим (жадно, от сильнейшего).
#   r² — в панели набора (S2: поток 2, все S2; IT: поток 2, IT; RU и мета RU + F2: поток 1, S1; F2: генотипы F2; мета IT + RU + F2: поток 2, S2).
#  Общие локусы: интервалы разных наборов на одной хромосоме, перекрывающиеся или ближе 5 Мб, объединяются. Локус 1B (Q1B.6, 650.2–660.8 Мб) добавляется обязательно.
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"; WORK <- "/work/scratch"
rd_gemma <- function(f, set, scale, scheme, n) { r <- fread(f); data.table(set, scale, scheme, id = r$rs, chr = sub(":.*", "", r$rs), pos = r$ps, p = r$p_wald, beta = r$beta,
  ea = r$allele1, maf = pmin(r$af, 1 - r$af), n = n) }
R9 <- function(s) fread(file.path(WORK, "gwas09", sprintf("all_results_09_%s.csv.gz", s)))[model == "GEMMA"]
from9 <- function(x, set_, scale_, n_) x[, .(set = set_, scale = scale_, scheme, id = SNP, chr, pos, p, beta = effect, ea = effect_allele, maf = MAF, n = n_)]   # имена аргументов не совпадают со столбцами
A <- list()
s2 <- R9("S2"); A$S2a <- from9(s2[trait == "FT12j_ang"], "S2", "ang", 273); A$S2r <- from9(s2[trait == "FT12j_rint"], "S2", "rint", 273); rm(s2)
it <- R9("IT"); A$ITa <- from9(it[trait == "ITFT_m12_ang"], "IT", "ang", 182)
A$ITr <- rd_gemma(file.path(WORK, "gwas12/c_IT_rint/output/lmm.assoc.txt"), "IT", "rint", "A", 182)
for (sc in c("B", "Bn")) for (tr in c("ang", "rint")) A[[paste0("IT", tr, sc)]] <- rd_gemma(file.path(WORK, sprintf("gwas14/IT_%s_%s/output/lmm.assoc.txt", tr, sc)), "IT", tr, sc, 182)
ru <- R9("S1"); A$RUa <- from9(ru[trait == "FT_all_ang"], "RU (S1)", "ang", 144); rm(ru)
A$RUaBn <- rd_gemma(file.path(WORK, "gwas14/RU_ang_Bn/output/lmm.assoc.txt"), "RU (S1)", "ang", "Bn", 144)
for (sc in c("A", "B", "Bn")) A[[paste0("RUr", sc)]] <- rd_gemma(file.path(WORK, sprintf("gwas14/RU_rint_%s/output/lmm.assoc.txt", sc)), "RU (S1)", "rint", sc, 144)
f2 <- fread(file.path(WORK, "f2/08_f2_scan_hk.csv.gz"))
f2raw <- fread(file.path(WORK, "f2/f2.raw")); fm <- colMeans(as.matrix(f2raw[, -(1:6)]), na.rm = TRUE) / 2; names(fm) <- sub("_[^_]*$", "", names(f2raw)[-(1:6)])
f2[, `:=`(id = sub("_[^_]*$", "", snp), ea = sub(".*_", "", snp))]; f2[, maf := pmin(fm[id], 1 - fm[id])]
for (cv in c("без Mor", "с Mor")) for (sc in c("ang", "rint")) A[[paste0("F2", sc, cv)]] <- f2[covariates == cv, .(set = "F2", scale = sc, scheme = ifelse(cv == "без Mor", "A", "B"), id, chr, pos, p = p_add, beta = beta_add, ea, maf, n = 171)]
D <- rbindlist(A, use.names = TRUE)
D[, alt := sub("^[^:]+:[0-9]+:[^:]+:", "", id)]
D[, z := sign(beta) * qnorm(pmax(p, 1e-300) / 2, lower.tail = FALSE) * fifelse(ea == alt, 1, -1)]
## мета-анализ
meta <- function(sets, label) rbindlist(lapply(c("ang", "rint"), function(sc) {
  x <- D[scheme == "A" & scale == sc & set %in% sets]; common <- x[, .N, by = id][N == length(sets), id]; x <- x[id %in% common]
  x[, .(set = label, scale = sc, scheme = "A", chr = chr[1], pos = pos[1], z = sum(sqrt(n) * z) / sqrt(sum(n)), maf = mean(maf), n = sum(n),
        dirs = paste(sign(z)[order(set)], collapse = "")), by = id][, p := 2 * pnorm(-abs(z))][, `:=`(beta = z, ea = sub("^[^:]+:[0-9]+:[^:]+:", "", id))] }))
M1 <- meta(c("RU (S1)", "F2"), "мета RU + F2"); M2 <- meta(c("IT", "RU (S1)", "F2"), "мета IT + RU + F2")
cnt <- D[scheme == "A" & scale == "ang", .(n_snp = uniqueN(id)), by = set][match(c("S2", "IT", "RU (S1)", "F2"), set)]
CS <- rbind(data.table(meta = "RU + F2", common = uniqueN(M1$id)), data.table(meta = "IT + RU + F2", common = uniqueN(M2$id)))
for (k in seq_len(nrow(cnt))) CS[, (paste0("share_of_", cnt$set[k])) := round(common / cnt$n_snp[k], 3)]
fwrite(cnt, file.path(OUT, "tables/14_snp_counts_by_set.csv")); fwrite(CS, file.path(OUT, "tables/14_meta_common_snps.csv")); print(cnt); print(CS)
ALL <- rbind(D[, .(set, scale, scheme, id, chr, pos, p, z, maf)], M1[, .(set, scale, scheme, id, chr, pos, p, z, maf)], M2[, .(set, scale, scheme, id, chr, pos, p, z, maf)])
ALL[, mlp := -log10(p)]
fwrite(ALL[scheme == "A" & mlp >= 3], file.path(WORK, "gwas14", "all_A_mlp3.csv.gz"))
fwrite(rbind(M1, M2), file.path(WORK, "gwas14", "meta_results.csv.gz"))
## пороги
f2X <- as.matrix(f2raw[, -(1:6)]); colnames(f2X) <- names(fm)
simpleM <- function(X) { ch <- sub(":.*", "", colnames(X)); sum(sapply(unique(ch), function(c) { Y <- X[, ch == c, drop = FALSE]; Y <- Y[, apply(Y, 2, function(v) sd(v, na.rm = TRUE)) > 0, drop = FALSE]
  Y <- apply(Y, 2, function(v) { v[is.na(v)] <- mean(v, na.rm = TRUE); v }); ev <- eigen(cor(Y), symmetric = TRUE, only.values = TRUE)$values; ev <- pmax(ev, 0)
  which(cumsum(ev) / sum(ev) >= 0.995)[1] })) }
meff_f2 <- simpleM(f2X)
TH <- data.table(set = c("S2", "IT", "RU (S1)", "F2", "мета RU + F2", "мета IT + RU + F2"),
                 meff_n = c(3717, 3717, 1731, meff_f2, 1731, 3717),
                 M = c(cnt[set == "S2", n_snp], cnt[set == "IT", n_snp], cnt[set == "RU (S1)", n_snp], cnt[set == "F2", n_snp], uniqueN(M1$id), uniqueN(M2$id)))
TH[, `:=`(meff_thr = round(-log10(0.05 / meff_n), 2), bonf_thr = round(-log10(0.05 / M), 2))]
TH[, meff_source := c("поток 2 (раздел 3)", "поток 2 (раздел 3)", "поток 1 (раздел 3)", "simpleM по генотипам F2 (раздел 14)", "поток 1", "поток 2")]
fwrite(TH, file.path(OUT, "tables/14_thresholds.csv")); print(TH)
## панели генотипов для r²
ldpanel <- list(); for (s in c("stream1", "stream2")) { raw <- fread(file.path(WORK, "qc", s, "final.raw")); X <- as.matrix(raw[, -(1:6)]); rownames(X) <- raw$IID
  colnames(X) <- sub("_[^_]*$", "", colnames(X)); ldpanel[[s]] <- X; rm(raw) }
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8")
S2ids <- intersect(readLines(file.path(OUT, "tables/02_samples_stream2.txt")), rownames(ldpanel$stream2)); ITids <- intersect(OR[collection == "italy", TP_MO], rownames(ldpanel$stream2))
S1ids <- intersect(readLines(file.path(OUT, "tables/02_samples_stream1.txt")), rownames(ldpanel$stream1))
panel_of <- function(set) switch(set, "S2" = ldpanel$stream2[S2ids, ], "IT" = ldpanel$stream2[ITids, ], "RU (S1)" = ldpanel$stream1[S1ids, ], "F2" = f2X,
                                 "мета RU + F2" = ldpanel$stream1[S1ids, ], "мета IT + RU + F2" = ldpanel$stream2[S2ids, ])
## пики и локусы по наборам
loci <- list()
for (s in TH$set) {
  x <- ALL[set == s & scheme == "A", .SD[which.max(mlp)], by = id][mlp >= 4][order(-mlp)]
  P <- panel_of(s); used <- rep(FALSE, nrow(x)); k <- 0
  while (any(!used)) { i <- which(!used)[1]; k <- k + 1; ld <- x$id[i]
    cand <- which(!used & x$chr == x$chr[i]); r2 <- rep(NA_real_, length(cand))
    if (ld %in% colnames(P)) { ok <- x$id[cand] %in% colnames(P); if (any(ok)) r2[ok] <- suppressWarnings(cor(P[, ld], P[, x$id[cand][ok], drop = FALSE], use = "pairwise.complete.obs"))[1, ]^2 }
    mem <- cand[abs(x$pos[cand] - x$pos[i]) <= 5e6 | (!is.na(r2) & r2 > 0.2)]
    used[mem] <- TRUE
    loci[[length(loci) + 1]] <- data.table(set = s, local_id = k, chr = x$chr[i], start = min(x$pos[mem]), end = max(x$pos[mem]), lead = ld, lead_pos = x$pos[i],
                                           lead_mlp = round(x$mlp[i], 2), lead_scale = x$scale[i], lead_maf = round(x$maf[i], 3), n_snp_ge4 = length(mem)) } }
LOC <- rbindlist(loci)
## общие локусы
LOC <- rbind(LOC, data.table(set = "обязательный", local_id = 0L, chr = "1B", start = 650196761, end = 660809546, lead = "1B:653803203", lead_pos = 653803203, lead_mlp = NA, lead_scale = NA, lead_maf = NA, n_snp_ge4 = 0L))
setorder(LOC, chr, start)
LOC[, g := { cl <- integer(.N); cur <- 0; e <- -Inf; for (j in seq_len(.N)) { if (start[j] > e + 5e6) cur <- cur + 1; cl[j] <- cur; e <- max(e, end[j]) }; cl }, by = chr]
GL <- LOC[, .(start = min(start), end = max(end), sets = paste(unique(setdiff(set, "обязательный")), collapse = ", ")), by = .(chr, g)]
GL[, locus := sprintf("%s %.1f–%.1f", chr, start / 1e6, end / 1e6)]
LOC <- merge(LOC, GL[, .(chr, g, locus)], by = c("chr", "g"))
fwrite(LOC[order(chr, start)], file.path(OUT, "tables/14_loci_by_set_raw.csv"))
## длинная таблица: локус × набор
MK <- data.table(marker = c("S1298957", "S2269949", "S1862541", "S1077313", "S1051014", "Vrn1_new/VRN-A1", "Vrn-B1", "Mor534", "Mor587", "Mor866"),
                 chr = c("5A", "5A", "5A", "5A", "5A", "5A", "5B", "5B", "1B", "7B"), pos = c(513240733, 524288011, 530965192, 537231599, 549361380, 598268462, 605129460, 550225979, 54736990, 610783371))
best_in <- function(s, sch, ch, a, b) { x <- ALL[set == s & scheme == sch & chr == ch & pos >= a & pos <= b]; if (!nrow(x)) return(NULL); x[which.max(mlp)] }
long <- list()
for (j in seq_len(nrow(GL))) { gl <- GL[j]
  for (s in TH$set) { th <- TH[set == s]; xa <- best_in(s, "A", gl$chr, gl$start, gl$end)
    if (is.null(xa)) { long[[length(long) + 1]] <- data.table(locus = gl$locus, chr = gl$chr, set = s, note = "нет SNP в интервале"); next }
    xb <- best_in(s, "B", gl$chr, gl$start, gl$end); xbn <- best_in(s, "Bn", gl$chr, gl$start, gl$end)
    lc <- LOC[locus == gl$locus & set == s]
    mk_in <- MK[chr == gl$chr & pos >= gl$start - 5e6 & pos <= gl$end + 5e6, paste(marker, collapse = ", ")]
    long[[length(long) + 1]] <- data.table(locus = gl$locus, chr = gl$chr, set = s,
      start_Mb = if (nrow(lc)) round(min(lc$start) / 1e6, 2) else NA, end_Mb = if (nrow(lc)) round(max(lc$end) / 1e6, 2) else NA,
      lead = xa$id, lead_Mb = round(xa$pos / 1e6, 3), mlog10p = round(xa$mlp, 2), scale = xa$scale, maf = round(xa$maf, 3), n_snp_ge4 = if (nrow(lc)) sum(lc$n_snp_ge4) else 0L,
      meff_thr = th$meff_thr, bonf_thr = th$bonf_thr, pass_meff = xa$mlp >= th$meff_thr, pass_bonf = xa$mlp >= th$bonf_thr,
      schemeB = if (!is.null(xb)) round(xb$mlp, 2) else NA, schemeBn = if (!is.null(xbn)) round(xbn$mlp, 2) else NA, markers_in_locus = mk_in) } }
LONG <- rbindlist(long, fill = TRUE)
# «да» — пересекает порог (Meff или 4) и падает >= 1 единицы −log10 p; «частично» — пересекает порог при падении < 1
LONG[, cross := (mlog10p >= meff_thr & schemeB < meff_thr) | (mlog10p >= 4 & schemeB < 4)]
LONG[, removed_by_markers := fifelse(is.na(schemeB), NA_character_, fifelse(cross & mlog10p - schemeB >= 1, "да", fifelse(cross, "частично", "нет")))][, cross := NULL]
## прежние стадии (без поправки на происхождение, этап 4/7): MLM/GEMMA и FarmCPU/BLINK, все признаки и шкалы
st4 <- rbindlist(lapply(c("stream1", "stream2"), function(s) fread(file.path(WORK, "gwas", sprintf("all_results_%s.csv.gz", s)), select = c("trait", "scheme", "model", "chr", "pos", "p"))[, stream := s]))
st4[, mg := model %in% c("MLM", "GEMMA")]
EARLY <- rbindlist(lapply(seq_len(nrow(GL)), function(j) { gl <- GL[j]; x <- st4[chr == gl$chr & pos >= gl$start & pos <= gl$end]
  if (!nrow(x)) return(data.table(locus = gl$locus))
  data.table(locus = gl$locus, early_MLM_GEMMA = round(-log10(min(x[mg == TRUE]$p)), 2), early_MLM_GEMMA_at = x[mg == TRUE][which.min(p), paste(stream, trait, scheme, model)],
             early_FarmCPU_BLINK = round(-log10(min(x[mg == FALSE]$p)), 2), early_FB_at = x[mg == FALSE][which.min(p), paste(stream, trait, scheme, model)]) }), fill = TRUE)
## статус
indep <- c("IT", "RU (S1)", "F2")
STAT <- rbindlist(lapply(GL$locus, function(l) { x <- LONG[locus == l & !is.na(mlog10p)]; e <- EARLY[locus == l]
  nb <- x[set %in% indep & pass_bonf == TRUE, .N]; nm <- x[set %in% indep & pass_meff == TRUE, .N]
  any_meff <- any(x$pass_meff); any_bonf <- any(x$pass_bonf); mx <- max(x$mlog10p); n4 <- x[mlog10p >= 4, .N]; nsnp <- max(x$n_snp_ge4)
  rem_mk <- x[removed_by_markers == "да" & markers_in_locus == "", .N] > 0
  early_m <- !is.na(e$early_MLM_GEMMA) && e$early_MLM_GEMMA >= 4.54
  st <- if (nb >= 1 && nm >= 2) "подтверждён" else if (nb >= 1 || (any_meff && n4 >= 2)) "вероятный" else if (any_meff || (mx >= 4 && n4 >= 2)) "предположительный" else
        if ((mx < 4 || !any_meff) && (early_m || (!is.na(e$early_FarmCPU_BLINK) && e$early_FarmCPU_BLINK >= 5.78)) ) "снят поправкой" else if (mx >= 4 && nsnp <= 1) "шум" else if (mx >= 4) "предположительный" else "ниже 4 во всех наборах"
  if (rem_mk && st %in% c("вероятный", "предположительный")) st <- paste0(st, "; снимается маркерами вне локуса")
  data.table(locus = l, status = st, n_indep_bonf = nb, n_indep_meff = nm, max_mlog10p = mx, n_sets_ge4 = n4) }))
STAT <- merge(STAT, EARLY, by = "locus", all.x = TRUE)
LONG <- merge(LONG, STAT[, .(locus, status)], by = "locus")
setorder(LONG, chr, locus, set)
fwrite(LONG, file.path(OUT, "tables/14_loci_long.csv")); fwrite(STAT, file.path(OUT, "tables/14_loci_status.csv"))
W <- dcast(LONG[!is.na(mlog10p)], locus + chr + status ~ set, value.var = "mlog10p")
W <- merge(W, dcast(LONG[!is.na(mlog10p)], locus ~ set, value.var = "lead_Mb"), by = "locus", suffixes = c("", "_leadMb"))
fwrite(W, file.path(OUT, "tables/14_loci_wide.csv"))
print(STAT[order(-max_mlog10p)], nrows = 200)
## 1B подробно
B1 <- ALL[chr == "1B" & pos >= 650196761 & pos <= 660809546, .(mlog10p = round(max(mlp), 2), best = id[which.max(mlp)]), by = .(set, scheme, scale)][order(set, scheme, scale)]
B1e <- st4[chr == "1B" & pos >= 650196761 & pos <= 660809546, .(mlog10p = round(-log10(min(p)), 2)), by = .(stream, trait, scheme, model)][order(-mlog10p)]
fwrite(B1, file.path(OUT, "tables/14_1B_by_set_scheme.csv")); fwrite(B1e, file.path(OUT, "tables/14_1B_stage4_no_origin.csv"))
print(B1, nrows = 100)
