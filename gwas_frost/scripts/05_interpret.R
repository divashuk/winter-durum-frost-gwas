# Этап 5.2. Интерпретация: QTL-интервалы, воспроизводимость, PVE, схема A против B, известные локусы.
# Значимые SNP: Бонферрони (0.05/M), Бонферрони по Meff (0.05/Meff) или FDR(BH) < 0.05 — в любой модели.
# Суггестивные: -log10 p >= 4 (отдельный уровень).
# Интервал QTL: от ведущего SNP — крайние SNP той же хромосомы в пределах ±D (расстояние спада LD до критического r2,
#   по потоку и субгеному), у которых r2 с ведущим SNP >= 0.2 (минимум ±0.5 Мб). Пересекающиеся интервалы
#   (в т. ч. из разных признаков/моделей/потоков) объединяются в один QTL.
suppressPackageStartupMessages(library(data.table))
OUT  <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
R <- rbindlist(lapply(c("stream1", "stream2"), function(s) fread(file.path(WORK, "gwas", sprintf("all_results_%s.csv.gz", s)))))
R[, sig := sig_bonf | sig_meff | sig_fdr]
R[, level := fifelse(sig, "значимый", fifelse(suggestive, "суггестивный", "нет"))]
H <- R[level != "нет"]
cat("значимых записей:", H[level == "значимый", .N], "; суггестивных:", H[level == "суггестивный", .N], "\n")

## аллель эффекта: GAPIT — суффикс в имени SNP (столбец .raw), GEMMA — allele1 = 5-й столбец bim
bims <- rbindlist(lapply(c("stream1", "stream2"), function(s) fread(file.path(WORK, "qc", s, "final.bim"), header = FALSE)[, stream := s]))
bims[, key := paste(V1, V4)]
H[, key := paste(chr, pos)]
H[, effect_allele := fifelse(model == "GEMMA", bims$V5[match(paste(stream, key), paste(bims$stream, bims$key))], sub(".*_", "", SNP))]
H[, snp_id := paste(chr, pos, sep = ":")]

## PVE = var(доза) * beta^2 / var(y) (по признаку и образцам прогона; beta на копию аллеля эффекта)
B <- fread(file.path(OUT, "tables/01_blues.csv"))
geno <- list()
for (s in c("stream1", "stream2")) {
  raw <- fread(file.path(WORK, "qc", s, "final.raw"))
  geno[[s]] <- list(ids = raw$IID, X = as.matrix(raw[, -(1:6)]), key = sub(":[^:]*:[^:]*$", "", sub("_[^_]*$", "", names(raw)[-(1:6)])))
  rm(raw)
}
traitY <- function(s, trait) {
  tn <- sub("_(ang|raw)$", "", trait); sc <- sub(".*_", "", trait)
  tr <- switch(tn, FT12j = c("FT_m12", "ITFT_m12"), FTcol = c("FT_all", "ITFT_m12"), tn)
  B[trait %in% tr & scale == sc, .(TP_MO, y = BLUE)]
}
H[, pve := {
  g <- geno[[stream[1]]]; y <- traitY(stream[1], trait[1]); i <- match(y$TP_MO, g$ids); ok <- !is.na(i)
  j <- match(snp_id, g$key); vy <- var(y$y[ok])
  sapply(seq_len(.N), function(k) if (is.na(j[k]) || is.na(effect[k])) NA_real_ else var(g$X[i[ok], j[k]]) * effect[k]^2 / vy)
}, by = .(stream, trait)]
H[pve > 1, pve := NA]   # FarmCPU/BLINK псевдо-QTN могут давать завышенные эффекты; PVE > 1 не интерпретируем

## интервалы по локальному LD
LDD <- rbindlist(lapply(c("stream1", "stream2"), function(s) fread(file.path(OUT, sprintf("tables/03_ld_decay_%s.csv", s)))))
H[, bin2Mb := round(pos / 2e6)]   # отдельный столбец: столбцы из выражения by выпадают из .SD
lead <- H[level == "значимый", .SD[which.min(p)], by = .(stream, trait, scheme, model, chr, bin2Mb)]
interval <- function(s, ch, ps) {
  g <- geno[[s]]; D <- LDD[stream == s & subgenome == substr(ch, 2, 2), dist_critical_Mb] * 1e6
  kk <- which(startsWith(g$key, paste0(ch, ":"))); pp <- as.numeric(sub(".*:", "", g$key[kk]))
  j0 <- kk[which.min(abs(pp - ps))]; w <- kk[abs(pp - ps) <= D]
  r2 <- suppressWarnings(cor(g$X[, j0], g$X[, w])^2)[1, ]
  inn <- as.numeric(sub(".*:", "", g$key[w[!is.na(r2) & r2 >= 0.2]]))
  c(max(1, min(c(inn, ps - 5e5))), max(c(inn, ps + 5e5)))
}
if (nrow(lead)) {
  iv <- t(mapply(interval, lead$stream, lead$chr, lead$pos)); lead[, `:=`(start = iv[, 1], end = iv[, 2])]
  setorder(lead, chr, start)
  lead[, qtl := {cl <- integer(.N); cur <- 0; e <- -Inf; for (k in seq_len(.N)) { if (start[k] > e) cur <- cur + 1; cl[k] <- cur; e <- max(e, end[k]) }; cl}, by = chr]
  lead[, QTL := sprintf("Q%s.%d", chr, qtl)]
  Q <- lead[, .(start = min(start), end = max(end), lead_snp = snp_id[which.min(p)], min_p = min(p),
                n_models = uniqueN(model), models = paste(sort(unique(model)), collapse = ","),
                n_streams = uniqueN(stream), streams = paste(sort(unique(stream)), collapse = ","),
                n_traits = uniqueN(sub("_(ang|raw)$", "", trait)), traits = paste(sort(unique(trait)), collapse = ","),
                schemes = paste(sort(unique(scheme)), collapse = ",")), by = .(QTL, chr)]
  # градация: воспроизводимый = >=2 модели, среди которых MLM или GEMMA; или >=2 потока; или >=2 разных признака
  Q[, reliability := fifelse((n_models >= 2 & grepl("MLM|GEMMA", models)) | n_streams >= 2 | n_traits >= 2,
                             "воспроизводимый", "единичный")]
} else { Q <- data.table(); lead[, `:=`(start = numeric(), end = numeric(), QTL = character())] }

## таблица всех значимых и суггестивных SNP с привязкой к QTL
H[, QTL := NA_character_]
if (nrow(Q)) for (k in seq_len(nrow(Q))) H[chr == Q$chr[k] & pos >= Q$start[k] & pos <= Q$end[k], QTL := Q$QTL[k]]
## воспроизводимость на уровне SNP
H[, `:=`(rep_models = uniqueN(model[level == "значимый"]), rep_streams = uniqueN(stream[level == "значимый"]),
         rep_traits = uniqueN(sub("_(ang|raw)$", "", trait[level == "значимый"]))), by = snp_id]

## схема A против B: для каждого QTL и сочетания поток×признак×модель — min p в A и в B внутри интервала
AB <- data.table()
if (nrow(Q)) {
  AB <- rbindlist(lapply(seq_len(nrow(Q)), function(k) {
    d <- R[chr == Q$chr[k] & pos >= Q$start[k] & pos <= Q$end[k]]
    d[, .(minp = min(p), sig = any(sig), sugg = any(suggestive)), by = .(stream, trait, model, scheme)][, QTL := Q$QTL[k]]
  }))
  AB <- dcast(AB, QTL + stream + trait + model ~ scheme, value.var = c("minp", "sig", "sugg"))
  AB[, verdict := fifelse(sig_A & !sig_B & !sugg_B, "уходит при маркерах-ковариатах",
                  fifelse(sig_A & !sig_B, "ослабевает до суггестивного",
                  fifelse(sig_A & sig_B, "остаётся", fifelse(!sig_A & sig_B, "появляется только в B", "незначим"))))]
  AB <- AB[sig_A | sig_B]
}

## известные локусы (позиции на Langdon) — из 05_known_loci_positions.csv, если есть
kp <- file.path(OUT, "tables/05_known_loci_positions.csv")
if (file.exists(kp) && nrow(Q)) {
  K <- fread(kp)
  Q[, known_overlap := sapply(seq_len(.N), function(k) {
    D <- LDD[stream == "stream1" & subgenome == substr(chr[k], 2, 2), dist_critical_Mb] * 1e6
    kk <- K[chr == Q$chr[k]]; if (!nrow(kk)) return("")
    d <- pmax(0, pmax(kk$start - Q$end[k], Q$start[k] - kk$end))
    hit <- kk[d <= D]; if (!nrow(hit)) return("")
    paste(sprintf("%s (%.1f Мб)", hit$locus, d[d <= D] / 1e6), collapse = "; ")
  })]
}
fwrite(H[order(level, p)], file.path(OUT, "tables/05_significant_and_suggestive_snps.csv"))
fwrite(Q, file.path(OUT, "tables/05_qtl_intervals.csv"))
fwrite(AB, file.path(OUT, "tables/05_schemeA_vs_B.csv"))
print(Q); print(AB)
