# Этап 9.5. Кандидаты Q5B.1, Q5B.2, Q1B.6, Q4A.5 внутри коллекций.
#  - частоты аллелей ведущих SNP у российских (S1, 144) и итальянских (182) образцов;
#  - −log10 p, эффект и знак в GWAS внутри коллекций (09_gwas.R: RU — FT_m12 и S1 FT_all, схема A; IT — ITFT_m12, PC + родство);
#  - взаимодействие SNP × коллекция в S2: смешанная модель (rrBLUP::mixed.solve, родство VanRaden потока 2),
#    фиксированные: происхождение (дамми, включая итальянскую коллекцию), PC (число по BIC S2), SNP, SNP×коллекция; тест Вальда.
suppressPackageStartupMessages({library(data.table); library(rrBLUP)})
OUT <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
LEAD <- data.table(QTL = c("Q5B.1", "Q5B.2", "Q5B.2", "Q1B.6", "Q4A.5"), chr = c("5B", "5B", "5B", "1B", "4A"),
                   pos = c(521615359, 535618213, 532482598, 653803203, 642497807))
raw <- fread(file.path(WORK, "qc/stream2/final.raw")); ids <- raw$IID; X <- as.matrix(raw[, -(1:6)]); rm(raw)
key <- sub(":[^:]*:[^:]*_[^_]*$", "", colnames(X)); cnt <- sub(".*_", "", colnames(X))
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8")
coll <- OR$collection[match(ids, OR$TP_MO)]
LEAD[, j := match(paste(chr, pos, sep = ":"), key)]
FR <- LEAD[, .(QTL, snp = colnames(X)[j], counted_allele = cnt[j],
               freq_RU = round(mean(X[coll == "RU_winter", j], na.rm = TRUE) / 2, 3),
               freq_IT = round(mean(X[coll == "italy", j], na.rm = TRUE) / 2, 3)), by = seq_len(nrow(LEAD))][, seq_len := NULL]
## GWAS внутри коллекций и в S2 (схема A), MLM и GEMMA
rd <- function(s) { f <- file.path(WORK, "gwas09", sprintf("all_results_09_%s.csv.gz", s)); if (file.exists(f)) fread(f) else NULL }
RR <- rbindlist(lapply(c("S1", "RU", "IT", "S2"), rd), fill = TRUE)
W <- rbindlist(lapply(seq_len(nrow(LEAD)), function(k) RR[scheme == "A" & model %in% c("MLM", "GEMMA") & chr == LEAD$chr[k] & pos == LEAD$pos[k] &
      trait %in% c("FT_all_ang", "FT_m12_ang", "ITFT_m12_ang", "FT12j_rint", "FTcol_rint"),
  .(QTL = LEAD$QTL[k], snp_pos = paste(chr, pos, sep = ":"), set, trait, model, n, mlog10p = round(-log10(p), 2), effect = signif(effect, 3), effect_allele)]))
# GAPIT-эффект — на копию суффиксного аллеля SNP, GEMMA — на allele1; приводим знак к считаемому аллелю .raw (суффикс)
W[, effect_counted := fifelse(effect_allele == FR$counted_allele[match(QTL, FR$QTL)], effect, -effect)]
W[, sign := fifelse(effect_counted > 0, "+", "−")]
## взаимодействие SNP × коллекция
K <- as.matrix(fread(file.path(WORK, "kin/K_vanraden_stream2.csv"))[, -1]); rownames(K) <- colnames(K)
S2T <- fread(file.path(OUT, "tables/09_s2_traits.csv"))
bic <- function(t) { f <- file.path(OUT, sprintf("tables/09_bic_S2_%s.csv", t)); if (!file.exists(f)) return(1); b <- fread(f); as.integer(b[[1]][which.max(b[[2]])]) }
INT <- list()
for (tr in c("FT12j", "FTcol")) for (v in c("rint", "ang")) {
  d <- S2T[trait == tr & TP_MO %in% ids]; yv <- if (v == "rint") d$rint else d$y; ii <- match(d$TP_MO, ids)
  npc <- bic(paste(tr, v, sep = "_")); Xs <- X[ii, ]; Xs <- Xs[, apply(Xs, 2, sd, na.rm = TRUE) > 0]
  pcs <- if (npc > 0) prcomp(Xs, rank. = npc)$x else NULL
  o <- factor(OR$origin[match(d$TP_MO, OR$TP_MO)]); D <- model.matrix(~ o)
  it <- as.numeric(coll[ii] == "italy"); Kd <- K[d$TP_MO, d$TP_MO]
  for (k in seq_len(nrow(LEAD))) {
    g <- X[ii, LEAD$j[k]]; g[is.na(g)] <- mean(g, na.rm = TRUE)
    Xf <- cbind(D, pcs, snp = g, snp_x_italy = g * it)
    fit <- mixed.solve(yv, X = Xf, K = Kd, SE = TRUE)
    b <- fit$beta[c("snp", "snp_x_italy")]; se <- fit$beta.SE[c("snp", "snp_x_italy")]
    INT[[length(INT) + 1]] <- data.table(QTL = LEAD$QTL[k], snp_pos = paste(LEAD$chr[k], LEAD$pos[k], sep = ":"), trait = paste(tr, v, sep = "_"),
      beta_snp_RU = signif(b[1], 3), beta_snp_IT = signif(b[1] + b[2], 3), beta_interaction = signif(b[2], 3),
      p_interaction = signif(2 * pnorm(-abs(b[2] / se[2])), 3), p_snp_main = signif(2 * pnorm(-abs(b[1] / se[1])), 3))
  }
}
INT <- rbindlist(INT)
fwrite(FR, file.path(OUT, "tables/09_candidates_allele_freq.csv"))
fwrite(W, file.path(OUT, "tables/09_candidates_within_collections.csv"))
fwrite(INT, file.path(OUT, "tables/09_candidates_interaction.csv"))
print(FR); print(W[order(QTL, snp_pos, set, trait, model)], nrows = 200); print(INT)
