# Этап 9.9 (сводка). QTL после поправки на происхождение: интервалы по локальному LD (как этап 5), привязка к QTL этапа 6,
# минимальные p по схемам A / B / Bn, отдельно MLM/GEMMA и FarmCPU/BLINK, по наборам S1, S2, RU, IT.
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
R <- rbindlist(lapply(c("S1", "S2", "RU", "IT"), function(s) fread(file.path(WORK, "gwas09", sprintf("all_results_09_%s.csv.gz", s)))), fill = TRUE)
R[, sig := sig_bonf | sig_meff | sig_fdr]
R[, mg := fifelse(model %in% c("MLM", "GEMMA"), "MLM/GEMMA", "FarmCPU/BLINK")]
LDD <- rbindlist(lapply(c("stream1", "stream2"), function(s) fread(file.path(OUT, sprintf("tables/03_ld_decay_%s.csv", s)))))
geno <- list()
for (s in c("stream1", "stream2")) { raw <- fread(file.path(WORK, "qc", s, "final.raw"))
  geno[[s]] <- list(X = as.matrix(raw[, -(1:6)]), key = sub(":[^:]*:[^:]*$", "", sub("_[^_]*$", "", names(raw)[-(1:6)]))); rm(raw) }
str_of <- function(set) if (set %in% c("S1", "RU")) "stream1" else "stream2"
interval <- function(set, ch, ps) {
  s <- str_of(set); g <- geno[[s]]; D <- LDD[stream == s & subgenome == substr(ch, 2, 2), dist_critical_Mb] * 1e6
  kk <- which(startsWith(g$key, paste0(ch, ":"))); pp <- as.numeric(sub(".*:", "", g$key[kk]))
  j0 <- kk[which.min(abs(pp - ps))]; w <- kk[abs(pp - ps) <= D]
  r2 <- suppressWarnings(cor(g$X[, j0], g$X[, w])^2)[1, ]
  inn <- as.numeric(sub(".*:", "", g$key[w[!is.na(r2) & r2 >= 0.2]])); c(max(1, min(c(inn, ps - 5e5))), max(c(inn, ps + 5e5))) }
H <- R[sig == TRUE]; H[, bin2Mb := round(pos / 2e6)]
lead <- H[, .SD[which.min(p)], by = .(set, trait, scheme, model, chr, bin2Mb)]
iv <- t(mapply(interval, lead$set, lead$chr, lead$pos)); lead[, `:=`(start = iv[, 1], end = iv[, 2])]
setorder(lead, chr, start)
lead[, qtl := {cl <- integer(.N); cur <- 0; e <- -Inf; for (k in seq_len(.N)) { if (start[k] > e) cur <- cur + 1; cl[k] <- cur; e <- max(e, end[k]) }; cl}, by = chr]
Q6 <- fread(file.path(OUT, "tables/06_qtl_summary.csv"))
lead[, Q := sprintf("%s.%d", chr, qtl)]
Q <- lead[, .(start = min(start), end = max(end), lead_snp = SNP[which.min(p)], min_p = min(p),
              mg_models = paste(sort(unique(model[mg == "MLM/GEMMA"])), collapse = ","), fb_models = paste(sort(unique(model[mg != "MLM/GEMMA"])), collapse = ","),
              sets = paste(sort(unique(set)), collapse = ","), traits = paste(sort(unique(trait)), collapse = ","), schemes = paste(sort(unique(scheme)), collapse = ",")), by = .(Q, chr)]
qc <- Q$chr; qs <- Q$start; qe <- Q$end   # внутри j столбец Q затеняет таблицу Q
Q[, previous_QTL := sapply(seq_along(qc), function(k) paste(Q6[chr == qc[k] & start <= qe[k] & end >= qs[k], QTL], collapse = ","))]
Q[, mg_meff_any := sapply(seq_along(qc), function(k) R[chr == qc[k] & pos >= qs[k] & pos <= qe[k] & mg == "MLM/GEMMA" & (sig_meff | sig_bonf), .N] > 0)]
fwrite(Q, file.path(OUT, "tables/09_qtl_intervals.csv"))
## минимальные p (−log10) по ключевым интервалам и наборам/схемам
KEY <- data.table(locus = c("5A Fr-A2 (Q5A.1)", "5B Q5B.1", "5B Q5B.2", "1B Q1B.6", "4A Q4A.5", "5B Q5B.3 (VRN-B1)", "5A Q5A.4 (у VRN-A1)", "2A Q2A.1 (Ppd-A1)"),
                  chr = c("5A", "5B", "5B", "1B", "4A", "5B", "5A", "2A"),
                  start = c(466363940, 520973288, 525925602, 650196761, 617944415, 580408612, 609173502, 30767399),
                  end = c(545501466, 522742783, 539497626, 660809546, 650214537, 608478255, 622369126, 58925206))
KM <- rbindlist(lapply(seq_len(nrow(KEY)), function(k) R[chr == KEY$chr[k] & pos >= KEY$start[k] & pos <= KEY$end[k],
  .(mlp = round(-log10(min(p)), 2), meff_thr = round(-log10(thr_meff[1]), 2)), by = .(set, trait, scheme, mg)][, locus := KEY$locus[k]]))
W <- dcast(KM, locus + set + trait + mg ~ scheme, value.var = "mlp")
W <- merge(W, unique(KM[, .(locus, set, trait, mg, meff_thr)]), by = c("locus", "set", "trait", "mg"))
fwrite(W, file.path(OUT, "tables/09_key_loci_by_scheme.csv"))
print(Q[mg_meff_any == TRUE]); print(W[mg == "MLM/GEMMA" & trait %in% c("FT_all_ang", "FT12j_rint", "FTcol_rint", "FT_m12_ang", "ITFT_m12_ang", "FTcol_ang")], nrows = 200)
