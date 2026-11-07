# Этап 2.1. Гетерозиготность по образцу и IBS-дубликаты (по неимпутированным генотипам после QC).
# Гетерозиготность: доля гетерозигот среди генотипированных сайтов; «аномально высокая» = выше Q3 + 3·IQR
# внутри потока (твёрдая пшеница — самоопылитель). Помечаем, не исключаем.
# Дубликаты: IBS >= 0.98 (plink1.9 --distance ibs). Исключаем один образец из пары только если это один сорт под
# разными названиями (по названиям/родословной); иначе пара остаётся и помечается.
suppressPackageStartupMessages(library(data.table))
OUT  <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
idm  <- fread(file.path(OUT, "tables/01_id_map.csv"), encoding = "UTF-8")
res_het <- list(); res_ibs <- list()
for (s in c("stream1", "stream2")) {
  d  <- file.path(WORK, "qc", s)
  sc <- fread(file.path(d, "s1.scount"))
  setnames(sc, sub("^#", "", names(sc)))
  sc[, called := HOM_REF_CT + HOM_ALT_SNP_CT + HET_SNP_CT]
  sc[, het_rate := HET_SNP_CT / called]
  sm <- fread(file.path(d, "s1.smiss")); setnames(sm, sub("^#", "", names(sm)))
  sc <- merge(sc, sm[, .(IID, F_MISS)], by = "IID")
  q  <- quantile(sc$het_rate, c(.25, .75)); thr <- q[2] + 3 * diff(q)
  sc[, het_flag := het_rate > thr]
  # Распределение непрерывное до ~0.37, правило IQR слишком мягкое -> дополнительные фиксированные пороги:
  # > 0.125 — больше ожидаемой остаточной гетерозиготности линии F4; > 0.20 — аномально высокая (вероятна смесь/гетерогенный образец)
  sc[, het_class := fifelse(het_rate > 0.20, "аномально высокая (>0.20)", fifelse(het_rate > 0.125, "повышенная (0.125–0.20)", "норма (<=0.125)"))]
  sc[, `:=`(stream = s, het_threshold = thr)]
  res_het[[s]] <- sc
  ids <- fread(file.path(d, "s1_ibs.mibs.id"), header = FALSE)$V2
  M   <- as.matrix(fread(file.path(d, "s1_ibs.mibs"), header = FALSE))
  ut  <- which(upper.tri(M) & M >= 0.98, arr.ind = TRUE)
  if (nrow(ut)) res_ibs[[s]] <- data.table(stream = s, id1 = ids[ut[, 1]], id2 = ids[ut[, 2]], IBS = M[ut])
  # распределение максимального IBS на образец (для отчёта)
  diag(M) <- NA
  res_het[[s]][, maxIBS := apply(M, 1, max, na.rm = TRUE)[match(IID, ids)]]
}
H <- rbindlist(res_het)
H <- merge(H, idm[, .(TP_MO, name_pheno, collection)], by.x = "IID", by.y = "TP_MO", all.x = TRUE)
fwrite(H[order(stream, -het_rate)], file.path(OUT, "tables/02_sample_qc_het_missing.csv"))
I <- if (length(res_ibs)) rbindlist(res_ibs) else data.table(stream = character(), id1 = character(), id2 = character(), IBS = numeric())
I <- merge(I, idm[, .(TP_MO, name1 = name_pheno)], by.x = "id1", by.y = "TP_MO", all.x = TRUE)
I <- merge(I, idm[, .(TP_MO, name2 = name_pheno)], by.x = "id2", by.y = "TP_MO", all.x = TRUE)
fwrite(I[order(stream, -IBS)], file.path(OUT, "tables/02_ibs_pairs_ge098.csv"))
cat("Порог гетерозиготности:\n"); print(H[, .(thr = unique(het_threshold), median_het = median(het_rate), n_flag = sum(het_flag)), by = stream])
print(H[, .N, by = .(stream, collection, het_class)][order(stream, collection, het_class)])
cat("\nОбразцы с аномальной гетерозиготностью (IQR):\n"); print(H[het_flag == TRUE,.(stream, IID, name_pheno, collection, het_rate = round(het_rate, 3), F_MISS = round(F_MISS, 3))])
cat("\nПары IBS >= 0.98:\n"); print(I[, .(stream, id1, name1, id2, name2, IBS = round(IBS, 4))])
cat("\nМакс. IBS на образец, квантили:\n"); print(H[, as.list(round(quantile(maxIBS, c(0, .5, .9, .99, 1)), 4)), by = stream])
