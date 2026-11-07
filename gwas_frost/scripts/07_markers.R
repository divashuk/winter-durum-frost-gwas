# Этап 7 (ТЗ, раздел 2). Маркеры durum_frost как отдельные ковариаты.
#  1) опись всех маркеров: локус, тип, кодировка, число типированных в S1 (161 / 144 в GWAS) и у итальянцев (200 / 182);
#  2) импутация пропусков KASP по фланкирующим SNP GBS (±5 Мб от позиции маркера на Langdon):
#     логистическая регрессия elastic net (glmnet, alpha = 0.5), обучение на образцах S1 с вызовом, 5 × 10-кратная CV;
#     дополнительно для маркеров, типированных в обеих коллекциях, — перенос «обучение на РФ → проверка на итальянцах»;
#     импутированное значение используется только при точности CV >= 0.9 (порог ТЗ);
#  3) итоговые матрицы ковариат схемы B (каждый маркер отдельно, без гаплотипов) для S1 и S2, r2 и VIF (контроль).
suppressPackageStartupMessages({library(data.table); library(glmnet)})
OUT  <- "/work/project/gwas_frost"
SRC  <- "/work/phenotypes"
WORK <- "/work/scratch"
set.seed(20261025)
idm <- fread(file.path(OUT, "tables/01_id_map.csv"), encoding = "UTF-8")
g161 <- fread(file.path(SRC, "08_deposit/data/genotypes_161x11.csv"), encoding = "UTF-8")
gv2  <- fread(file.path(SRC, "01_data/markers/genotypes_clean_v2.csv"), encoding = "UTF-8")

## ---------- 1. опись и кодировка ----------
# marker, локус, тип, класс -> 2, позиция на Langdon (по праймерам или гену; * — по максимуму LD), источник РФ, источник ИТ
INV <- rbindlist(list(
  list("S2269949", "Fr-A2, CBF-A14 CNV-тег (5A)", "KASP", "FAM", "5A", 524288011, "праймеры", "Табл. S1", "genotypes_clean_v2"),
  list("S1077313", "Fr-A2, CBF-A14 CNV-тег (5A)", "KASP", "HEX", "5A", 537231599, "праймеры", "Табл. S1", "genotypes_clean_v2"),
  list("S1862541", "Fr-A2, CBF-A14 CNV-тег (5A)", "KASP", "HEX", "5A", 530965192, "праймеры", "Табл. S1", "genotypes_clean_v2"),
  list("S1298957", "Fr-A2, CBF-A14 CNV-тег (5A)", "KASP", "HEX", "5A", 513240733, "праймеры", "Табл. S1", "genotypes_clean_v2"),
  list("S1051014", "Fr-A2, не-CNV компонент (5A)", "KASP", "FAM", "5A", 546128492, "LD*", "Табл. S1", "genotypes_clean_v2"),
  list("Vrn1_new", "Fr-A1/VRN-A1, экзон 4 (5A)", "KASP", "2147", "5A", 598268462, "ген VRN-A1", "Табл. S1", "genotypes_clean_v2"),
  list("VrnA1_9K0001", "VRN-A1, тип развития (5A)", "KASP", "vrn-A1 Winter", "5A", 598268462, "ген VRN-A1", "Табл. S1", "genotypes_clean_v2"),
  list("Vrn-A1b-Marq", "VRN-A1, Vrn-A1b Marquis (5A); в статью не вошёл", "KASP", "Vrn-A1b Marquis", "5A", 598268462, "ген VRN-A1", "genotypes_clean_v2", "genotypes_clean_v2"),
  list("Exon7", "VRN-A1, экзон 7 C/T (5A); в статью не вошёл", "KASP", "Hereward", "5A", 598268462, "ген VRN-A1", "genotypes_clean_v2", "genotypes_clean_v2"),
  list("Vrn-B1", "VRN-B1 (5B)", "ПЦР + форез", "vrn-B1", "5B", 605129460, "ген VRN-B1", "Табл. S1", "genotypes_clean_v2"),
  list("Mor534", "не-Fr локус 5B (5B_517276534, Svevo v2)", "KASP", "G", "5B", 551172879, "LD* (праймеры пары не дали)", "Табл. S1", "не типирован"),
  list("Mor587", "не-Fr локус 1B (1B_41099587, Svevo v2)", "KASP", "T", "1B", 54736990, "праймеры", "Табл. S1", "не типирован"),
  list("Mor866", "не-Fr локус 7B (7B_598228866, Svevo v2)", "KASP", "A", "7B", 610783371, "праймеры", "Табл. S1", "не типирован")))
setnames(INV, c("marker", "locus", "type", "class_coded_2", "chr", "pos_langdon", "pos_method", "src_RU", "src_IT"))
code <- function(v, cls) {
  s <- trimws(as.character(v)); out <- rep(NA_real_, length(s))
  het <- grepl("heterozyg", s, ignore.case = TRUE); bad <- is.na(v) | s %in% c("", "nan", "NA") | grepl("heterogen", s, ignore.case = TRUE)
  out[!bad] <- ifelse(startsWith(toupper(s[!bad]), toupper(cls)), 2, 0); out[het] <- 1; out[bad] <- NA; out
}
ru <- idm[collection == "RU_winter" & in_article == TRUE, .(TP_MO, sid = sample_id_TP)]
it <- idm[collection == "italy", .(TP_MO, sid = sample_id_TP)]
w2 <- gv2[collection == "winter"]; i2 <- gv2[collection == "italy"]
M <- rbind(ru[, .(TP_MO, collection = "RU")], it[, .(TP_MO, collection = "IT")])
for (k in seq_len(nrow(INV))) {
  mk <- INV$marker[k]; cls <- INV$class_coded_2[k]
  vr <- if (mk %in% names(g161)) g161[[mk]][match(ru$sid, g161$sample_id)] else w2[[mk]][match(ru$sid, w2$sample_id)]
  vi <- if (mk %in% names(i2)) i2[[mk]][match(it$sid, i2$sample_id)] else rep(NA, nrow(it))
  M[, (mk) := c(code(vr, cls), code(vi, cls))]
}
M[, in_gwas := TP_MO %in% idm[!grepl("^исключён", role), TP_MO]]
cnt <- function(x) sum(!is.na(x))
poly <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(NA_real_); p <- mean(x) / 2; min(p, 1 - p) }
INV[, `:=`(
  n_S1_161 = sapply(marker, function(m) cnt(M[collection == "RU", get(m)])),
  n_S1_gwas144 = sapply(marker, function(m) cnt(M[collection == "RU" & in_gwas, get(m)])),
  MAF_S1 = round(sapply(marker, function(m) poly(M[collection == "RU" & in_gwas, get(m)])), 3),
  n_IT_200 = sapply(marker, function(m) cnt(M[collection == "IT", get(m)])),
  n_IT_gwas182 = sapply(marker, function(m) cnt(M[collection == "IT" & in_gwas, get(m)])),
  MAF_IT = round(sapply(marker, function(m) poly(M[collection == "IT" & in_gwas, get(m)])), 3))]
INV[, previous_schemeB := fcase(
  marker %in% c("S2269949", "S1077313", "S1862541", "S1298957"), "да, но сведены в FrA2_score (среднее 4 доз) — S1 и S2; в S1 FT_m15/16/17 выпал (мономорфен)",
  marker %in% c("S1051014", "Vrn1_new"), "да — S1 и S2",
  marker %in% c("Mor534", "Mor587", "Mor866"), "да — только S1",
  default = "нет")]

## ---------- 2. импутация по GBS ----------
raw <- fread(file.path(WORK, "qc/stream2/final.raw"))
gid <- raw$IID; X <- as.matrix(raw[, -(1:6)]); key <- names(raw)[-(1:6)]; rm(raw)
kchr <- sub(":.*", "", key); kpos <- as.numeric(sapply(strsplit(key, ":"), `[`, 2))
knn3 <- function(xtr, ytr, xte) {   # запасной вариант для редких классов: 3 ближайших по манхэттенскому расстоянию
  apply(xte, 1, function(z) { dd <- colSums(abs(t(xtr) - z)); nn <- ytr[order(dd)[1:3]]; as.numeric(names(which.max(table(nn)))) })
}
fitpred <- function(xtr, ytr, xte) {
  yb <- factor(ytr >= 1, levels = c(FALSE, TRUE))
  if (min(table(yb)) < 3) return(rep(as.numeric(names(which.max(table(ytr)))), nrow(xte)))  # почти мономорфен
  tryCatch({
    fid <- integer(length(yb)); for (lv in levels(yb)) { i <- which(yb == lv); fid[i] <- sample(rep(1:5, length.out = length(i))) }
    cv <- cv.glmnet(xtr, yb, family = "binomial", alpha = 0.5, foldid = fid)
    ifelse(predict(cv, xte, s = "lambda.1se", type = "response")[, 1] > 0.5, 2, 0)
  }, error = function(e) knn3(xtr, ytr, xte))
}
ACC <- list(); IMP <- copy(M)
for (k in seq_len(nrow(INV))) {
  mk <- INV$marker[k]
  win <- which(kchr == INV$chr[k] & abs(kpos - INV$pos_langdon[k]) <= 5e6)
  y <- M[[mk]][match(gid, M$TP_MO)]; coll <- M$collection[match(gid, M$TP_MO)]
  typed_IT <- sum(!is.na(y[coll == "IT"])) >= 20
  # обучение: только РФ (S1), если у итальянцев маркер не типирован (Mor); иначе — все типированные (РФ + ИТ)
  tr <- which(coll %in% (if (typed_IT) c("RU", "IT") else "RU") & !is.na(y) & y != 1)
  tr_ru <- which(coll == "RU" & !is.na(y) & y != 1)
  if (length(tr) < 30 || length(unique(y[tr])) < 2 || length(win) < 5) {
    ACC[[mk]] <- data.table(marker = mk, n_snp_window = length(win), n_train_RU = length(tr), note = "не импутируется: мономорфен в S1 или мало данных"); next
  }
  xw <- X[, win, drop = FALSE]
  cvacc <- c(); minor <- c()
  for (rep in 1:5) {
    fold <- integer(length(tr)); for (lv in unique(y[tr])) { i <- which(y[tr] == lv); fold[i] <- sample(rep(1:10, length.out = length(i))) }  # стратифицированно
    pr <- rep(NA, length(tr))
    for (f in 1:10) pr[fold == f] <- fitpred(xw[tr[fold != f], ], y[tr[fold != f]], xw[tr[fold == f], , drop = FALSE])
    cvacc <- c(cvacc, mean(pr == y[tr]))
    mc <- as.numeric(names(which.min(table(y[tr])))); minor <- c(minor, mean(pr[y[tr] == mc] == mc))
  }
  base <- max(table(y[tr])) / length(tr)
  # перенос РФ -> ИТ (только если маркер типирован у итальянцев)
  te <- which(coll == "IT" & !is.na(y) & y != 1); tacc <- NA; tbase <- NA
  if (length(te) >= 20) { p <- fitpred(xw[tr_ru, ], y[tr_ru], xw[te, , drop = FALSE]); tacc <- mean(p == y[te]); tbase <- max(table(y[te])) / length(te) }
  acc <- mean(cvacc); use <- acc >= 0.9 && mean(minor) >= 0.5   # доп. условие: модель лучше «всегда мажорный класс»
  ACC[[mk]] <- data.table(marker = mk, n_snp_window = length(win), train_set = if (typed_IT) "РФ + ИТ" else "только РФ (S1)", n_train = length(tr), cv_accuracy = round(acc, 3),
                          cv_accuracy_sd = round(sd(cvacc), 3), cv_minor_class_recall = round(mean(minor), 3),
                          majority_baseline = round(base, 3), check_transfer_RU_to_IT_accuracy = round(tacc, 3),
                          check_majority_baseline_IT = round(tbase, 3), use_imputation = use,
                          note = fifelse(acc >= 0.9 & !use, "точность >= 0.9 только за счёт мажорного класса (полнота минорного < 0.5) — не используется; пропуски = среднее внутри коллекции", ""))
  if (use) {   # заполняем все пропуски этого маркера (РФ и ИТ) предсказанием по модели на всех типированных РФ
    miss <- which(is.na(y) & coll %in% c("RU", "IT"))
    if (length(miss)) { p <- fitpred(xw[tr, ], y[tr], xw[miss, , drop = FALSE]); IMP[match(gid[miss], TP_MO), (mk) := p]
      IMP[match(gid[miss], TP_MO), paste0(mk, "_src") := "GBS-импутация"] }
  }
}
A <- rbindlist(ACC, fill = TRUE)
fwrite(A, file.path(OUT, "tables/07_marker_imputation_accuracy.csv")); print(A)
for (mk in INV$marker) { s <- paste0(mk, "_src"); if (!s %in% names(IMP)) IMP[, (s) := NA_character_]
  IMP[is.na(get(s)) & !is.na(M[[mk]]), (s) := "KASP/ПЦР"] }
fwrite(IMP, file.path(OUT, "tables/07_markers_final.csv"))

## ---------- 3. наборы ковариат схемы B, r2 и VIF ----------
B1 <- c("S2269949", "S1077313", "S1862541", "S1298957", "S1051014", "Vrn1_new", "Vrn-A1b-Marq", "Exon7", "Mor534", "Mor587", "Mor866")
B2 <- c(B1, "VrnA1_9K0001")
s1 <- readLines(file.path(OUT, "tables/02_samples_stream1.txt")); s2 <- readLines(file.path(OUT, "tables/02_samples_stream2.txt"))
covset <- function(ids, mks) {
  d <- IMP[match(ids, TP_MO), ..mks]
  keep <- sapply(d, function(x) sum(!is.na(x)) >= 10 && sd(x, na.rm = TRUE) > 0)
  d <- d[, ..keep]
  nas <- sapply(d, function(x) sum(is.na(x)))
  cl <- IMP$collection[match(ids, IMP$TP_MO)]
  for (m in names(d)) for (cc in unique(cl)) {   # оставшиеся пропуски — среднее внутри коллекции
    i <- which(is.na(d[[m]]) & cl == cc); if (length(i)) set(d, i, m, mean(d[[m]][cl == cc], na.rm = TRUE)) }
  for (m in names(d)) set(d, which(is.na(d[[m]])), m, mean(d[[m]], na.rm = TRUE))
  list(d = d, nas = nas, dropped = mks[!keep])
}
ctl <- list(); r2l <- list()
for (S in list(list("S1", s1, B1), list("S2 (все 326)", s2, B2))) {
  cs <- covset(S[[2]], S[[3]]); d <- as.matrix(cs$d)
  fwrite(data.table(TP_MO = S[[2]], cs$d), file.path(OUT, sprintf("tables/07_covariates_%s.csv", substr(S[[1]], 1, 2))))   # для 04_gwas.R
  if (S[[1]] != "S1") d <- cbind(d, italy = as.numeric(IMP$collection[match(S[[2]], IMP$TP_MO)] == "IT"))
  r2 <- cor(d)^2; r2l[[S[[1]]]] <- data.table(set = S[[1]], m1 = rownames(r2)[row(r2)[upper.tri(r2)]], m2 = colnames(r2)[col(r2)[upper.tri(r2)]], r2 = round(r2[upper.tri(r2)], 3))
  vif <- sapply(colnames(d), function(m) { o <- setdiff(colnames(d), m); rs <- summary(lm(d[, m] ~ d[, o]))$r.squared; if (rs >= 1) Inf else 1 / (1 - rs) })
  ctl[[S[[1]]]] <- data.table(set = S[[1]], marker = colnames(d), n_missing_mean_imputed = c(cs$nas, italy = 0)[colnames(d)], VIF = round(vif, 2), VIF_gt_10 = vif > 10)
  cat(S[[1]], ": исключены (мономорфны/нет данных):", paste(cs$dropped, collapse = ", "), "\n")
}
V <- rbindlist(ctl); R2 <- rbindlist(r2l)
fwrite(V, file.path(OUT, "tables/07_marker_vif.csv")); fwrite(R2[order(set, -r2)], file.path(OUT, "tables/07_marker_r2.csv"))
fwrite(INV, file.path(OUT, "tables/07_markers_inventory.csv"))
print(INV[, .(marker, n_S1_161, n_S1_gwas144, MAF_S1, n_IT_200, n_IT_gwas182, MAF_IT)]); print(V); print(R2[r2 > 0.5])
