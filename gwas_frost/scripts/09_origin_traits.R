# Этап 9.0–9.3. Происхождение (ковариата), Cramér's V «происхождение × маркеры», VIF с происхождением,
# новые признаки: LT50 (S1), варианты S2 — ранговое обратное нормальное преобразование внутри коллекции и
# бинарный «выжило хотя бы одно растение».
# Происхождение S1 — Табл. S4 статьи (originator, country, origin group); итальянская коллекция — одна категория
# (поимённого происхождения в источниках нет: итальянский xlsx, абоминация, Табл. S3 Habib 2026 недоступна).
# Правило объединения групп < 5: присоединить к ближайшей по географии группе с >= 5 образцами.
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"
idm <- fread(file.path(OUT, "tables/01_id_map.csv"), encoding = "UTF-8")
idm[, origin := fcase(
  collection == "italy", "итальянская коллекция",
  origin_group == "NCZ Krasnodar", "НЦЗ им. Лукьяненко (Краснодар)",
  origin_group == "ANC Donskoy", "АНЦ Донской (Ростовская обл.)",
  origin_group == "other Russian" | country == "Azerbaijan", "другие российские и Закавказье",
  country %in% c("Ukraine", "Russia/Ukraine"), "Украина",
  country %in% c("Hungary", "Romania", "France", "Austria", "Germany"), "Центральная и Западная Европа",
  default = NA_character_)]
s1 <- readLines(file.path(OUT, "tables/02_samples_stream1.txt")); s2 <- readLines(file.path(OUT, "tables/02_samples_stream2.txt"))
OR <- idm[TP_MO %in% s2, .(TP_MO, collection, name = name_pheno, country, originator, origin_group, origin, in_S1 = TP_MO %in% s1)]
fwrite(OR, file.path(OUT, "tables/09_origin_by_sample.csv"))
print(OR[, .N, by = .(origin, in_S1)][order(-in_S1, -N)])
RULE <- data.table(source_group = c("Romania (4)", "France (2)", "Austria (1)", "Germany (1)", "Hungary (7)", "Azerbaijan (1)", "Russia/Ukraine (1)",
                                    "Prikumskaya, Dagestan, Koshelevsky Posad (6)"),
                   merged_into = c(rep("Центральная и Западная Европа", 5), "другие российские и Закавказье", "Украина", "другие российские и Закавказье"))
fwrite(RULE, file.path(OUT, "tables/09_origin_merge_rule.csv"))

## ---------- Cramér's V: происхождение × аллели маркеров (только вызовы KASP/ПЦР, без импутированных) ----------
mk <- fread(file.path(OUT, "tables/09_markers_final.csv"))
MK <- c("S2269949", "S1077313", "S1862541", "S1298957", "S1051014", "Vrn1_new", "VrnA1_9K0001", "Vrn-A1b-Marq", "Exon7", "Vrn-B1", "Mor534", "Mor587", "Mor866")
cv <- list()
for (S in list(list("S1", s1), list("S2", s2))) for (m in MK) {
  x <- mk[match(S[[2]], TP_MO)]; g <- x[[m]]; src <- x[[paste0(m, "_src")]]; o <- OR$origin[match(S[[2]], OR$TP_MO)]
  ok <- !is.na(g) & src == "KASP/ПЦР" & !is.na(o)
  if (sum(ok) < 10 || length(unique(g[ok])) < 2) { cv[[length(cv) + 1]] <- data.table(set = S[[1]], marker = m, n = sum(ok), note = "мономорфен или нет данных"); next }
  tb <- table(o[ok], g[ok]); ch <- suppressWarnings(chisq.test(tb))
  V <- sqrt(ch$statistic / (sum(tb) * (min(dim(tb)) - 1)))
  cv[[length(cv) + 1]] <- data.table(set = S[[1]], marker = m, n = sum(ok), chi2 = round(ch$statistic, 2), df = ch$parameter, p = signif(ch$p.value, 3),
                                     cramers_V = round(V, 3), table = paste(apply(tb, 1, function(r) paste(r, collapse = "/")), rownames(tb), sep = ": ", collapse = " | "),
                                     classes = paste(colnames(tb), collapse = "/"), note = "")
}
CV <- rbindlist(cv, fill = TRUE); fwrite(CV, file.path(OUT, "tables/09_origin_x_markers.csv")); print(CV[, .(set, marker, n, cramers_V, p)])

## ---------- VIF: маркеры + происхождение (дамми) ----------
vif <- function(d) sapply(colnames(d), function(m) { rs <- summary(lm(d[, m] ~ d[, setdiff(colnames(d), m)]))$r.squared; if (rs >= 1) Inf else 1 / (1 - rs) })
VV <- list()
for (S in list(list("S1", s1, "S1"), list("S2", s2, "S2"))) {
  cov <- fread(file.path(OUT, sprintf("tables/09_covariates_%s.csv", S[[3]])))[match(S[[2]], TP_MO)]
  o <- factor(OR$origin[match(S[[2]], OR$TP_MO)]); D <- model.matrix(~ o)[, -1, drop = FALSE]; colnames(D) <- paste0("orig:", levels(o)[-1])
  X <- cbind(as.matrix(cov[, -1]), D); X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
  VV[[S[[1]]]] <- data.table(set = S[[1]], covariate = colnames(X), VIF = round(vif(X), 2))[, VIF_gt_10 := VIF > 10]
}
V <- rbindlist(VV); fwrite(V, file.path(OUT, "tables/09_vif_markers_origin.csv")); print(V[VIF_gt_10 == TRUE])

## ---------- признаки ----------
B <- fread(file.path(OUT, "tables/01_blues.csv"))
pl <- fread(file.path(OUT, "tables/01_pheno_long.csv"), encoding = "UTF-8")[use == "да"]
# LT50 по BLUE на исходной шкале для каждой температуры (среды внутри температуры — случайный эффект года)
tmp <- c(FT_m12 = -12, FT_m14 = -14, FT_m15 = -15, FT_m16 = -16, FT_m17 = -17)
L <- B[trait %in% names(tmp) & scale == "raw", .(TP_MO, temp = tmp[trait], p = pmin(pmax(BLUE / 100, 0.005), 0.995))]
LT <- L[, { if (uniqueN(temp) >= 3) { m <- tryCatch(glm(p ~ temp, quasibinomial), error = function(e) NULL)
              if (!is.null(m) && coef(m)[2] > 0) .(n_temp = uniqueN(temp), slope = coef(m)[2], LT50 = -coef(m)[1] / coef(m)[2]) else .(n_temp = uniqueN(temp), slope = if (is.null(m)) NA_real_ else coef(m)[2], LT50 = NA_real_) }
            else .(n_temp = uniqueN(temp), slope = NA_real_, LT50 = NA_real_) }, by = TP_MO]
LT[, extrapolated := !is.na(LT50) & (LT50 < -17 | LT50 > -12)]
fwrite(LT, file.path(OUT, "tables/09_lt50.csv"))
cat("\nLT50: генотипов с >= 3 температурами:", LT[n_temp >= 3, .N], "; из них с положительным наклоном (LT50 оценён):", LT[!is.na(LT50), .N],
    "; в S1 (с генотипами):", LT[!is.na(LT50) & TP_MO %in% s1, .N], "; LT50 вне диапазона −17…−12 °C:", LT[extrapolated == TRUE, .N], "\n")
print(summary(LT[!is.na(LT50), LT50]))
print(LT[n_temp >= 3, .N, by = n_temp])
# S2: ранговое обратное нормальное преобразование внутри коллекции и бинарный признак
rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
mkS2 <- function(tr_ru, envs_ru) {
  ru <- B[trait == tr_ru & scale == "ang", .(TP_MO, y = BLUE, coll = "RU")]; it <- B[trait == "ITFT_m12" & scale == "ang", .(TP_MO, y = BLUE, coll = "IT")]
  d <- rbind(ru, it); d[, rint := rint(y), by = coll]
  any_s <- pl[(collection == "RU_winter" & env %in% envs_ru) | (collection == "italy"), .(bin = as.integer(any(value > 0))), by = TP_MO]
  merge(d, any_s, by = "TP_MO", all.x = TRUE)
}
envs_all <- unique(pl[collection == "RU_winter", env]); envs_12 <- c("RU_2021_-12", "RU_2022_-12", "RU_2023_-12")
T12 <- mkS2("FT_m12", envs_12)[, trait := "FT12j"]; TC <- mkS2("FT_all", envs_all)[, trait := "FTcol"]
TT <- rbind(T12, TC)
fwrite(TT, file.path(OUT, "tables/09_s2_traits.csv"))
print(TT[TP_MO %in% s2, .(n = .N, bin_share1 = round(mean(bin), 3), rint_sd = round(sd(rint), 2)), by = .(trait, coll)])
