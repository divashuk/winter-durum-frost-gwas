# Этап 14.1. Недостающие прогоны (все остальные взяты из разделов 9 и 12): GEMMA (-lmm 1), настройки раздела 9.
#  RU (S1, FT_all, 144 образца, поток 1): rint (внутри S1) — схемы A, B, Bn; ang — схема Bn (в 9 для S1 считались только A и B).
#  IT (182, поток 2): ang и rint (внутри IT) — схемы B и Bn (в 9/12 — только A).
#  A = PC (BIC раздела 9: S1 0, IT 0) + родство + происхождение (у IT одно происхождение — не входит); B = A + маркеры статьи (09_covariates_*);
#  Bn = B без Mor534/Mor587/Mor866. Мономорфные и линейно зависимые ковариаты отбрасываются (QR), как в 09_gwas.R.
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"; WORK <- "/work/scratch"; BIN <- "/opt/micromamba/envs/gwas_frost/bin"
RD <- file.path(WORK, "gwas14"); dir.create(RD, showWarnings = FALSE)
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8"); B <- fread(file.path(OUT, "tables/01_blues.csv"))
rint <- function(x) qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
gemma <- function(tag, st, ids, y, C) { d <- file.path(RD, tag); dir.create(d, showWarnings = FALSE); old <- setwd(d); on.exit(setwd(old))
  if (file.exists("output/lmm.assoc.txt")) return(invisible(tag))
  fwrite(data.table(0, ids), "keep.txt", sep = "\t", col.names = FALSE)
  system2(file.path(BIN, "plink2"), c("--bfile", file.path(WORK, "qc", st, "final"), "--allow-extra-chr", "--chr-set", "14", "no-xy", "--keep", "keep.txt", "--maf", "0.05", "--make-bed", "--out", "g"), stdout = FALSE, stderr = FALSE)
  fam <- fread("g.fam", header = FALSE); fam$V6 <- y[match(fam$V2, ids)]; fwrite(fam, "g.fam", sep = " ", col.names = FALSE)
  fwrite(as.data.table(C[match(fam$V2, ids), , drop = FALSE]), "cov.txt", sep = "\t", col.names = FALSE)
  system2(file.path(BIN, "gemma"), c("-bfile", "g", "-gk", "1", "-o", "kin"), stdout = FALSE, stderr = FALSE)
  system2(file.path(BIN, "gemma"), c("-bfile", "g", "-k", "output/kin.cXX.txt", "-c", "cov.txt", "-lmm", "1", "-o", "lmm"), stdout = FALSE, stderr = FALSE)
  invisible(tag) }
covm <- function(ids, origin, COV, scheme) { C <- cbind(int = rep(1, length(ids)))
  if (origin) { o <- factor(OR$origin[match(ids, OR$TP_MO)]); D <- model.matrix(~ o)[, -1, drop = FALSE]; C <- cbind(C, D[, apply(D, 2, sd) > 0, drop = FALSE]) }
  if (scheme != "A") { m <- setdiff(names(COV), "TP_MO"); if (scheme == "Bn") m <- m[!grepl("^Mor", m)]
    M <- as.matrix(COV[match(ids, TP_MO), ..m]); M <- M[, apply(M, 2, function(z) sd(z, na.rm = TRUE)) > 0, drop = FALSE]; C <- cbind(C, M) }
  q <- qr(C); dropped <- setdiff(colnames(C), colnames(C)[q$pivot[seq_len(q$rank)]]); if (length(dropped)) cat("отброшены (ранг):", paste(dropped, collapse = ","), "\n")
  C[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE] }
log <- list()
## RU (S1)
s1 <- readLines(file.path(OUT, "tables/02_samples_stream1.txt")); C1 <- fread(file.path(OUT, "tables/09_covariates_S1.csv"))
d <- B[trait == "FT_all" & scale == "ang" & TP_MO %in% s1, .(TP_MO, y = BLUE)]
for (sc in c("A", "B", "Bn")) { gemma(paste0("RU_rint_", sc), "stream1", d$TP_MO, rint(d$y), covm(d$TP_MO, TRUE, C1, sc)); log[[length(log) + 1]] <- data.table(set = "RU (S1)", trait = "FT_all_rint", scheme = sc, n = nrow(d)) }
gemma("RU_ang_Bn", "stream1", d$TP_MO, d$y, covm(d$TP_MO, TRUE, C1, "Bn")); log[[length(log) + 1]] <- data.table(set = "RU (S1)", trait = "FT_all_ang", scheme = "Bn", n = nrow(d))
## IT
C2 <- fread(file.path(OUT, "tables/09_covariates_S2.csv")); itx <- OR[collection == "italy", TP_MO]
d <- B[trait == "ITFT_m12" & scale == "ang" & TP_MO %in% itx, .(TP_MO, y = BLUE)]
raw_ids <- fread(file.path(WORK, "qc", "stream2", "final.fam"), header = FALSE)$V2; d <- d[TP_MO %in% raw_ids]
for (sc in c("B", "Bn")) for (tr in c("ang", "rint")) { y <- if (tr == "ang") d$y else rint(d$y)
  gemma(sprintf("IT_%s_%s", tr, sc), "stream2", d$TP_MO, y, covm(d$TP_MO, FALSE, C2, sc)); log[[length(log) + 1]] <- data.table(set = "IT", trait = paste0("ITFT_m12_", tr), scheme = sc, n = nrow(d)) }
L <- rbindlist(log); L[, note := "посчитано заново в разделе 14 (GEMMA, настройки раздела 9)"]
fwrite(L, file.path(OUT, "tables/14_new_runs.csv")); print(L)
