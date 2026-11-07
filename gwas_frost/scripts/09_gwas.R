# Этап 9.1–9.5. Полный пересчёт GWAS с происхождением как обязательной ковариатой.
# Использование: Rscript 09_gwas.R S1|S2|RU|IT
#  S1 (144, генотипы потока 1): FT_all_ang (A, B; MLM/FarmCPU/BLINK/GEMMA), LT50 (A, B; MLM, GEMMA)
#  S2 (326, генотипы потока 2): FT12j и FTcol — ang (исходный), rint (ранговое обратное нормальное внутри коллекции; основной),
#     bin (выжило хотя бы одно растение, 0/1, линейная модель); схемы A, B, Bn; 4 модели
#  RU (91 с данными −12 °C): FT_m12_ang — схема A (PC + родство + происхождение), MLM и GEMMA (для сравнения внутри российской коллекции; FT_all — это S1 A)
#  IT (182, только итальянские; генотипы потока 2): ITFT_m12_ang — PC + родство (поимённого происхождения нет), MLM и GEMMA
# Схемы: A = PC + родство + происхождение (дамми; в S2 итальянская коллекция — одна из категорий, «коллекция» отдельно не входит);
#        B = A + все маркеры статьи (09_covariates_*, каждый отдельно, VIF > 10 не трогаются); Bn (S2) = B без Mor.
# Число PC: BIC GAPIT (0..6, только PC). Для GEMMA PC считаются тем же числом внутри анализируемой подвыборки (prcomp).
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
set_ <- commandArgs(TRUE)[1]
OUT  <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
BIN  <- "/opt/micromamba/envs/gwas_frost/bin"
suppressWarnings(suppressMessages(source(file.path(WORK, "soft/gapit_functions.txt"))))
CHR <- c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B")
dir.create(file.path(OUT, "figures/gwas09"), showWarnings = FALSE, recursive = TRUE)
gstream <- if (set_ %in% c("S1", "RU")) "stream1" else "stream2"

raw <- fread(file.path(WORK, "qc", gstream, "final.raw"))
GD  <- data.frame(taxa = raw$IID, raw[, -(1:6)], check.names = FALSE); rm(raw)
snp <- names(GD)[-1]
GM  <- data.frame(SNP = snp, Chromosome = match(sub(":.*", "", snp), CHR), Position = as.integer(sapply(strsplit(snp, ":"), `[`, 2)))
o   <- order(GM$Chromosome, GM$Position); GM <- GM[o, ]; GD <- GD[, c(1, o + 1)]

B   <- fread(file.path(OUT, "tables/01_blues.csv"))
OR  <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8")
LT  <- fread(file.path(OUT, "tables/09_lt50.csv"))
S2T <- fread(file.path(OUT, "tables/09_s2_traits.csv"))
getB <- function(tr) B[trait == tr & scale == "ang", .(taxa = TP_MO, y = BLUE)]
ALL4 <- c("MLM", "FarmCPU", "BLINK"); MLMONLY <- "MLM"
if (set_ == "S1") {
  COV <- fread(file.path(OUT, "tables/09_covariates_S1.csv")); SCH <- c("A", "B")
  # LT50 ограничен диапазоном −20…−10 °C (оценки вне — экстраполяция при 3–5 температурах)
  TR <- list(FT_all_ang = list(getB("FT_all"), ALL4),
             LT50 = list(LT[!is.na(LT50) & LT50 >= -20 & LT50 <= -10, .(taxa = TP_MO, y = LT50)], MLMONLY))
} else if (set_ == "S2") {
  COV <- fread(file.path(OUT, "tables/09_covariates_S2.csv")); SCH <- c("A", "B", "Bn")
  TR <- list()
  for (t in c("FT12j", "FTcol")) for (v in c("ang", "rint", "bin")) {
    d <- S2T[trait == t, .(taxa = TP_MO, y = switch(v, ang = y, rint = rint, bin = as.numeric(bin)))]
    TR[[paste(t, v, sep = "_")]] <- list(d, ALL4) }
} else if (set_ == "RU") {
  COV <- NULL; SCH <- "A"; TR <- list(FT_m12_ang = list(getB("FT_m12"), MLMONLY))
} else {
  COV <- NULL; SCH <- "A"; TR <- list(ITFT_m12_ang = list(getB("ITFT_m12"), MLMONLY))
}
keep_taxa <- switch(set_, S1 = readLines(file.path(OUT, "tables/02_samples_stream1.txt")), RU = readLines(file.path(OUT, "tables/02_samples_stream1.txt")),
                    S2 = readLines(file.path(OUT, "tables/02_samples_stream2.txt")), IT = OR[collection == "italy", TP_MO])
meff <- fread(file.path(OUT, sprintf("tables/03_meff_%s.csv", gstream)))[chr == "ИТОГО", Meff]
lam <- function(p) median(qchisq(1 - p[is.finite(p) & p > 0], 1)) / qchisq(.5, 1)
allres <- list(); runinfo <- list()

for (trait in names(TR)) {
  Y <- TR[[trait]][[1]][taxa %in% GD$taxa & taxa %in% keep_taxa & !is.na(y)]; setnames(Y, "y", trait)
  gmodels <- TR[[trait]][[2]]; n <- nrow(Y)
  orig <- OR$origin[match(Y$taxa, OR$TP_MO)]
  CVa <- data.table(taxa = Y$taxa)
  if (set_ %in% c("S1", "S2", "RU") && length(unique(orig)) > 1) { Dm <- model.matrix(~ factor(orig))[, -1, drop = FALSE]   # LT50: все 47 — НЦЗ, ковариата не нужна
    colnames(Dm) <- paste0("orig", seq_len(ncol(Dm))); CVa <- cbind(CVa, Dm[, apply(Dm, 2, sd) > 0, drop = FALSE]) }
  CVb <- copy(CVa)
  if (!is.null(COV)) for (m in setdiff(names(COV), "TP_MO")) { x <- COV[[m]][match(Y$taxa, COV$TP_MO)]; if (sd(x) > 0) CVb[, (m) := x] }
  setnames(CVb, make.names(names(CVb))); CVbn <- CVb[, !grepl("^Mor", names(CVb)), with = FALSE]
  for (scheme in SCH) {
    if (trait == "LT50" && scheme == "Bn") next
    CV <- switch(scheme, A = CVa, B = CVb, Bn = CVbn)
    if (ncol(CV) > 1) {   # ковариаты, дающие вырожденную матрицу (малые подвыборки, напр. LT50 n = 47), отбрасываются
      Xc <- cbind(1, as.matrix(CV[, -1])); qx <- qr(Xc); keepc <- sort(qx$pivot[seq_len(qx$rank)]) - 1; keepc <- keepc[keepc > 0]
      dropped <- setdiff(names(CV)[-1], names(CV)[-1][keepc]); if (length(dropped)) cat("отброшены (ранг):", trait, scheme, paste(dropped, collapse = ","), "
")
      CV <- CV[, c(1, keepc + 1), with = FALSE] }
    rd <- file.path(WORK, "gwas09", set_, paste(trait, scheme, sep = "_")); dir.create(rd, recursive = TRUE, showWarnings = FALSE); setwd(rd)
    if (file.exists("done.rds")) { dn <- readRDS("done.rds"); allres <- c(allres, dn$res); runinfo[[length(runinfo) + 1]] <- dn$info
      if (scheme == "A") npc <- dn$info$nPC; cat("пропуск (готово):", set_, trait, scheme, "\n"); setwd(WORK); next }
    n_before <- length(allres); t0 <- Sys.time()
    gd <- GD[GD$taxa %in% Y$taxa, ]; y <- as.data.frame(Y[match(gd$taxa, taxa)])
    cvdf <- if (ncol(CV) > 1) as.data.frame(CV[match(gd$taxa, taxa)]) else NULL
    if (scheme == "A") {
      dir.create("bic", showWarnings = FALSE); setwd("bic")
      invisible(capture.output(suppressWarnings(GAPIT(Y = y, GD = gd, GM = GM, CV = NULL, PCA.total = 6, model = "MLM",
                     Model.selection = TRUE, kinship.algorithm = "VanRaden", SNP.MAF = 0.05, file.output = TRUE))))
      bf <- list.files(".", "BIC.*Model.Selection.*csv$", recursive = TRUE, full.names = TRUE); npc <- 3
      if (length(bf)) { bt <- fread(bf[1]); npc <- as.integer(bt[[grep("PC", names(bt), value = TRUE)[1]]][which.max(bt[[grep("BIC", names(bt), value = TRUE)[1]]])])
        fwrite(bt, file.path(OUT, sprintf("tables/09_bic_%s_%s.csv", set_, trait))) }
      setwd(rd)
    }
    invisible(capture.output(suppressWarnings(GAPIT(Y = y, GD = gd, GM = GM, CV = cvdf, PCA.total = npc, model = gmodels,
                     kinship.algorithm = "VanRaden", SNP.MAF = 0.05, file.output = TRUE))))
    for (mod in gmodels) {
      f <- list.files(".", sprintf("GWAS_Results\\.%s\\..*csv$", mod), full.names = TRUE); if (!length(f)) next
      r <- fread(f[1])
      allres[[length(allres) + 1]] <- data.table(set = set_, trait = trait, scheme = scheme, model = mod, n = n, nPC = npc,
        SNP = r$SNP, chr = CHR[r$Chr], pos = r$Pos, MAF = r$MAF, p = r$P.value, effect = r$Effect, effect_allele = sub(".*_", "", r$SNP))
    }
    ## GEMMA: 1 + PC (внутри подвыборки) + CV
    gdir <- file.path(rd, "gemma"); dir.create(gdir, showWarnings = FALSE); setwd(gdir)
    fwrite(data.table(0, gd$taxa), "keep.txt", sep = "\t", col.names = FALSE)
    system2(file.path(BIN, "plink2"), c("--bfile", file.path(WORK, "qc", gstream, "final"), "--allow-extra-chr", "--chr-set", "14", "no-xy",
            "--keep", "keep.txt", "--maf", "0.05", "--make-bed", "--out", "g"), stdout = FALSE, stderr = FALSE)
    fam <- fread("g.fam", header = FALSE); fam$V6 <- y[[trait]][match(fam$V2, y$taxa)]; fwrite(fam, "g.fam", sep = " ", col.names = FALSE)
    cov <- data.table(int = rep(1, nrow(fam)))
    if (npc > 0) { Xg <- as.matrix(gd[match(fam$V2, gd$taxa), -1]); Xg <- Xg[, apply(Xg, 2, sd) > 0]
      cov <- cbind(cov, prcomp(Xg, center = TRUE, scale. = FALSE, rank. = npc)$x) }
    if (ncol(CV) > 1) cov <- cbind(cov, CV[match(fam$V2, taxa), -1])
    fwrite(cov, "cov.txt", sep = "\t", col.names = FALSE)
    system2(file.path(BIN, "gemma"), c("-bfile", "g", "-gk", "1", "-o", "kin"), stdout = FALSE, stderr = FALSE)
    system2(file.path(BIN, "gemma"), c("-bfile", "g", "-k", "output/kin.cXX.txt", "-c", "cov.txt", "-lmm", "1", "-o", "lmm"), stdout = FALSE, stderr = FALSE)
    if (file.exists("output/lmm.assoc.txt")) { r <- fread("output/lmm.assoc.txt")
      allres[[length(allres) + 1]] <- data.table(set = set_, trait = trait, scheme = scheme, model = "GEMMA", n = n, nPC = npc,
        SNP = r$rs, chr = sub(":.*", "", r$rs), pos = r$ps, MAF = r$af, p = r$p_wald, effect = r$beta, effect_allele = r$allele1) }
    runinfo[[length(runinfo) + 1]] <- data.table(set = set_, trait = trait, scheme = scheme, n = n, nPC = npc,
      covariates = paste(names(CV)[-1], collapse = "+"), sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    saveRDS(list(res = allres[(n_before + 1):length(allres)], info = runinfo[[length(runinfo)]]), file.path(rd, "done.rds"))
    cat(sprintf("%s %s %s: n=%d, PC=%d, CV=%s\n", set_, trait, scheme, n, npc, paste(names(CV)[-1], collapse = "+")))
    setwd(WORK)
  }
}
R <- rbindlist(allres, fill = TRUE); R[MAF > 0.5, MAF := 1 - MAF]
R[, `:=`(M = .N, q_BH = p.adjust(p, "BH"), lambda = lam(p)), by = .(trait, scheme, model)]
R[, `:=`(thr_bonf = 0.05 / M, thr_meff = 0.05 / meff)]
R[, `:=`(sig_bonf = p < thr_bonf, sig_meff = p < thr_meff, sig_fdr = q_BH < 0.05, suggestive = -log10(p) >= 4)]
fwrite(R, file.path(WORK, "gwas09", sprintf("all_results_09_%s.csv.gz", set_)))
L <- R[, .(n = n[1], nPC = nPC[1], M = M[1], Meff = meff, lambda = round(lambda[1], 3), n_bonf = sum(sig_bonf), n_meff = sum(sig_meff),
           n_fdr = sum(sig_fdr), n_sugg = sum(suggestive), min_p = min(p), top_snp = SNP[which.min(p)]), by = .(set, trait, scheme, model)]
L <- merge(L, rbindlist(runinfo)[, .(trait, scheme, covariates)], by = c("trait", "scheme"))
fwrite(L, file.path(OUT, sprintf("tables/09_gwas_run_summary_%s.csv", set_)))
print(L[, .(trait, scheme, model, n, nPC, lambda, n_bonf, n_meff, n_fdr, n_sugg, mlp = round(-log10(min_p), 2))])
## графики Manhattan + QQ
cum <- GM; off <- tapply(cum$Position, cum$Chromosome, max); off <- c(0, cumsum(as.numeric(off)))[1:14]; names(off) <- CHR
for (k in unique(R[, paste(trait, scheme)])) {
  d <- R[paste(trait, scheme) == k]; d[, x := pos + off[chr]]; d[, lp := pmin(-log10(p), 30)]
  ax <- data.table(chr = CHR, mid = off + tapply(cum$Position, cum$Chromosome, max) / 2)
  thr <- unique(d[, .(model, b = -log10(thr_bonf), m = -log10(thr_meff))])
  pm <- ggplot(d, aes(x, lp, colour = factor(match(chr, CHR) %% 2))) + geom_point(size = .5) +
    geom_hline(data = thr, aes(yintercept = b), colour = "red", linetype = 2) + geom_hline(data = thr, aes(yintercept = m), colour = "orange", linetype = 2) +
    geom_hline(yintercept = 4, colour = "grey50", linetype = 3) + facet_wrap(~model, ncol = 1) +
    scale_colour_manual(values = c("grey20", "steelblue"), guide = "none") + scale_x_continuous(breaks = ax$mid, labels = ax$chr) +
    labs(x = NULL, y = "-log10 p", title = sprintf("09 %s — %s (n=%d). Красная: Бонферрони, оранжевая: Meff, серая: 1e-4", set_, k, d$n[1])) + theme_bw(base_size = 9)
  q <- d[, .(obs = sort(-log10(p), decreasing = TRUE), exp = -log10(ppoints(.N))), by = model]
  lab <- unique(d[, .(model, lambda)])
  pq <- ggplot(q, aes(exp, pmin(obs, 30), colour = model)) + geom_point(size = .6) + geom_abline(linetype = 2) +
    geom_text(data = lab, aes(x = .5, y = 28 - 2.5 * (as.integer(factor(model)) - 1), label = sprintf("%s λ=%.2f", model, lambda)), hjust = 0, show.legend = FALSE, size = 3) +
    labs(x = "ожидаемое -log10 p", y = "наблюдаемое -log10 p") + theme_bw(base_size = 9)
  ggsave(file.path(OUT, sprintf("figures/gwas09/09_%s_%s.png", set_, gsub(" ", "_", k))), gridExtra::arrangeGrob(pm, pq, ncol = 2, widths = c(3, 1.2)),
         width = 15, height = if (length(unique(d$model)) > 2) 9 else 6, dpi = 120)
}
writeLines(capture.output(sessionInfo()), file.path(OUT, sprintf("sessionInfo_09_gwas_%s.txt", set_)))
