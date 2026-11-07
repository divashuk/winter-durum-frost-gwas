# Этап 4. GWAS: GAPIT (MLM, FarmCPU, BLINK) + GEMMA LMM, схемы ковариат A и B.
# Использование: Rscript 04_gwas.R stream1|stream2
#   stream1 = поток 1а: 144 озимых образца статьи (поток 1б невозможен — яровых в статье нет).
#   stream2 = поток 2б: образцы статьи + итальянская коллекция, ковариата «коллекция» (заменяет тип развития:
#             поимённого типа развития итальянских образцов нет; поток 2а «только озимые» невыполним).
# Схема A: PC + родство (+ коллекция в потоке 2). Схема B: A + маркеры durum_frost как фиксированные эффекты.
# Число PC: GAPIT Model.selection (BIC) по схеме A, 0..6 PC; выбранное число используется во всех моделях.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
s    <- commandArgs(TRUE)[1]
OUT  <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
BIN  <- "/opt/micromamba/envs/gwas_frost/bin"
suppressWarnings(suppressMessages(source(file.path(WORK, "soft/gapit_functions.txt"))))
CHR <- c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B")
dir.create(file.path(OUT, "figures/gwas"), showWarnings = FALSE, recursive = TRUE)

## ---------- генотипы ----------
raw <- fread(file.path(WORK, "qc", s, "final.raw"))
GD  <- data.frame(taxa = raw$IID, raw[, -(1:6)], check.names = FALSE); rm(raw)
snp <- names(GD)[-1]
GM  <- data.frame(SNP = snp, Chromosome = match(sub(":.*", "", snp), CHR), Position = as.integer(sapply(strsplit(snp, ":"), `[`, 2)))
o   <- order(GM$Chromosome, GM$Position); GM <- GM[o, ]; GD <- GD[, c(1, o + 1)]
cat(s, ": генотипы", nrow(GD), "x", ncol(GD) - 1, "\n")

## ---------- фенотипы ----------
B  <- fread(file.path(OUT, "tables/01_blues.csv"))
mk <- fread(file.path(OUT, "tables/01_markers_012.csv"))
idm <- fread(file.path(OUT, "tables/01_id_map.csv"), encoding = "UTF-8")
getB <- function(tr, sc) B[trait == tr & scale == sc, .(taxa = TP_MO, y = BLUE)]
if (s == "stream1") {
  TR <- list(FT_all = "FT_all", FT_m12 = "FT_m12", FT_m15 = "FT_m15", FT_m16 = "FT_m16", FT_m17 = "FT_m17")
  COV <- fread(file.path(OUT, "tables/07_covariates_S1.csv"))   # каждый маркер — отдельная ковариата (этап 7)
  SCHEMES <- c("A", "B")
} else {
  TR <- list(FT12j = c("FT_m12", "ITFT_m12"), FTcol = c("FT_all", "ITFT_m12"))
  COV <- fread(file.path(OUT, "tables/07_covariates_S2.csv"))   # Mor* у итальянцев импутированы по GBS (этап 7)
  SCHEMES <- c("A", "B", "Bn")   # Bn — проверка: схема B без Mor534/587/866
}
if (length(commandArgs(TRUE)) > 1) TR <- TR[commandArgs(TRUE)[2]]   # отладка: один признак
meff <- fread(file.path(OUT, sprintf("tables/03_meff_%s.csv", s)))[chr == "ИТОГО", Meff]
pcs_all <- fread(file.path(WORK, "pca", s, "pca.eigenvec")); if ("IID" %in% names(pcs_all)) pcs_all[, 1 := NULL]
setnames(pcs_all, 1, "taxa")

lam <- function(p) median(qchisq(1 - p[is.finite(p) & p > 0], 1)) / qchisq(.5, 1)
allres <- list(); runinfo <- list()

for (tn in names(TR)) for (sc in c("ang", "raw")) {
  Y <- rbindlist(lapply(TR[[tn]], getB, sc = sc))[taxa %in% GD$taxa]
  trait <- paste(tn, sc, sep = "_")
  Y <- Y[!is.na(y)]; setnames(Y, "y", trait)
  n <- nrow(Y)
  # ковариаты
  CVa <- data.table(taxa = Y$taxa)
  if (s == "stream2") CVa[, italy := as.integer(idm$collection[match(taxa, idm$TP_MO)] == "italy")]
  CVb <- copy(CVa)
  for (m in setdiff(names(COV), "TP_MO")) { x <- COV[[m]][match(Y$taxa, COV$TP_MO)]; if (sd(x) > 0) CVb[, (m) := x] }
  setnames(CVb, make.names(names(CVb)))   # «Vrn-A1b-Marq» -> «Vrn.A1b.Marq» (GAPIT)
  CVbn <- CVb[, !grepl("^Mor", names(CVb)), with = FALSE]
  for (scheme in SCHEMES) {
    CV <- switch(scheme, A = CVa, B = CVb, Bn = CVbn)
    rd <- file.path(WORK, "gwas", s, paste(trait, scheme, sep = "_")); dir.create(rd, recursive = TRUE, showWarnings = FALSE)
    setwd(rd)
    if (file.exists("done.rds")) {   # возобновление после сбоя: прогон уже выполнен
      dn <- readRDS("done.rds"); allres <- c(allres, dn$res); runinfo[[length(runinfo) + 1]] <- dn$info
      if (scheme == "A") npc <- dn$info$nPC
      cat("пропуск (готово):", s, trait, scheme, "\n"); setwd(WORK); next
    }
    n_before <- length(allres)
    gd <- GD[GD$taxa %in% Y$taxa, ]; y <- as.data.frame(Y[match(gd$taxa, taxa)])
    cvdf <- if (ncol(CV) > 1) as.data.frame(CV[match(gd$taxa, taxa)]) else NULL
    ## --- выбор числа PC по BIC (только схема A; B использует то же число) ---
    if (scheme == "A") {
      dir.create("bic", showWarnings = FALSE); setwd("bic")
      # BIC только по числу PC (CV = NULL): в X0 GAPIT смешивает PC и CV, порядок неочевиден
      invisible(capture.output(suppressWarnings(GAPIT(Y = y, GD = gd, GM = GM, CV = NULL, PCA.total = 6, model = "MLM",
                     Model.selection = TRUE, kinship.algorithm = "VanRaden", SNP.MAF = 0.05, file.output = TRUE))))
      bf <- list.files(".", "BIC.*Model.Selection.*csv$", recursive = TRUE, full.names = TRUE)
      npc <- 3
      if (length(bf)) {
        bt <- fread(bf[1]); bcol <- grep("BIC", names(bt), value = TRUE)[1]; pcol <- grep("PC", names(bt), value = TRUE)[1]
        npc <- as.integer(bt[[pcol]][which.max(bt[[bcol]])])
        fwrite(bt, file.path(OUT, sprintf("tables/04_bic_%s_%s_%s.csv", s, trait, "A")))
      }
      setwd(rd)
    }
    ## --- GAPIT: MLM, FarmCPU, BLINK ---
    t0 <- Sys.time()
    invisible(capture.output(suppressWarnings(GAPIT(Y = y, GD = gd, GM = GM, CV = cvdf, PCA.total = npc,
                     model = c("MLM", "FarmCPU", "BLINK"), kinship.algorithm = "VanRaden", SNP.MAF = 0.05, file.output = TRUE))))
    for (mod in c("MLM", "FarmCPU", "BLINK")) {
      f <- list.files(".", sprintf("GWAS_Results\\.%s\\..*csv$", mod), full.names = TRUE)
      if (!length(f)) { cat("нет результата", mod, trait, scheme, "\n"); next }
      r <- fread(f[1])
      pc <- grep("^P.value$", names(r), value = TRUE); ec <- grep("^Effect$|effect", names(r), value = TRUE, ignore.case = TRUE)[1]
      allres[[length(allres) + 1]] <- data.table(stream = s, trait = trait, scheme = scheme, model = mod, n = n, nPC = npc,
        SNP = r$SNP, chr = CHR[r$Chr], pos = r$Pos, MAF = r$MAF, p = r[[pc]], effect = if (!is.na(ec)) r[[ec]] else NA_real_)
    }
    ## --- GEMMA LMM (Wald), ковариаты: 1 + PC(npc) + CV ---
    gdir <- file.path(rd, "gemma"); dir.create(gdir, showWarnings = FALSE); setwd(gdir)
    fwrite(data.table(0, gd$taxa), "keep.txt", sep = "\t", col.names = FALSE)   # FID = 0 в final.fam
    system2(file.path(BIN, "plink2"), c("--bfile", file.path(WORK, "qc", s, "final"), "--allow-extra-chr", "--chr-set", "14", "no-xy",
            "--keep", "keep.txt", "--maf", "0.05", "--make-bed", "--out", "g"), stdout = FALSE, stderr = FALSE)
    fam <- fread("g.fam", header = FALSE); fam$V6 <- y[[trait]][match(fam$V2, y$taxa)]
    fwrite(fam, "g.fam", sep = " ", col.names = FALSE)
    cov <- data.table(int = 1)[rep(1, nrow(fam))]
    if (npc > 0) cov <- cbind(cov, pcs_all[match(fam$V2, taxa), paste0("PC", seq_len(npc)), with = FALSE])
    if (ncol(CV) > 1) cov <- cbind(cov, CV[match(fam$V2, taxa), -1])
    fwrite(cov, "cov.txt", sep = "\t", col.names = FALSE)
    system2(file.path(BIN, "gemma"), c("-bfile", "g", "-gk", "1", "-o", "kin"), stdout = FALSE, stderr = FALSE)
    system2(file.path(BIN, "gemma"), c("-bfile", "g", "-k", "output/kin.cXX.txt", "-c", "cov.txt", "-lmm", "1", "-o", "lmm"), stdout = FALSE, stderr = FALSE)
    if (file.exists("output/lmm.assoc.txt")) {
      r <- fread("output/lmm.assoc.txt")
      allres[[length(allres) + 1]] <- data.table(stream = s, trait = trait, scheme = scheme, model = "GEMMA", n = n, nPC = npc,
        SNP = r$rs, chr = sub(":.*", "", r$rs), pos = r$ps, MAF = r$af, p = r$p_wald, effect = r$beta)
    }
    runinfo[[length(runinfo) + 1]] <- data.table(stream = s, trait = trait, scheme = scheme, n = n, nPC = npc,
      covariates = paste(names(CV)[-1], collapse = "+"), sec = round(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    saveRDS(list(res = allres[(n_before + 1):length(allres)], info = runinfo[[length(runinfo)]]), file.path(rd, "done.rds"))
    cat(sprintf("%s %s %s: n=%d, PC=%d, CV=%s\n", s, trait, scheme, n, npc, paste(names(CV)[-1], collapse = "+")))
    setwd(WORK)
  }
}
R <- rbindlist(allres)
R[MAF > 0.5, MAF := 1 - MAF]
R[, `:=`(M = .N, q_BH = p.adjust(p, "BH"), lambda = lam(p)), by = .(trait, scheme, model)]
R[, `:=`(thr_bonf = 0.05 / M, thr_meff = 0.05 / meff)]
R[, `:=`(sig_bonf = p < thr_bonf, sig_meff = p < thr_meff, sig_fdr = q_BH < 0.05, suggestive = -log10(p) >= 4)]
fwrite(R, file.path(WORK, "gwas", sprintf("all_results_%s.csv.gz", s)))
fwrite(R[p < 1e-3], file.path(OUT, sprintf("tables/04_gwas_hits_p1e-3_%s.csv", s)))
L <- unique(R[, .(stream, trait, scheme, model, n, nPC, M, Meff = meff, lambda,
                  n_bonf = sum(sig_bonf), n_meff = sum(sig_meff), n_fdr = sum(sig_fdr), n_sugg = sum(suggestive), min_p = min(p)),
              by = .(trait, scheme, model)])
L <- L[, .SD[1], by = .(trait, scheme, model)]
L <- merge(L, rbindlist(runinfo)[, .(trait, scheme, covariates)], by = c("trait", "scheme"))
fwrite(L, file.path(OUT, sprintf("tables/04_gwas_run_summary_%s.csv", s)))
print(L[, .(trait, scheme, model, n, nPC, lambda = round(lambda, 2), n_bonf, n_meff, n_fdr, n_sugg, minlogp = round(-log10(min_p), 2))])

## ---------- графики: Manhattan (4 модели) + QQ на сочетание признак×схема ----------
cum <- GM; cum$chr <- CHR[cum$Chromosome]
off <- tapply(cum$Position, cum$Chromosome, max); off <- c(0, cumsum(as.numeric(off)))[1:14]; names(off) <- CHR
for (k in unique(R[, paste(trait, scheme)])) {
  d <- R[paste(trait, scheme) == k]; d[, x := pos + off[chr]]; d[, lp := -log10(p)]
  d[lp > 30, lp := 30]
  ax <- data.table(chr = CHR, mid = off + tapply(cum$Position, cum$Chromosome, max) / 2)
  thr <- unique(d[, .(model, b = -log10(thr_bonf), m = -log10(thr_meff))])
  pm <- ggplot(d, aes(x, lp, colour = factor(match(chr, CHR) %% 2))) + geom_point(size = .5) +
    geom_hline(data = thr, aes(yintercept = b), colour = "red", linetype = 2) +
    geom_hline(data = thr, aes(yintercept = m), colour = "orange", linetype = 2) +
    geom_hline(yintercept = 4, colour = "grey50", linetype = 3) +
    facet_wrap(~model, ncol = 1) + scale_colour_manual(values = c("grey20", "steelblue"), guide = "none") +
    scale_x_continuous(breaks = ax$mid, labels = ax$chr) + labs(x = NULL, y = "-log10 p",
      title = sprintf("%s — %s, схема %s (n=%d). Красная: Бонферрони, оранжевая: Meff, серая: 1e-4", s, sub(" .*", "", k), sub(".* ", "", k), d$n[1])) +
    theme_bw(base_size = 9)
  q <- d[, .(obs = sort(-log10(p), decreasing = TRUE), exp = -log10(ppoints(.N))), by = model]
  lab <- unique(d[, .(model, lambda)])
  pq <- ggplot(q, aes(exp, pmin(obs, 30), colour = model)) + geom_point(size = .6) + geom_abline(linetype = 2) +
    geom_text(data = lab, aes(x = 0.5, y = 28 - 2.5 * (as.integer(factor(model)) - 1), label = sprintf("%s λ=%.2f", model, lambda)),
              hjust = 0, show.legend = FALSE, size = 3) + labs(x = "ожидаемое -log10 p", y = "наблюдаемое -log10 p") + theme_bw(base_size = 9)
  g <- gridExtra::arrangeGrob(pm, pq, ncol = 2, widths = c(3, 1.2))
  ggsave(file.path(OUT, sprintf("figures/gwas/%s_%s.png", s, gsub(" ", "_", k))), g, width = 15, height = 9, dpi = 130)
}
writeLines(capture.output(sessionInfo()), file.path(OUT, sprintf("sessionInfo_04_gwas_%s.txt", s)))
