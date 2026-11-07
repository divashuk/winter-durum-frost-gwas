# Этап 10.7 (часть 1). Отбор SNP для KASP и фланки ±100 п.н. из Langdon.
#  Кандидаты на локус: ведущий SNP, до 5 лучших SNP MLM/GEMMA этапа 9 в узком интервале, все SNP GBS внутри генов-кандидатов (10_candidate_genes.csv).
#  Для каждого: лучший −log10 p MLM/GEMMA (наборы S1, S2, RU, IT), набор/признак/схема, эффект GEMMA и эффектный аллель.
#  Соседние варианты: все SNP обоих исходных VCF (коллекция и F2), до КК — для правила «нет другого SNP GBS в ±50 п.н.».
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"
WORK <- "/work/scratch"; K <- file.path(WORK, "kasp"); dir.create(K, showWarnings = FALSE)
L <- fread(file.path(OUT, "tables/10_loci.csv")); CG <- fread(file.path(OUT, "tables/10_candidate_genes.csv"))
R9 <- rbindlist(lapply(c("S1", "S2", "RU", "IT"), function(s) fread(file.path(WORK, "gwas09", sprintf("all_results_09_%s.csv.gz", s)),
                select = c("set", "trait", "scheme", "model", "SNP", "chr", "pos", "p", "effect", "effect_allele"))))[model %in% c("MLM", "GEMMA")]
R9[, alleles := sub("_[^_]*$", "", sub("^[^:]*:[0-9]+:", "", SNP))]
best <- R9[order(p)][, .SD[1], by = .(chr, pos)][, .(chr, pos, alleles, best_mlp = round(-log10(p), 2), best_at = paste(set, trait, scheme, model, sep = "/"))]
gem <- R9[model == "GEMMA"][order(p)][, .SD[1], by = .(chr, pos)][, .(chr, pos, gemma_effect = signif(effect, 3), gemma_effect_allele = effect_allele, gemma_at = paste(set, trait, scheme, sep = "/"))]
best <- merge(best, gem, by = c("chr", "pos"), all.x = TRUE)
sel <- list()
for (r in seq_len(nrow(L))) { x <- L[r]
  b <- best[chr == x$chr & pos >= x$narrow_start & pos <= x$narrow_end]
  ld <- b[which.min(abs(pos - x$lead_pos))][, reason := "ведущий SNP (ближайший к позиции)"]
  tp <- b[order(-best_mlp)][1:min(5, .N)][, reason := "топ-5 MLM/GEMMA в узком интервале"]
  cg <- CG[locus_id == x$locus_id & nzchar(gbs_snps_in_gene)]
  gs <- if (nrow(cg)) rbindlist(lapply(seq_len(nrow(cg)), function(i) { p <- as.numeric(sub("\\(.*", "", sub("^[^:]+:", "", strsplit(cg$gbs_snps_in_gene[i], ";")[[1]])))
          best[chr == x$chr & pos %in% p][, reason := paste("в гене-кандидате", cg$svevo_gene[i])] })) else NULL
  s <- rbindlist(list(ld, tp, gs), fill = TRUE)
  if (nrow(s)) sel[[r]] <- s[, .(reason = paste(unique(reason), collapse = "; ")), by = setdiff(names(s), "reason")][, locus_id := x$locus_id] }
S <- rbindlist(sel, fill = TRUE)[!is.na(pos)]   # в узком окне Fr-B2 (L22) SNP GBS нет
S[, id := sprintf("%s_%s_%d", locus_id, chr, pos)]
fwrite(S, file.path(K, "kasp_candidates_raw.tsv"), sep = "\t")
writeLines(sprintf("%s:%d-%d", S$chr, S$pos - 100, S$pos + 100), file.path(K, "regions.txt"))
cat("SNP-кандидатов KASP:", nrow(S), "; по локусам:\n"); print(S[, .N, by = locus_id])
