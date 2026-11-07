# Этап 10.7 (часть 3). Таблица KASP: последовательность [REF/ALT] с фланками 100 п.н., соседние варианты в ±50 п.н.,
# число хитов фланка по геному (длина >= 150 п.н.; идентичность >= 0.90 — «хит»; отдельно — лучший чужой хит, обычно гомеолог).
# Годен (passes): REF во фланке совпадает с аллелем SNP, нет других вариантов в ±50 п.н., N нет, ровно 1 хит с идентичностью >= 0.90 (сам себя).
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"
K <- "/work/scratch/kasp"
S <- fread(file.path(K, "kasp_candidates_raw.tsv"))
fa <- readLines(file.path(K, "flanks.fa")); hi <- grep("^>", fa)
seqs <- data.table(reg = sub("^>", "", fa[hi]), seq = sapply(seq_along(hi), function(i) paste(fa[(hi[i] + 1):(if (i < length(hi)) hi[i + 1] - 1 else length(fa))], collapse = "")))
S[, reg := sprintf("%s:%d-%d", chr, pos - 100, pos + 100)]; S <- merge(S, seqs, by = "reg", all.x = TRUE)
S[, `:=`(ref = sub(":.*", "", alleles), alt = sub(".*:", "", alleles))]
S[, ref_in_genome := toupper(substr(seq, 101, 101))]
S[, kasp_sequence := paste0(substr(seq, 1, 100), "[", ref, "/", alt, "]", substr(seq, 102, 201))]
S[, n_N := nchar(gsub("[^Nn]", "", seq))]
V <- fread(file.path(K, "all_variants_named.tsv"), header = FALSE, col.names = c("chr", "pos")); setkey(V, chr, pos)
S[, other_var_50bp := mapply(function(c, p) V[chr == c & pos >= p - 50 & pos <= p + 50 & pos != p, paste(pos - p, collapse = ",")], chr, pos)]
S[, other_var_100bp_n := mapply(function(c, p) V[chr == c & pos >= p - 100 & pos <= p + 100 & pos != p, .N], chr, pos)]
P <- fread(file.path(K, "flank_hits.paf"), header = FALSE, fill = TRUE, col.names = c("q", "qlen", "qs", "qe", "strand", "t", "tlen", "ts", "te", "nm", "alen", "mapq"))
P <- P[alen >= 150][, ident := nm / alen]
P[, `:=`(qchr = sub(":.*", "", q), qstart = as.numeric(sub(".*:([0-9]+)-.*", "\\1", q)))]
P[, self := t == qchr & abs(ts + 1 - qstart) <= 5]
H <- P[, .(n_hits_id90 = sum(ident >= 0.90), self_found = any(self), best_other_ident = round(max(c(0, ident[!self])), 3),
           best_other_locus = if (any(!self)) paste0(t[!self][which.max(ident[!self])], ":", ts[!self][which.max(ident[!self])]) else ""), by = .(reg = q)]
S <- merge(S, H, by = "reg", all.x = TRUE); S[is.na(n_hits_id90), n_hits_id90 := 0L]
S[, passes := ref_in_genome == ref & other_var_50bp == "" & n_N == 0 & n_hits_id90 == 1 & self_found %in% TRUE]
S[, fail_reason := paste0(fifelse(ref_in_genome != ref, "REF не совпадает с геномом; ", ""), fifelse(other_var_50bp != "", "другой SNP в ±50 п.н.; ", ""),
                          fifelse(n_N > 0, "N во фланке; ", ""), fifelse(n_hits_id90 != 1, sprintf("хитов с идентичностью >= 0.90: %d; ", n_hits_id90), ""))]
L <- fread(file.path(OUT, "tables/10_loci.csv"))[, .(locus_id, locus_name = name, category)]
S <- merge(S, L, by = "locus_id")
res <- S[order(locus_id, -passes, -best_mlp), .(locus_id, locus_name, category, snp = sprintf("%s:%d", chr, pos), chr, pos, ref, alt, reason, best_mlp, best_at,
          gemma_effect, gemma_effect_allele, gemma_at, kasp_sequence, other_var_50bp, other_var_100bp_n, n_hits_id90, best_other_ident, best_other_locus,
          passes, fail_reason)]
fwrite(res, file.path(OUT, "tables/10_kasp_candidates.csv"))
print(res[, .(n = .N, passes = sum(passes)), by = .(locus_id, locus_name)])
print(table(res$fail_reason))
