# Этап 10.2 (R-часть 10_genes.sh): гены Svevo в широких интервалах; узкий интервал; расстояние до ведущего SNP;
# SNP GBS (после КК, поток 1 ∪ поток 2) внутри гена и их лучший −log10 p MLM/GEMMA этапа 9 (наборы S1, S2, RU, IT).
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"
WORK <- "/work/scratch"; G <- file.path(WORK, "genes10")
L <- fread(file.path(OUT, "tables/10_loci.csv"))
p <- fread(file.path(G, "hits.paf"), header = FALSE, fill = TRUE, select = 1:12,
           col.names = c("tx", "qlen", "qs", "qe", "strand", "target", "tlen", "ts", "te", "nmatch", "alnlen", "mapq"))
p[, `:=`(identity = nmatch / alnlen, qcov = (qe - qs) / qlen)]
p <- p[identity >= 0.98 & qcov >= 0.90]
p[, gene := sub("\\.[0-9]+$", "", tx)]
p <- p[order(-nmatch)][, .SD[1], by = .(gene, target)]
p[, `:=`(locus_id = sub("\\|.*", "", target), region = sub(".*\\|", "", target))]
p[, int_start := as.numeric(sub(".*:([0-9]+)-.*", "\\1", region))]
p[, `:=`(langdon_start = int_start + ts - 1, langdon_end = int_start + te - 1)]
h <- fread(file.path(G, "cdna_headers.txt"), header = FALSE, sep = "\n", col.names = "line")
h[, tx := sub(" .*", "", line)][, svevo_pos := sub(".*chromosome:Svevo.v1:([^ ]+).*", "\\1", line)][, biotype := sub(".*gene_biotype:([^ ]+).*", "\\1", line)]
p <- merge(p, h[, .(tx, svevo_pos, biotype)], by = "tx", all.x = TRUE)
p <- merge(p, L[, .(locus_id, name, chr, category, lead_pos, narrow_start, narrow_end, qtl_start, qtl_end)], by = "locus_id")
p[, svevo_chr := sub(":.*", "", svevo_pos)][, svevo_start := as.numeric(sub("^[^:]+:([0-9]+):.*", "\\1", svevo_pos))]
p[, dist_to_lead := fifelse(lead_pos >= langdon_start & lead_pos <= langdon_end, 0, pmin(abs(langdon_start - lead_pos), abs(langdon_end - lead_pos)))]
p[, in_narrow := langdon_end >= narrow_start & langdon_start <= narrow_end]
## SNP внутри гена
bim <- unique(rbindlist(lapply(c("stream1", "stream2"), function(s) fread(file.path(WORK, "qc", s, "final.bim"), header = FALSE)[, .(chr = V1, pos = V4, snp = V2)])), by = c("chr", "pos"))
bim[, chr := as.character(chr)]
R9 <- rbindlist(lapply(c("S1", "S2", "RU", "IT"), function(s) fread(file.path(WORK, "gwas09", sprintf("all_results_09_%s.csv.gz", s)),
                select = c("set", "trait", "scheme", "model", "chr", "pos", "p"))))
R9 <- R9[model %in% c("MLM", "GEMMA"), .(best_mlp = round(-log10(min(p)), 2), best_at = paste(set, trait, scheme, model, sep = "/")[which.min(p)]), by = .(chr, pos)]
bim <- merge(bim, R9, by = c("chr", "pos"), all.x = TRUE)
setkey(bim, chr, pos)
snp_in <- function(ch, s, e) { b <- bim[chr == ch & pos >= s & pos <= e]
  if (!nrow(b)) return(list("", NA_real_, ""))
  list(paste(sprintf("%s:%d(%s)", b$chr, b$pos, ifelse(is.na(b$best_mlp), "NA", b$best_mlp)), collapse = ";"), max(b$best_mlp, na.rm = TRUE) |> suppressWarnings(), b$best_at[which.max(b$best_mlp)] %||% "") }
`%||%` <- function(a, b) if (length(a) == 0 || is.na(a)) b else a
p[, c("gbs_snps_in_gene", "gbs_best_mlp_in_gene", "gbs_best_at") := snp_in(chr, langdon_start, langdon_end), by = .(locus_id, gene)]
p[!is.finite(gbs_best_mlp_in_gene), gbs_best_mlp_in_gene := NA]
res <- p[order(locus_id, langdon_start), .(locus_id, locus_name = name, category, chr, langdon_start, langdon_end, strand, svevo_gene = gene, svevo_transcript = tx,
           svevo_chr, svevo_start, biotype, identity = round(identity, 4), coverage = round(qcov, 3), lead_pos, dist_to_lead, in_narrow,
           gbs_snps_in_gene, gbs_best_mlp_in_gene, gbs_best_at)]
fwrite(res, file.path(G, "genes_wide.csv"))
print(res[, .(n_wide = .N, n_narrow = sum(in_narrow), n_svevo_other_chr = sum(svevo_chr != chr), n_with_snp = sum(nzchar(gbs_snps_in_gene))), by = .(locus_id, locus_name)])
cat("уникальных генов:", uniqueN(res$svevo_gene), "\n")
