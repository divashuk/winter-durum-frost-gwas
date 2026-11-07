# Этап 14.4 (часть 3, R-часть 14_genes.sh). Фильтр переноса и положение генов.
#  Хит: идентичность >= 0.98, покрытие транскрипта >= 0.90, лучший хит на ген в интервале (как в разделах 5 и 10).
#  Коллинеарность: смещение (Langdon − Svevo) гена отличается от скользящей медианы смещения соседних генов (31 ген) не более чем на 5 Мб (как в 10.5).
#  Для гена: расстояние до ведущего SNP (0 — SNP внутри гена); r² ближайшего SNP панели интервала (IT или S1) с ведущим SNP.
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"; WORK <- "/work/scratch"; G <- file.path(WORK, "genes14")
IV <- fread(file.path(OUT, "tables/14_gene_intervals.csv"), encoding = "UTF-8")
p <- fread(file.path(G, "hits.paf"), header = FALSE, fill = TRUE, select = 1:12, col.names = c("tx", "qlen", "qs", "qe", "strand", "target", "tlen", "ts", "te", "nmatch", "alnlen", "mapq"))
p[, `:=`(identity = nmatch / alnlen, qcov = (qe - qs) / qlen)]
n_raw <- p[, .(n_hits_genes = uniqueN(sub("\\.[0-9]+$", "", tx))), by = .(interval_id = sub("\\|.*", "", target))]
p <- p[identity >= 0.98 & qcov >= 0.90]; p[, gene := sub("\\.[0-9]+$", "", tx)]
p <- p[order(-nmatch)][, .SD[1], by = .(gene, target)]
p[, `:=`(interval_id = sub("\\|.*", "", target), int_start = as.numeric(sub(".*:([0-9]+)-.*", "\\1", target)))]
p[, `:=`(langdon_start = int_start + ts - 1, langdon_end = int_start + te - 1)]
h <- fread(file.path(G, "cdna_headers.txt"), header = FALSE, sep = "\n", col.names = "line")
h[, tx := sub(" .*", "", line)][, svevo_pos := sub(".*chromosome:Svevo.v1:([^ ]+).*", "\\1", line)]
p <- merge(p, h[, .(tx, svevo_pos)], by = "tx")
p[, `:=`(svevo_chr = sub(":.*", "", svevo_pos), svevo_start = as.numeric(sub("^[^:]+:([0-9]+):.*", "\\1", svevo_pos)))]
p <- merge(p, IV[, .(interval_id, locus, chr, lead, lead_pos, panel, r2_thr)], by = "interval_id")
setorder(p, interval_id, langdon_start)
p[, offset := langdon_start - svevo_start]
rmed <- function(v, h = 15) sapply(seq_along(v), function(i) median(v[max(1, i - h):min(length(v), i + h)]))   # центрированное окно 31 ген (у краёв короче)
p[, roll_med := rmed(offset), by = interval_id]
p[, collinear := svevo_chr == chr & abs(offset - roll_med) <= 5e6]
CNT <- merge(n_raw, p[, .(n_pass_identity = .N, n_collinear = sum(collinear)), by = interval_id], by = "interval_id")
CNT[, n_HC := n_collinear][, note := "все гены Svevo v1 в Ensembl Plants — HC (logic_name=gff3_hc_genes); LC-набор в Ensembl отсутствует"]
g <- p[collinear == TRUE]
g[, dist_to_lead := fifelse(lead_pos >= langdon_start & lead_pos <= langdon_end, 0, pmin(abs(langdon_start - lead_pos), abs(langdon_end - lead_pos)))]
# r² ближайшего SNP с ведущим
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8")
P <- list(); for (s in c("stream1", "stream2")) { raw <- fread(file.path(WORK, "qc", s, "final.raw")); X <- as.matrix(raw[, -(1:6)]); rownames(X) <- raw$IID
  colnames(X) <- sub(":[^:]*:[^:]*_[^_]*$", "", colnames(X)); P[[s]] <- X; rm(raw) }
PN <- list(IT = P$stream2[intersect(OR[collection == "italy", TP_MO], rownames(P$stream2)), ], S1 = P$stream1[intersect(readLines(file.path(OUT, "tables/02_samples_stream1.txt")), rownames(P$stream1)), ])
g[, c("nearest_snp", "nearest_snp_dist", "r2_nearest_with_lead") := {
  X <- PN[[panel[1]]]; key <- colnames(X); k <- which(sub(":.*", "", key) == chr[1]); ps <- as.numeric(sub(".*:", "", key[k]))
  mid <- (langdon_start + langdon_end) / 2; j <- sapply(mid, function(m) k[which.min(abs(ps - m))])
  r2 <- suppressWarnings(cor(X[, lead[1]], X[, key[j], drop = FALSE], use = "pairwise.complete.obs"))[1, ]^2
  list(key[j], round(abs(as.numeric(sub(".*:", "", key[j])) - mid)), round(r2, 3)) }, by = interval_id]
fwrite(g[, .(interval_id, locus, chr, langdon_start, langdon_end, strand, svevo_gene = gene, svevo_transcript = tx, svevo_pos, identity = round(identity, 4), coverage = round(qcov, 3),
             lead, dist_to_lead, nearest_snp, nearest_snp_dist, r2_nearest_with_lead)], file.path(G, "genes.csv"))
fwrite(CNT, file.path(OUT, "tables/14_gene_counts_by_interval.csv")); print(CNT)
