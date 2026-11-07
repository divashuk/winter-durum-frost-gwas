# Этап 5.1б. Сводная таблица позиций известных локусов на Langdon: tables/05_known_loci_positions.csv
#  гены — лучший хит minimap2 на ожидаемой хромосоме (идентичность, покрытие запроса);
#  KASP — ампликон по праймерам (аллель-специфичный + обратный на разных цепях, <= 1 кб); если праймеры не дали
#         однозначной пары — позиция SNP с максимальным r2 к маркеру (tables/01_marker_snp_collinearity.csv, поток 1).
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"
q  <- fread(file.path(OUT, "tables/05_known_genes_queries.tsv"), header = FALSE, col.names = c("locus", "acc", "exp_chr", "note"), sep = "\t")
h  <- fread(file.path(OUT, "tables/05_known_genes_hits_raw.tsv"))
h[, locus := sub("\\|.*", "", query)]
h[startsWith(locus, "Rht-A1"), locus := q[startsWith(locus, "Rht-A1"), locus]]   # имя с пробелом обрезано minimap2
h <- merge(h, q, by = "locus")
h <- h[!(locus == "Fr-Am2_BAC" & qcov < 0.3)]   # короткие частичные хиты BAC (повторы) не считаем позицией локуса
h[, on_exp := mapply(function(c, e) c %in% strsplit(e, ",")[[1]], chr, exp_chr)]
G <- h[on_exp == TRUE][order(-nmatch)][, .SD[1], by = .(locus, chr)]
G <- G[, .(locus, chr, start = tstart, end = tend, method = "minimap2 asm20", identity, qcov, source = paste(acc, note))]
# Rht-A1 запрос: хит на 4B = гомеолог Rht-B1
G[locus == "Rht-A1(→Rht-B1 гомеолог)" & chr == "4B", locus := "Rht-B1 (по гомеологу Rht-A1)"]
G[locus == "Rht-A1(→Rht-B1 гомеолог)" & chr == "4A", locus := "Rht-A1"]
G[locus %in% c("Fr-Am2_BAC", "CBF12", "CBF15") & chr == "5A", locus := paste0(locus, " → Fr-A2")]
G[locus %in% c("Fr-Am2_BAC", "CBF12", "CBF15") & chr == "5B", locus := paste0(locus, " → Fr-B2")]
# KASP по праймерам
kb <- fread(file.path(OUT, "tables/05_kasp_primer_hits_raw.tsv"), header = FALSE,
            col.names = c("q", "chr", "pident", "len", "mm", "gap", "qs", "qe", "ss", "se", "ev", "bit", "qlen"))
kb[, `:=`(marker = sub("_(ASfam|AShex|R|A|C)$", "", q), role = sub(".*_", "", q), strand = fifelse(se > ss, "+", "-"), pos = pmin(ss, se))]
KS <- rbindlist(lapply(unique(kb$marker), function(m) {
  a <- kb[marker == m & role %in% c("ASfam", "AShex", "A")]; r <- kb[marker == m & role %in% c("R", "C")]
  pr <- merge(a, r, by = "chr", allow.cartesian = TRUE)[strand.x != strand.y & abs(pos.x - pos.y) <= 1000]
  if (!nrow(pr)) return(data.table(locus = m, chr = NA, start = NA, end = NA, method = "праймеры: пары нет", identity = NA, qcov = NA, source = "KASP"))
  pr[, .(locus = m, chr, start = pmin(pos.x, pos.y), end = pmax(pos.x, pos.y) + 30, method = sprintf("праймеры: %d пар(ы)", .N),
         identity = NA, qcov = NA, source = "KASP")][1]
}))
# KASP по эмпирическому LD
cl <- fread(file.path(OUT, "tables/01_marker_snp_collinearity.csv"))[stream == "stream1"]
cl[, `:=`(chr = sub(":.*", "", top_snp), pos = as.numeric(sapply(strsplit(top_snp, ":"), `[`, 2)))]
KL <- cl[, .(locus = paste0(marker, " (LD)"), chr, start = pos, end = pos, method = sprintf("max r2 с SNP = %.2f", max_r2),
             identity = NA, qcov = NA, source = "поток 1, r2 доза маркера ~ SNP")]
P <- rbind(G, KS, KL, fill = TRUE)
fwrite(P, file.path(OUT, "tables/05_known_loci_positions.csv"))
print(P)
