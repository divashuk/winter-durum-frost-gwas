# Раздел 17, часть A, шаг 1. Интервалы поиска генов для локусов 1A 565.9, 1A 597.2–597.6, 4B 28.7, 5B 502.8–504.3 — правилом 14.4
#  (scripts/14_gene_intervals.R, функция block без изменений): ведущий SNP и панель r² — та из коллекций IT / RU (S1), где у локуса
#  наибольший −log10 p (14_loci_long.csv); блок — пока среди следующих 50 SNP есть SNP с r² > 0.2 к ведущему; границы — крайние такие SNP.
#  Остальные локусы раздела 17 (7B, 1B, 2B 767, 2B 776, 4A) — интервалы 14.4 (tables/14_gene_intervals.csv, I01, I02, I05–I07) без пересчёта.
# Выход: tables/17_gene_intervals.csv (все 9 локусов; новые — I09–I12).
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"; WORK <- "/work/scratch"
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8")
P <- list(); for (s in c("stream1", "stream2")) { raw <- fread(file.path(WORK, "qc", s, "final.raw")); X <- as.matrix(raw[, -(1:6)]); rownames(X) <- raw$IID
  colnames(X) <- sub(":[^:]*:[^:]*_[^_]*$", "", colnames(X)); P[[s]] <- X; rm(raw) }
P$S1 <- P$stream1[intersect(readLines(file.path(OUT, "tables/02_samples_stream1.txt")), rownames(P$stream1)), ]
P$IT <- P$stream2[intersect(OR[collection == "italy", TP_MO], rownames(P$stream2)), ]
LONG <- fread(file.path(OUT, "tables/14_loci_long.csv"), encoding = "UTF-8"); STAT <- fread(file.path(OUT, "tables/14_loci_status.csv"), encoding = "UTF-8")
block <- function(pn, lead, thr, gap) { X <- P[[pn]]; key <- colnames(X); ch <- sub(":.*", "", lead); k <- which(sub(":.*", "", key) == ch)
  ps <- as.numeric(sub(".*:", "", key[k])); o <- order(ps); k <- k[o]; ps <- ps[o]; i0 <- which(key[k] == lead)
  r2 <- suppressWarnings(cor(X[, lead], X[, k], use = "pairwise.complete.obs"))[1, ]^2; r2[is.na(r2)] <- 0; hit <- r2 > thr
  lo <- i0; j <- i0; while (j > 1) { w <- max(1, j - gap):(j - 1); if (any(hit[w])) { j <- max(w[hit[w]]); lo <- j } else break }
  hi <- i0; j <- i0; while (j < length(k)) { w <- (j + 1):min(length(k), j + gap); if (any(hit[w])) { j <- min(w[hit[w]]); hi <- j } else break }
  data.table(start = ps[lo], end = ps[hi], n_snp_block = hi - lo + 1, n_snp_r2 = sum(hit[lo:hi])) }
loci_new <- c("1A 565.9–565.9", "1A 597.2–597.6", "4B 28.7–28.7", "5B 502.8–504.3")
NEW <- rbindlist(lapply(seq_along(loci_new), function(i) { l <- loci_new[i]
  x <- LONG[locus == l & set %in% c("IT", "RU (S1)") & !is.na(mlog10p)][which.max(mlog10p)]
  pn <- ifelse(x$set == "IT", "IT", "S1"); lead <- sub(":[^:]*:[^:]*$", "", x$lead)
  cbind(data.table(interval_id = sprintf("I%02d", 8 + i), locus = l, chr = x$chr, lead_set = x$set, lead = lead, lead_pos = as.numeric(sub(".*:", "", lead)),
                   lead_mlp = x$mlog10p, status = STAT[locus == l, status], panel = pn, r2_thr = 0.2), block(pn, lead, 0.2, 50)) }))
NEW[, `:=`(start_Mb = round(start / 1e6, 3), end_Mb = round(end / 1e6, 3), length_Mb = round((end - start) / 1e6, 2))]
OLD <- fread(file.path(OUT, "tables/14_gene_intervals.csv"), encoding = "UTF-8")[interval_id %in% c("I01", "I02", "I05", "I06", "I07", "I08")]
IV <- rbind(OLD, NEW, fill = TRUE)
IV[, source := fifelse(interval_id %in% NEW$interval_id, "раздел 17 (правило 14.4)", "раздел 14.4 (tables/14_gene_intervals.csv)")]
fwrite(IV, file.path(OUT, "tables/17_gene_intervals.csv")); print(IV[, .(interval_id, locus, lead_set, lead, lead_mlp, panel, start_Mb, end_Mb, length_Mb, n_snp_r2)])
