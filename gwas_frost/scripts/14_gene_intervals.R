# Этап 14.4 (часть 1). Интервалы поиска генов: блок LD вокруг ведущего SNP.
#  Локусы: общие локусы 14.2, у которых есть пик (>= 4) в S2 или в мета IT + RU + F2 и статус «подтверждён / вероятный / предположительный», плюс 1B (Q1B.6)
#  и 4A (Q4A.5, для проверки 14.5) независимо от статуса.
#  Ведущий SNP и панель r²: та из коллекций IT / RU (S1), где у локуса наибольший −log10 p (14_loci_long.csv); её ведущий SNP; r² — внутри этой коллекции.
#   Причина: в объединённой панели S2 частные гаплотипы двух коллекций дают r² > 0.2 на сотни Мб (структура) — первый вариант с панелью S2 давал 5A 52–629 Мб и 7B 108–587 Мб.
#   F2 как панель не используется (LD в F2 тянется на всю хромосому). 1B и 4A — RU (S1).
#  Блок: от ведущего SNP в обе стороны по SNP хромосомы (по порядку позиций), пока среди следующих 50 SNP есть хотя бы один с r² > 0.2;
#  границы — крайние SNP с r² > 0.2 (участки без SNP разрывом не считаются). Узкое окно 7B: то же для r² > 0.6 вокруг Q7B.1 (окно — 20 SNP).
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"; WORK <- "/work/scratch"
LOC <- fread(file.path(OUT, "tables/14_loci_by_set_raw.csv"), encoding = "UTF-8"); STAT <- fread(file.path(OUT, "tables/14_loci_status.csv"), encoding = "UTF-8")
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8")
P <- list(); for (s in c("stream1", "stream2")) { raw <- fread(file.path(WORK, "qc", s, "final.raw")); X <- as.matrix(raw[, -(1:6)]); rownames(X) <- raw$IID
  colnames(X) <- sub(":[^:]*:[^:]*_[^_]*$", "", colnames(X)); P[[s]] <- X; rm(raw) }
P$S2 <- P$stream2[intersect(readLines(file.path(OUT, "tables/02_samples_stream2.txt")), rownames(P$stream2)), ]
P$S1 <- P$stream1[intersect(readLines(file.path(OUT, "tables/02_samples_stream1.txt")), rownames(P$stream1)), ]
P$IT <- P$stream2[intersect(OR[collection == "italy", TP_MO], rownames(P$stream2)), ]
LONG <- fread(file.path(OUT, "tables/14_loci_long.csv"), encoding = "UTF-8")
ok_status <- STAT[grepl("^(подтверждён|вероятный|предположительный)", status), locus]
loci_sel <- unique(LOC[set %in% c("мета IT + RU + F2", "S2") & locus %in% ok_status, locus])
loci_sel <- unique(c(loci_sel, LOC[locus %like% "^1B 650", locus][1], LOC[locus %like% "^4A 63", locus][1]))
sel <- rbindlist(lapply(loci_sel, function(l) { x <- LONG[locus == l & set %in% c("IT", "RU (S1)") & !is.na(mlog10p)]
  if (grepl("^(1B 650|4A 63)", l)) x <- x[set == "RU (S1)"]
  x <- x[which.max(mlog10p)]; data.table(locus = l, chr = x$chr, set = x$set, lead = x$lead, lead_pos = round(x$lead_Mb * 1e6), lead_mlp = x$mlog10p,
                                         panel = ifelse(x$set == "IT", "IT", "S1")) }))
block <- function(pn, lead, thr, gap) { X <- P[[pn]]; key <- colnames(X); ch <- sub(":.*", "", lead); k <- which(sub(":.*", "", key) == ch)
  ps <- as.numeric(sub(".*:", "", key[k])); o <- order(ps); k <- k[o]; ps <- ps[o]; i0 <- which(key[k] == lead)
  r2 <- suppressWarnings(cor(X[, lead], X[, k], use = "pairwise.complete.obs"))[1, ]^2; r2[is.na(r2)] <- 0; hit <- r2 > thr
  lo <- i0; j <- i0; while (j > 1) { w <- max(1, j - gap):(j - 1); if (any(hit[w])) { j <- max(w[hit[w]]); lo <- j } else break }
  hi <- i0; j <- i0; while (j < length(k)) { w <- (j + 1):min(length(k), j + gap); if (any(hit[w])) { j <- min(w[hit[w]]); hi <- j } else break }
  data.table(start = ps[lo], end = ps[hi], n_snp_block = hi - lo + 1, n_snp_r2 = sum(hit[lo:hi])) }
IV <- rbindlist(lapply(seq_len(nrow(sel)), function(i) { x <- sel[i]; lead <- sub(":[^:]*:[^:]*$", "", x$lead)
  b <- block(x$panel, lead, 0.2, 50)
  cbind(data.table(interval_id = sprintf("I%02d", i), locus = x$locus, chr = x$chr, lead_set = x$set, lead = lead, lead_pos = x$lead_pos, lead_mlp = x$lead_mlp,
                   status = STAT[locus == x$locus, status], panel = x$panel, r2_thr = 0.2), b) }))
q7 <- block("IT", "7B:348683836", 0.6, 20)
IV <- rbind(IV, cbind(data.table(interval_id = sprintf("I%02d", nrow(IV) + 1), locus = IV[chr == "7B"]$locus[1], chr = "7B", lead_set = "IT (Q7B.1, узкое окно)", lead = "7B:348683836",
                                 lead_pos = 348683836, lead_mlp = NA, status = IV[chr == "7B"]$status[1], panel = "IT", r2_thr = 0.6), q7), fill = TRUE)
IV[, `:=`(start_Mb = round(start / 1e6, 3), end_Mb = round(end / 1e6, 3), length_Mb = round((end - start) / 1e6, 2))]
fwrite(IV, file.path(OUT, "tables/14_gene_intervals.csv")); print(IV[, .(interval_id, locus, lead_set, lead, lead_mlp, status, panel, r2_thr, start_Mb, end_Mb, length_Mb, n_snp_r2)])
