# Этап 10.1. Локусы для поиска генов-кандидатов и границы поиска.
#  широкий интервал  = интервал QTL (этап 9 после поправки на происхождение; если в 09 его нет — этап 6; для маркеров Mor — LD r2 >= 0.2 как в этапе 5);
#  узкий интервал    = LD-блок ведущего SNP (r2 >= 0.5 с ведущим SNP внутри широкого интервала) ± 1 Мб, обрезанный широким интервалом,
#                      но не уже ведущий SNP ± 1 Мб.
#  LD считается в той выборке, где найден ведущий сигнал: S1 (поток 1), S2 (поток 2) или только итальянские образцы потока 2.
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
L <- fread(text = "
locus_id|name|chr|category|lead_pos|ld_set|qtl_start|qtl_end|qtl_source
L01|Q5A.1 (пик 5A, Fr-A2)|5A|известный локус; подтверждён|514681016|S1|475360000|545501466|09_qtl_intervals 5A.3
L02|Fr-A2, кластер CBF (по последовательности)|5A|известный локус (положение гена)|534103897|S1|533103897|535181703|05_known_loci_positions: Fr-Am2 BAC 534.10–534.18 ± 1 Мб
L03|Q5B.1|5B|вероятный|521642333|S2|520973288|522742783|06_qtl_summary Q5B.1
L04|Q5B.2|5B|вероятный|535618213|S2|525925602|539497626|06_qtl_summary Q5B.2
L05|Q4A.5|4A|вероятный (внутри RU)|642497807|S1|617944415|650214537|06_qtl_summary Q4A.5
L06|Q7B.1|7B|вероятный (внутри IT)|348683836|IT|348183836|349183836|09_qtl_intervals 7B.3
L07|Q2B.4|2B|вероятный (S2, ниже Бонферрони)|775875145|S2|763241794|782150000|09_qtl_intervals 2B.4
L08|Q1B.6|1B|потенциальный|653803203|S1|650196761|660809546|06_qtl_summary Q1B.6
L09|Q5B.3 (VRN-B1)|5B|потенциальный|598637874|S1|580408612|608478255|06_qtl_summary Q5B.3 (в 09 значимых SNP нет)
L10|Q5A.4 (у VRN-A1)|5A|потенциальный|615862573|S2|609173502|622369126|06_qtl_summary Q5A.4
L11|Q2A.1 (Ppd-A1)|2A|потенциальный|41862224|S2|30767399|58925206|06_qtl_summary Q2A.1
L12|Mor534 (5B:550.23)|5B|потенциальный (маркер)|550225979|S1|NA|NA|LD r2>=0.2 вокруг маркера
L13|Mor587 (1B:54.74)|1B|потенциальный (маркер)|54736990|S1|NA|NA|LD r2>=0.2 вокруг маркера
L14|Mor866 (7B:610.78)|7B|потенциальный (маркер)|610783371|S1|NA|NA|LD r2>=0.2 вокруг маркера
L15|1A 563.8–569.6 (Q1A.1)|1A|новый, > Meff MLM/GEMMA (IT)|565864449|IT|563798332|569596347|09_qtl_intervals 1A.2
L16|6A 550.5–557.6|6A|новый, > Meff MLM/GEMMA (S2 bin)|557124931|S2|550470000|557620000|09_qtl_intervals 6A.4
L17|6B 570.5–574.1|6B|новый, > Meff MLM/GEMMA (S2 bin)|573622970|S2|570530000|574120000|09_qtl_intervals 6B.4
L18|6B 710.5–720.1|6B|новый, > Meff MLM/GEMMA (S2 bin)|715018217|S2|710460000|720130000|09_qtl_intervals 6B.7
L19|5B 516.9–521.0|5B|новый, > Meff MLM/GEMMA (S2 Bn)|519634162|S2|516890000|520960000|09_qtl_intervals 5B.3
L20|7B 139.7–154.7|7B|новый, > Meff MLM/GEMMA (S2 Bn)|148831119|S2|139730000|154660000|09_qtl_intervals 7B.2
L21|4B 28.2–30.2|4B|новый, > Meff MLM/GEMMA (S1)|28712352|S1|28210000|30230000|09_qtl_intervals 4B.1
L22|Fr-B2, CBF12/CBF15 (по последовательности)|5B|известный локус (положение гена)|524074424|S2|523074424|525097128|05_known_loci_positions: CBF15/CBF12 → Fr-B2 524.07–524.10 ± 1 Мб (между Q5B.1 и Q5B.2)
L23|VRN-A1 (по последовательности)|5A|известный локус (положение гена)|598268462|S2|597268462|599273527|05_known_loci_positions: VRN-A1 598.27 ± 1 Мб (вне Q5A.4)
", sep = "|")
# точные границы интервалов 09
Q9 <- fread(file.path(OUT, "tables/09_qtl_intervals.csv"))
fix9 <- function(id, q) { r <- Q9[Q == q]; L[locus_id == id, `:=`(qtl_start = r$start, qtl_end = r$end)] }
fix9("L01", "5A.3"); fix9("L06", "7B.3"); fix9("L07", "2B.4"); fix9("L15", "1A.2"); fix9("L16", "6A.4"); fix9("L17", "6B.4")
fix9("L18", "6B.7"); fix9("L19", "5B.3"); fix9("L20", "7B.2"); fix9("L21", "4B.1")
S2T <- fread(file.path(OUT, "tables/09_s2_traits.csv")); it_ids <- unique(S2T[coll == "IT", TP_MO])
geno <- list()
for (s in c("stream1", "stream2")) { raw <- fread(file.path(WORK, "qc", s, "final.raw"))
  geno[[s]] <- list(X = as.matrix(raw[, -(1:6)]), id = raw$IID, key = sub(":[^:]*:[^:]*$", "", sub("_[^_]*$", "", names(raw)[-(1:6)]))); rm(raw) }
LDD <- rbindlist(lapply(c("stream1", "stream2"), function(s) fread(file.path(OUT, sprintf("tables/03_ld_decay_%s.csv", s)))))
r2_with <- function(set, ch, ps, lo, hi) {
  g <- geno[[if (set == "S1") "stream1" else "stream2"]]
  rows <- if (set == "IT") which(g$id %in% it_ids) else seq_along(g$id)
  kk <- which(startsWith(g$key, paste0(ch, ":"))); pp <- as.numeric(sub(".*:", "", g$key[kk]))
  j0 <- kk[which.min(abs(pp - ps))]; w <- kk[pp >= lo & pp <= hi]
  r2 <- suppressWarnings(cor(g$X[rows, j0], g$X[rows, w, drop = FALSE], use = "pairwise.complete.obs")^2)[1, ]
  list(lead_snp = g$key[j0], pos = as.numeric(sub(".*:", "", g$key[w])), r2 = r2, n = length(rows)) }
res <- list()
for (k in seq_len(nrow(L))) { x <- L[k]
  if (is.na(x$qtl_start)) {   # маркер Mor: интервал по LD r2 >= 0.2 в пределах критической дистанции LD (этап 3)
    D <- LDD[stream == "stream1" & subgenome == substr(x$chr, 2, 2), dist_critical_Mb] * 1e6
    z <- r2_with(x$ld_set, x$chr, x$lead_pos, x$lead_pos - D, x$lead_pos + D); inn <- z$pos[!is.na(z$r2) & z$r2 >= 0.2]
    x$qtl_start <- max(1, min(c(inn, x$lead_pos - 5e5))); x$qtl_end <- max(c(inn, x$lead_pos + 5e5)) }
  z <- r2_with(x$ld_set, x$chr, x$lead_pos, x$qtl_start, x$qtl_end)
  b <- z$pos[!is.na(z$r2) & z$r2 >= 0.5]
  ns <- max(x$qtl_start, min(c(b, x$lead_pos)) - 1e6); ne <- min(x$qtl_end, max(c(b, x$lead_pos)) + 1e6)
  ns <- min(ns, x$lead_pos - 1e6); ne <- max(ne, x$lead_pos + 1e6)
  res[[k]] <- cbind(x, data.table(lead_snp_nearest = z$lead_snp, n_ld = z$n, n_snp_qtl = length(z$pos), n_snp_r2_05 = length(b),
                                  ld_block_start = if (length(b)) min(b) else NA, ld_block_end = if (length(b)) max(b) else NA,
                                  narrow_start = ns, narrow_end = ne)) }
R <- rbindlist(res)
R[, `:=`(wide_Mb = round((qtl_end - qtl_start) / 1e6, 2), narrow_Mb = round((narrow_end - narrow_start) / 1e6, 2))]
fwrite(R, file.path(OUT, "tables/10_loci.csv"))
print(R[, .(locus_id, name, chr, lead_snp_nearest, qtl_start, qtl_end, wide_Mb, n_snp_r2_05, ld_block_start, ld_block_end, narrow_start, narrow_end, narrow_Mb)])
