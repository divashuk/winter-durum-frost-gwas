# Раздел 17, шаг 7а. Сведения о SNP для KASP типа (б).
#  1) SNP-заместители 5B:535618213: SNP панели IT (поток 2, только итальянские образцы, генотипы после QC/импутации — как r² в 14.4)
#     в ± 2 Мб с r² ≥ 0.8 к 5B:535618213, по возрастанию расстояния -> s17/kasp/proxies_535618213.tsv;
#  2) список SNP типа (б) раздела 17 -> s17/kasp/snps_b.tsv: ведущие SNP локусов части A (7B — оба: 7B:348683836 и 7B:233439923),
#     5B:535618213 и его заместители.
suppressPackageStartupMessages(library(data.table))
OUT <- "/work/project/gwas_frost"; WORK <- "/work/scratch"
K <- file.path(WORK, "s17/kasp"); dir.create(K, showWarnings = FALSE, recursive = TRUE)
OR <- fread(file.path(OUT, "tables/09_origin_by_sample.csv"), encoding = "UTF-8")
raw <- fread(file.path(WORK, "qc/stream2/final.raw")); X <- as.matrix(raw[, -(1:6)]); rownames(X) <- raw$IID
colnames(X) <- sub(":[^:]*:[^:]*_[^_]*$", "", colnames(X)); X <- X[intersect(OR[collection == "italy", TP_MO], rownames(X)), ]
key <- colnames(X); k <- which(sub(":.*", "", key) == "5B"); ps <- as.numeric(sub(".*:", "", key[k]))
k <- k[abs(ps - 535618213) <= 2e6]
r2 <- suppressWarnings(cor(X[, "5B:535618213"], X[, k], use = "pairwise.complete.obs"))[1, ]^2
P <- data.table(snp = key[k], pos = as.numeric(sub(".*:", "", key[k])), r2 = round(r2, 3))[snp != "5B:535618213" & r2 >= 0.8]
P[, dist := abs(pos - 535618213)]; setorder(P, dist); fwrite(P, file.path(K, "proxies_535618213.tsv"), sep = "\t"); print(P)
leads <- data.table(snp = c("7B:348683836", "7B:233439923", "1B:653803203", "2B:767264315", "2B:775875145", "4A:645178456",
                            "1A:565864449", "1A:597601155", "4B:28712352", "5B:504253234"),
                    locus = c("7B 145.8–348.7", "7B 145.8–348.7", "1B 650.2–660.8", "2B 767.1–769.5", "2B 775.9–776.2", "4A 631.5–645.2",
                              "1A 565.9", "1A 597.2–597.6", "4B 28.7", "5B 502.8–504.3"), role = "ведущий SNP локуса")
S <- rbind(leads, data.table(snp = "5B:535618213", locus = "5B Fr-B2 (IT, rint / бинарный)", role = "повтор 16.7 (Tm 57–63, ампликон ≤ 150)"),
           P[, .(snp, locus = "5B Fr-B2 (IT, rint / бинарный)", role = paste0("заместитель 5B:535618213 (r² = ", r2, ", ", round(dist / 1e3, 1), " кб)"))])
S[, `:=`(chr = sub(":.*", "", snp), pos = as.numeric(sub(".*:", "", snp)))]
fwrite(S, file.path(K, "snps_b.tsv"), sep = "\t"); print(S)
