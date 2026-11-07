# Этап 9.6 (R-часть): лучший гомеолог генов Svevo 4A на 4B / 5B / 7B Langdon; профиль вдоль Svevo 4A; где лежит Q4A.5.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
OUT <- "/work/project/gwas_frost"
T <- "/work/scratch/transloc"
p <- fread(file.path(T, "genes_hits.paf"), header = FALSE, fill = TRUE, col.names = c("tx","qlen","qs","qe","strand","t","tlen","ts","te","nm","alen","mapq"))
p <- p[(qe - qs) / qlen >= 0.5]
p[, ident := nm / alen]
best <- p[order(-nm)][, .SD[1], by = tx]
h <- data.table(line = readLines(file.path(T, "cdna4A_headers.txt")))
h[, tx := sub(" .*", "", line)][, svevo_pos := as.numeric(sub(".*chromosome:Svevo.v1:4A:([0-9]+):.*", "\\1", line))]
best <- merge(best, h[, .(tx, svevo_pos)], by = "tx")
best[, gene := sub("\\.[0-9]+$", "", tx)]; best <- best[order(-nm)][, .SD[1], by = gene]
best[, bin := floor(svevo_pos / 10e6) * 10]
prof <- best[, .N, by = .(bin, t)][, share := N / sum(N), by = bin]
g <- fread(file.path(OUT, "tables/05_genes_in_qtl.csv"))[QTL == "Q4A.5"]
g[, svevo_start := as.numeric(sub("^4A:([0-9]+):.*", "\\1", svevo_pos))]
q <- best[gene %in% g$svevo_gene]
QS <- q[, .N, by = t][, share := round(N / sum(N), 3)]
cat("гены Q4A.5 (Svevo 4A", round(min(g$svevo_start, na.rm = TRUE) / 1e6, 1), "–", round(max(g$svevo_start, na.rm = TRUE) / 1e6, 1), "Мб):\n"); print(QS)
fwrite(QS[, region := "Q4A.5 (Langdon 4A 618–650 Мб)"], file.path(OUT, "tables/09_4A_gene_synteny_Q4A5.csv"))
fwrite(prof[order(bin, -N)], file.path(OUT, "tables/09_4A_gene_synteny_profile.csv"))
lo <- quantile(g$svevo_start, .05, na.rm = TRUE) / 1e6; hi <- quantile(g$svevo_start, .95, na.rm = TRUE) / 1e6   # без единичных выбросов
pl <- ggplot(prof, aes(bin + 5, share, fill = t)) + geom_col(width = 9) + annotate("rect", xmin = lo, xmax = hi, ymin = -0.05, ymax = 1.05, alpha = 0, colour = "black") +
  scale_fill_manual(values = c("4B" = "#2a78d6", "5B" = "#eb6834", "7B" = "#1baf7a"), name = "лучший гомеолог") +
  labs(x = "позиция гена на Svevo 4A, Мб (окна 10 Мб)", y = "доля генов", title = "Гены Svevo 4A: где лучший гомеолог на Langdon (4B — исходный 4AL, 5B — сегмент 5AL, 7B — сегмент 7BS); рамка — гены Q4A.5") +
  theme_bw(base_size = 9)
ggsave(file.path(OUT, "figures/09_4A_gene_synteny.png"), pl, width = 11, height = 4, dpi = 170)
