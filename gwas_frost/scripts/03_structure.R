# Этап 3. Структура популяции: PCA-графики, scree, родство VanRaden, спад LD по субгеномам A/B, Meff (simpleM).
# Использование: Rscript 03_structure.R stream1|stream2
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
s    <- commandArgs(TRUE)[1]
OUT  <- "/work/project/gwas_frost"
WORK <- "/work/scratch"
set.seed(20261025)
idm <- fread(file.path(OUT, "tables/01_id_map.csv"), encoding = "UTF-8")

## ---------- PCA ----------
ev  <- scan(file.path(WORK, "pca", s, "pca.eigenval"), quiet = TRUE)
pcs <- fread(file.path(WORK, "pca", s, "pca.eigenvec"))
if ("IID" %in% names(pcs)) pcs[, 1 := NULL] else setnames(pcs, 1, "IID")
pcs <- merge(pcs, idm[, .(IID = TP_MO, name_pheno, collection, origin3, origin_group, growth_habit)], by = "IID")
pcs[, habit := fifelse(collection == "RU_winter", "озимый", "неизвестен (итал. колл.)")]
pve <- 100 * ev / sum(ev)   # доля среди 20 вычисленных PC (plink2 --pca 20)
fwrite(data.table(PC = seq_along(ev), eigenval = ev, pct_of_20 = pve), file.path(OUT, sprintf("tables/03_pca_eigen_%s.csv", s)))
fwrite(pcs, file.path(OUT, sprintf("tables/03_pca_scores_%s.csv", s)))
pal <- c("Россия" = "#1f77b4", "прочее" = "#ff7f0e", "Италия (итал. коллекция)" = "#2ca02c", "Италия" = "#2ca02c", "нет данных" = "grey50")
mk <- function(a, b) {
  ggplot(pcs, aes(.data[[a]], .data[[b]], colour = origin3, shape = habit)) + geom_point(size = 2, alpha = .8) +
    scale_colour_manual(values = pal, name = "Происхождение") + scale_shape_manual(values = c(16, 17), name = "Тип развития") +
    labs(x = sprintf("%s (%.1f%%)", a, pve[as.integer(sub("PC", "", a))]), y = sprintf("%s (%.1f%%)", b, pve[as.integer(sub("PC", "", b))])) +
    theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.box = "vertical")
}
g <- gridExtra::arrangeGrob(mk("PC1", "PC2") + ggtitle(sprintf("%s: PCA (%d образцов)", s, nrow(pcs))), mk("PC1", "PC3"), ncol = 2)
ggsave(file.path(OUT, sprintf("figures/03_pca_%s.png", s)), g, width = 13, height = 5, dpi = 200)
if (s == "stream1") {  # дополнительно: группы происхождения из Табл. S4
  p2 <- ggplot(pcs, aes(PC1, PC2, colour = origin_group)) + geom_point(size = 2) + theme_bw(base_size = 10) +
    labs(title = "stream1: PCA, группы происхождения (Табл. S4 статьи)", colour = "Группа")
  ggsave(file.path(OUT, "figures/03_pca_stream1_origin_group.png"), p2, width = 7, height = 5, dpi = 200)
}
sc <- ggplot(data.table(PC = seq_along(ev), pve = pve), aes(PC, pve)) + geom_line() + geom_point() + theme_bw() +
  labs(y = "% дисперсии (из 20 PC)", title = sprintf("%s: scree plot", s))
ggsave(file.path(OUT, sprintf("figures/03_scree_%s.png", s)), sc, width = 5, height = 3.5, dpi = 200)

## ---------- Родство VanRaden ----------
raw <- fread(file.path(WORK, "qc", s, "final.raw"))
ids <- raw$IID; X <- as.matrix(raw[, -(1:6)]); rm(raw)
p   <- colMeans(X) / 2
Z   <- sweep(X, 2, 2 * p)
K   <- tcrossprod(Z) / (2 * sum(p * (1 - p)))
dimnames(K) <- list(ids, ids)
fwrite(as.data.table(K, keep.rownames = "IID"), file.path(WORK, "kin", sprintf("K_vanraden_%s.csv", s)))
hc <- hclust(as.dist(max(K) - K), "average"); o <- hc$order
png(file.path(OUT, sprintf("figures/03_kinship_%s.png", s)), width = 1400, height = 1300, res = 160)
col_side <- c(RU_winter = "#1f77b4", italy = "#2ca02c")[idm$collection[match(ids, idm$TP_MO)]]
heatmap(K[o, o], Rowv = NA, Colv = NA, scale = "none", labRow = NA, labCol = NA, RowSideColors = col_side[o],
        col = hcl.colors(50, "YlOrRd", rev = TRUE), main = sprintf("%s: родство VanRaden (синий — РФ, зелёный — итал.)", s))
dev.off()
kd <- K[upper.tri(K)]
cat(sprintf("%s: K VanRaden — диагональ медиана %.2f; вне диагонали: медиана %.3f, 99%% %.3f\n", s, median(diag(K)), median(kd), quantile(kd, .99)))

## ---------- Спад LD ----------
ld <- fread(cmd = sprintf("zcat %s | awk 'NR==1 || NR%%4==1'", file.path(WORK, "ld", s, "r2.ld.gz")),
            select = c("CHR_A", "BP_A", "BP_B", "R2"))   # каждая 4-я пара — для скорости
ld[, dist := abs(BP_B - BP_A)]
chr <- c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B","7A","7B")
ld[, sub := ifelse(grepl("A$", CHR_A), "A", "B")]   # хромосомы в bim/ld — строки «1A»…«7B»
# критический r2: 95-й перцентиль r2 между несцепленными SNP (разные хромосомы), Breseghello & Sorrells 2006
bim <- fread(file.path(WORK, "qc", s, "final.bim"), header = FALSE)
i1 <- sample(ncol(X), 60000, TRUE); i2 <- sample(ncol(X), 60000, TRUE); keep <- bim$V1[i1] != bim$V1[i2]
r2u <- mapply(function(a, b) suppressWarnings(cor(X[, a], X[, b]))^2, i1[keep], i2[keep])
crit <- quantile(r2u, .95, na.rm = TRUE)
n <- nrow(X)
hw <- function(d, rho) { C <- rho * d; ((10 + C) / ((2 + C) * (11 + C))) * (1 + ((3 + C) * (12 + 12 * C + C^2)) / (n * (2 + C) * (11 + C))) }
res <- list(); curves <- list()
for (g in c("A", "B")) {
  dd <- ld[sub == g & dist <= 1e8]
  ds <- dd[sample(.N, min(.N, 300000))]
  sse <- function(lr) sum((ds$R2 - hw(ds$dist, 10^lr))^2)   # МНК по log10(rho)
  rho <- 10^optimize(sse, c(-12, 0))$minimum
  grid <- 10^seq(2, 8, length.out = 2000); pr <- hw(grid, rho)
  rmax <- max(pr); half <- grid[which(pr <= rmax / 2)[1]]; dc <- grid[which(pr <= crit)[1]]
  d02  <- grid[which(pr <= 0.2)[1]]
  bins <- dd[, .(r2_mean = mean(R2), r2_med = median(R2), n = .N), by = .(bin = cut(dist, c(0, 1e4, 1e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8), include.lowest = TRUE))]
  res[[g]] <- data.table(stream = s, subgenome = g, n_pairs_used = nrow(dd), rho = rho, r2_max_fit = rmax, critical_r2_unlinked95 = crit,
                         dist_half_max_Mb = half / 1e6, dist_critical_Mb = dc / 1e6, dist_r2_0.2_Mb = d02 / 1e6)
  curves[[g]] <- data.table(sub = g, dist = grid, r2 = pr)
  fwrite(bins[order(bin)], file.path(OUT, sprintf("tables/03_ld_bins_%s_%s.csv", s, g)))
}
R <- rbindlist(res); print(R)
fwrite(R, file.path(OUT, sprintf("tables/03_ld_decay_%s.csv", s)))
pl <- ggplot(ld[dist <= 1e8][sample(.N, min(.N, 150000))], aes(dist / 1e6, R2)) + geom_point(alpha = .03, size = .3) +
  geom_line(data = rbindlist(curves), aes(dist / 1e6, r2, colour = sub), linewidth = 1) +
  geom_hline(yintercept = crit, linetype = 2) + facet_wrap(~sub) + scale_x_log10() +
  labs(x = "Расстояние, Мб (log)", y = expression(r^2), colour = "Субгеном",
       title = sprintf("%s: спад LD (Hill–Weir); пунктир — критический r2 = %.3f", s, crit)) + theme_bw(base_size = 10)
ggsave(file.path(OUT, sprintf("figures/03_ld_decay_%s.png", s)), pl, width = 10, height = 4.5, dpi = 200)

## ---------- Meff (simpleM, Gao et al. 2008): по хромосомам, 99.5% дисперсии ----------
meff <- rbindlist(lapply(unique(bim$V1), function(cc) {
  idx <- which(bim$V1 == cc); Xc <- X[, idx, drop = FALSE]
  Xc <- Xc[, apply(Xc, 2, sd) > 0, drop = FALSE]
  ev <- eigen(cor(Xc), symmetric = TRUE, only.values = TRUE)$values
  data.table(chr = cc, n_snp = ncol(Xc), Meff = which(cumsum(ev) / sum(ev) >= 0.995)[1])
}))
meff <- rbind(meff, data.table(chr = "ИТОГО", n_snp = sum(meff$n_snp), Meff = sum(meff$Meff)))
fwrite(meff, file.path(OUT, sprintf("tables/03_meff_%s.csv", s)))
cat(sprintf("%s: M = %d, Meff = %d; Бонферрони 0.05/M: -log10 p = %.2f; 0.05/Meff: -log10 p = %.2f\n",
            s, ncol(X), tail(meff$Meff, 1), -log10(0.05 / ncol(X)), -log10(0.05 / tail(meff$Meff, 1))))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_03_structure.txt"))
