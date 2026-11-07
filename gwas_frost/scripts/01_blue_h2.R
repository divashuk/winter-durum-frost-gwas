# Этап 1.2. BLUE по генотипу и наследуемость (повторяемость) морозостойкости.
# Модель BLUE: y ~ 0 + G + (1|env), env = год × режим (одна коробка на образец на режим — повторностей внутри среды нет).
# Наследуемость: y ~ (1|G) + (1|env); H2_entry = Vg / (Vg + Ve / n_h), n_h — гармоническое среднее числа сред на образец;
# H2_Cullis = 1 - mean(PEV_diff) / (2 Vg).
# Шкалы: raw (%, исходная) и ang (asin(sqrt(p)), как в статье). Каждый признак = отдельный GWAS.
suppressPackageStartupMessages({library(lme4); library(data.table); library(ggplot2)})
OUT <- "/work/project/gwas_frost"
pl  <- fread(file.path(OUT, "tables/01_pheno_long.csv"), encoding = "UTF-8")
idm <- fread(file.path(OUT, "tables/01_id_map.csv"), encoding = "UTF-8")
pl  <- pl[use == "да"]
pl[, ang := asin(sqrt(pmin(pmax(value, 0), 100) / 100))]
pl[, raw := value]

TRAITS <- list(  # имя -> (коллекция, среды)
  FT_all = list("RU_winter", unique(pl[collection == "RU_winter", env])),
  FT_m12 = list("RU_winter", c("RU_2021_-12", "RU_2022_-12", "RU_2023_-12")),
  FT_m15 = list("RU_winter", c("RU_2020_-15", "RU_2021_-15", "RU_2022_-15", "RU_2023_-15")),
  FT_m16 = list("RU_winter", c("RU_2020_-16", "RU_2021_-16", "RU_2022_-16")),
  FT_m17 = list("RU_winter", c("RU_2021_-17")),
  FT_m14 = list("RU_winter", c("RU_2022_-14")),
  ITFT_m12 = list("italy", c("IT_2020_-12", "IT_2021_-12"))
)
skew <- function(x) { x <- x[is.finite(x)]; mean((x - mean(x))^3) / sd(x)^3 }

blues <- list(); h2 <- list()
for (tr in names(TRAITS)) {
  d <- pl[collection == TRAITS[[tr]][[1]] & env %in% TRAITS[[tr]][[2]]]
  d[, G := factor(TP_MO)]; d[, env := factor(env)]
  for (sc in c("raw", "ang")) {
    d[, y := get(sc)]
    nenv <- d[, .(n_env = .N), by = TP_MO]
    if (nlevels(d$env) == 1) {
      b <- d[, .(TP_MO, BLUE = y)]
      h2[[length(h2) + 1]] <- data.table(trait = tr, scale = sc, n_geno = nrow(b), n_env = 1, n_obs = nrow(d),
                                         Vg = NA, Venv = NA, Ve = NA, n_h = 1, H2_entry = NA, H2_Cullis = NA)
    } else {
      m  <- lmer(y ~ 0 + G + (1 | env), data = d, REML = TRUE)
      fe <- fixef(m); b <- data.table(TP_MO = sub("^G", "", names(fe)), BLUE = unname(fe))
      mr <- lmer(y ~ (1 | G) + (1 | env), data = d, REML = TRUE)
      vc <- as.data.frame(VarCorr(mr)); Vg <- vc$vcov[vc$grp == "G"]; Ve <- vc$vcov[vc$grp == "Residual"]
      Venv <- vc$vcov[vc$grp == "env"]
      nh <- 1 / mean(1 / nenv$n_env)
      pev <- attr(ranef(mr, condVar = TRUE)$G, "postVar")[1, 1, ]
      h2[[length(h2) + 1]] <- data.table(trait = tr, scale = sc, n_geno = nrow(b), n_env = nlevels(d$env), n_obs = nrow(d),
                                         Vg = Vg, Venv = Venv, Ve = Ve, n_h = nh,
                                         H2_entry = Vg / (Vg + Ve / nh), H2_Cullis = 1 - mean(2 * pev) / (2 * Vg))
    }
    b <- merge(b, nenv, by = "TP_MO")
    b[, `:=`(trait = tr, scale = sc)]
    blues[[length(blues) + 1]] <- b
  }
}
B <- rbindlist(blues); H <- rbindlist(h2)
B <- merge(B, idm[, .(TP_MO, collection, role, in_vcf)], by = "TP_MO")
# асимметрия и доля нулей
S <- B[, .(n = .N, n_genotyped = sum(role != "" & !grepl("^исключён", role)), mean = mean(BLUE), sd = sd(BLUE),
           skewness = skew(BLUE), share_zero_or_less = mean(BLUE <= 1e-8),
           shapiro_p = tryCatch(shapiro.test(BLUE)$p.value, error = function(e) NA)), by = .(trait, scale)]
H <- merge(H, S, by = c("trait", "scale"))
fwrite(B, file.path(OUT, "tables/01_blues.csv"))
fwrite(H, file.path(OUT, "tables/01_heritability_distribution.csv"))
print(H[, .(trait, scale, n_geno, n_genotyped, n_env, H2_entry = round(H2_entry, 2), H2_Cullis = round(H2_Cullis, 2),
            mean = round(mean, 3), skewness = round(skewness, 2), share_zero = round(share_zero_or_less, 2), shapiro_p = signif(shapiro_p, 2))])

# распределения
B[, lab := paste(trait, scale)]
p <- ggplot(B[!grepl("^исключён", role)], aes(BLUE)) + geom_histogram(bins = 30, fill = "grey40") +
  facet_wrap(~ lab, scales = "free", ncol = 4) + theme_bw(base_size = 9) +
  labs(x = "BLUE (raw: %, ang: радианы asin(sqrt(p)))", y = "число образцов",
       title = "Распределение BLUE морозостойкости (образцы с генотипами)")
ggsave(file.path(OUT, "figures/01_pheno_distributions.png"), p, width = 11, height = 7, dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_01_blue.txt"))
