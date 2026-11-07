# -*- coding: utf-8 -*-
"""Этап 6. Итоговые таблицы и сводные графики.
  tables/06_qtl_summary.csv          — одна строка на QTL (воспроизводимость, известные локусы, вердикт A/B)
  tables/06_qtl_hits_main.csv        — признак × поток × схема × модель: ведущий SNP, хр., позиция, MAF, эффект, p, PVE
  tables/06_excluded_unmatched_samples.csv — исключённые/несопоставленные/помеченные образцы с причинами
  figures/06_qtl_map.png             — карта QTL на хромосомах Langdon + известные локусы
  figures/06_schemeA_vs_B_heatmap.png — min -log10 p в интервалах воспроизводимых QTL: схема A против B
Входы читаются через \\wsl$ только для полных результатов не нужно — всё из tables/."""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
B = r"/work/project/gwas_frost"
T, F = os.path.join(B, "tables"), os.path.join(B, "figures")
INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e2dc"
BLUE, GREY = "#2a78d6", "#b9b7ae"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})

q = pd.read_csv(os.path.join(T, "05_qtl_intervals.csv"))
h = pd.read_csv(os.path.join(T, "05_significant_and_suggestive_snps.csv"))
ab = pd.read_csv(os.path.join(T, "05_schemeA_vs_B.csv"))
kp = pd.read_csv(os.path.join(T, "05_known_loci_positions.csv"))
genes = pd.read_csv(os.path.join(T, "05_genes_in_qtl.csv"))

# ---------- 06_qtl_hits_main ----------
s = h[h.level == "значимый"].copy()
best = s.sort_values("p").groupby(["QTL", "stream", "trait", "scheme", "model"], as_index=False).first()
best["minus_log10_p"] = -np.log10(best.p)
best = best.merge(q[["QTL", "reliability", "n_models", "n_streams", "n_traits"]], on="QTL", how="left")
cols = ["QTL", "reliability", "trait", "stream", "scheme", "model", "n", "snp_id", "chr", "pos", "MAF", "effect", "effect_allele",
        "p", "minus_log10_p", "pve", "sig_bonf", "sig_meff", "sig_fdr", "rep_models", "rep_streams", "rep_traits", "n_models", "n_streams", "n_traits"]
best[cols].sort_values(["QTL", "stream", "trait", "scheme", "model"]).to_csv(os.path.join(T, "06_qtl_hits_main.csv"), index=False, encoding="utf-8-sig")

# ---------- 06_qtl_summary ----------
def verdict_summary(Q):
    d = ab[ab.QTL == Q]
    if not len(d):
        return ""
    c = d.verdict.value_counts()
    return "; ".join(f"{k}: {v}" for k, v in c.items())
q["A_vs_B (число сочетаний поток×признак×модель)"] = q.QTL.map(verdict_summary)
q["n_genes_svevo"] = q.QTL.map(genes.groupby("QTL").size()).fillna(0).astype(int)
q["length_Mb"] = ((q.end - q.start) / 1e6).round(1)
q["min_minus_log10_p"] = (-np.log10(q.min_p)).round(2)
q.to_csv(os.path.join(T, "06_qtl_summary.csv"), index=False, encoding="utf-8-sig")

# ---------- исключённые / несопоставленные / помеченные ----------
idm = pd.read_csv(os.path.join(T, "01_id_map.csv"))
het = pd.read_csv(os.path.join(T, "02_sample_qc_het_missing.csv"))
ibs = pd.read_csv(os.path.join(T, "02_ibs_pairs_ge098.csv"))
mc = pd.read_csv(os.path.join(T, "01_manual_check.csv"))
rows = []
for _, r in idm[idm.role.str.startswith("исключён")].iterrows():
    rows.append(dict(TP_MO=r.TP_MO, name=r.name_pheno, collection=r.collection, status="исключён", stage="1 (ID)", reason=r.role.replace("исключён: ", "")))
for _, r in het[het.het_class != "норма (<=0.125)"].iterrows():
    rows.append(dict(TP_MO=r.IID, name=r.name_pheno, collection=r.collection, status="помечен (оставлен)", stage=f"2 (QC, {r.stream})",
                     reason=f"гетерозиготность {r.het_rate:.3f}: {r.het_class}"))
for _, r in ibs.iterrows():
    for a, b, na, nb in [(r.id1, r.id2, r.name1, r.name2), (r.id2, r.id1, r.name2, r.name1)]:
        rows.append(dict(TP_MO=a, name=na, collection="", status="помечен (оставлен)", stage=f"2 (QC, {r.stream})",
                         reason=f"IBS {r.IBS:.4f} с {b} ({nb}); разные названия, синонимия не установлена"))
for _, r in mc.iterrows():
    rows.append(dict(TP_MO=r.TP_MO, name=r.names, collection="", status="для ручной проверки", stage="1 (ID)",
                     reason=f"{r.reason}; в VCF: {r.in_vcf}; в GWAS: {r.used_in_gwas}"))
ex = pd.DataFrame(rows)
ex["collection"] = ex.apply(lambda r: r.collection or idm.set_index("TP_MO").collection.get(r.TP_MO.split(";")[0], ""), axis=1)
ex.to_csv(os.path.join(T, "06_excluded_unmatched_samples.csv"), index=False, encoding="utf-8-sig")
print("исключено/помечено строк:", len(ex), ex.status.value_counts().to_dict())

# ---------- карта QTL ----------
CHR = ["1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"]
LEN = dict(zip(CHR, [602.5, 712.9, 803.5, 823.4, 763.0, 863.8, 767.9, 697.5, 722.6, 739.0, 631.9, 744.9, 754.1, 751.3]))
genes_known = kp[~kp.locus.str.contains(r"\(LD\)") & kp.chr.notna()].copy()
genes_known["label"] = genes_known.locus.str.replace(r" → .*", "", regex=True).str.replace("_BAC", "", regex=False)
genes_known.loc[genes_known.locus.str.contains("→ Fr-A2"), "label"] = "Fr-A2 (CBF)"
genes_known.loc[genes_known.locus.str.contains("→ Fr-B2"), "label"] = "Fr-B2 (CBF)"
genes_known = genes_known.drop_duplicates(["label", "chr"])
genes_known["is_kasp"] = genes_known.source.eq("KASP")
fig, ax = plt.subplots(figsize=(11, 12))
SP = 1.7   # расстояние между хромосомами (место под подписи)
for i, c in enumerate(CHR):
    y = (len(CHR) - 1 - i) * SP
    last_lab = -1e9; lev = 0
    ax.plot([0, LEN[c]], [y, y], color=GRID, lw=6, solid_capstyle="round", zorder=1)
    for _, r in q[q.chr == c].iterrows():
        rel = r.reliability == "воспроизводимый"
        ax.add_patch(plt.Rectangle((r.start / 1e6, y - 0.28), max((r.end - r.start) / 1e6, 1.5), 0.56,
                                   color=BLUE if rel else GREY, zorder=3, lw=0))
        if rel:
            xm = (r.start + r.end) / 2e6
            lev = lev + 1 if xm - last_lab < 25 else 0   # близкие подписи — ступенькой выше
            ax.text(xm, y + 0.36 + 0.3 * lev, r.QTL, ha="center", va="bottom", fontsize=7.5, color=INK)
            last_lab = xm
    for _, r in genes_known[genes_known.chr == c].iterrows():
        x = (r.start + r.end) / 2e6
        if r.is_kasp:   # KASP-маркеры — пунктирный штрих без подписи (см. легенду и таблицу 05_known_loci_positions)
            ax.plot([x, x], [y - 0.3, y + 0.3], color=INK2, lw=0.9, ls=(0, (2, 1.5)), zorder=4)
            continue
        ax.plot([x, x], [y - 0.42, y + 0.42], color=INK, lw=1.2, zorder=4)
        ax.text(x, y - 0.48, r.label, ha="center", va="top", fontsize=7, color=INK2)
ax.set_yticks([k * SP for k in range(len(CHR))]); ax.set_yticklabels(CHR[::-1])
ax.set_xlabel("Позиция на Langdon (GCA_051168695.1), Мб"); ax.set_xlim(-5, 880); ax.set_ylim(-1.2, (len(CHR) - 1) * SP + 1.2)
for sp in ["top", "right", "left"]:
    ax.spines[sp].set_visible(False)
ax.tick_params(axis="y", length=0)
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
ax.legend(handles=[Patch(color=BLUE, label="QTL воспроизводимый (≥2 модели вкл. MLM/GEMMA, ≥2 потока или ≥2 признака)"),
                   Patch(color=GREY, label="QTL единичный (одна модель / один признак)"),
                   Line2D([0], [0], color=INK, lw=1.2, label="известный ген (выравнивание последовательности NCBI на Langdon)"),
                   Line2D([0], [0], color=INK2, lw=0.9, ls=(0, (2, 1.5)), label="KASP-маркер durum_frost (по праймерам)")],
          loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=2, frameon=False, fontsize=8)
ax.set_title("QTL морозостойкости (оба потока, все модели и схемы) и известные локусы", loc="left", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(F, "06_qtl_map.png"), dpi=200); plt.close(fig)

# ---------- тепловая карта A против B ----------
rel = q[q.reliability == "воспроизводимый"].QTL.tolist()
full = []
for Q in rel:
    d = ab[ab.QTL == Q]
    for _, r in d.iterrows():
        full.append(dict(QTL=Q, col=f"{r.stream}\n{r.trait}", A=-np.log10(r.minp_A), B=-np.log10(r.minp_B)))
fd = pd.DataFrame(full)
COLS = [("stream1", "FT_all_ang"), ("stream1", "FT_m12_ang"), ("stream1", "FT_m15_ang"), ("stream2", "FT12j_ang"), ("stream2", "FTcol_ang")]
# min p по всем моделям внутри интервала; берём из полной таблицы A/B (там только сочетания со значимостью хотя бы в одной схеме)
M = np.full((len(rel), len(COLS) * 2), np.nan)
for i, Q in enumerate(rel):
    for j, (s_, t_) in enumerate(COLS):
        d = ab[(ab.QTL == Q) & (ab.stream == s_) & ab.trait.isin([t_, t_.replace("_ang", "_raw")])]   # минимум по обеим шкалам
        if len(d):
            M[i, 2 * j] = -np.log10(d.minp_A.min()); M[i, 2 * j + 1] = -np.log10(d.minp_B.min())
fig, ax = plt.subplots(figsize=(10, 6))
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("blue", ["#f4f8fd", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])
im = ax.imshow(np.nan_to_num(M, nan=0), cmap=cmap, vmin=0, vmax=12, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        ax.text(j, i, "—" if np.isnan(v) else f"{v:.1f}", ha="center", va="center", fontsize=7.5,
                color=INK2 if np.isnan(v) or v < 6 else "white")
ax.set_yticks(range(len(rel))); ax.set_yticklabels(rel)
ax.set_xticks(range(M.shape[1])); ax.set_xticklabels([f"{t}\n{s}\n{sch}" for s, t in COLS for sch in ("A", "B")], fontsize=7)
for j in range(1, len(COLS)):
    ax.axvline(2 * j - 0.5, color="white", lw=3)
cb = fig.colorbar(im, ax=ax, fraction=0.03); cb.set_label("min −log10 p в интервале (по 4 моделям)")
ax.set_title("Воспроизводимые QTL: схема A (PC + родство) против B (+ маркеры durum_frost)\n"
             "min по 4 моделям и 2 шкалам (угловая, %); «—» — в этом сочетании нет значимости ни в одной схеме", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(F, "06_schemeA_vs_B_heatmap.png"), dpi=200); plt.close(fig)
print(q[["QTL", "chr", "start", "end", "reliability", "known_overlap", "A_vs_B (число сочетаний поток×признак×модель)"]].to_string())
