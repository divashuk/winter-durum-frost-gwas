import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description="Reproduce manuscript Figures 1 and 2 from association results.")
parser.add_argument("--data-dir", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
args = parser.parse_args()
D, OUT = args.data_dir, args.output_dir
OUT.mkdir(parents=True, exist_ok=True)
CHR = [f"{i}{g}" for i in range(1, 8) for g in "AB"]
MEFF, BONF = 4.54, 5.78
INK, MUTED, GRID = "#1f2328", "#6b7280", "#e5e7eb"
DOT = ["#9aa3ad", "#3b5b7a"]  # чередование хромосом: два нейтральных тона
plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK})


def load(name):
    return pd.read_csv(D / name.removesuffix(".gz"))


def manhattan(ax, d, title, labels):
    d = d[d.chr.isin(CHR)].copy()
    d["mlp"] = -np.log10(d.p)
    off, ticks, x0 = {}, [], 0
    for c in CHR:
        L = d.loc[d.chr == c, "pos"].max()
        off[c] = x0
        ticks.append(x0 + L / 2)
        x0 += L + 20e6
    d["x"] = d.pos + d.chr.map(off)
    for i, c in enumerate(CHR):
        s = d[d.chr == c]
        ax.scatter(s.x / 1e6, s.mlp, s=3, c=DOT[i % 2], linewidths=0, rasterized=True)
    ax.axhline(BONF, color="#b42318", lw=0.8, ls="--")
    ax.axhline(MEFF, color="#b54708", lw=0.8, ls=":")
    ax.text(x0 / 1e6, BONF, " Bonferroni", va="center", ha="left", fontsize=7, color="#b42318")
    ax.text(x0 / 1e6, MEFF, " Meff", va="center", ha="left", fontsize=7, color="#b54708")
    for lab, c, pos, dx in labels:
        s = d[(d.chr == c) & (abs(d.pos - pos) < 15e6)]
        r = s.loc[s.mlp.idxmax()]
        ax.annotate(lab, (r.x / 1e6, r.mlp), xytext=(dx, 6), textcoords="offset points", ha="center", fontsize=7, color=INK)
    ax.set_xticks(np.array(ticks) / 1e6)
    ax.set_xticklabels(CHR)
    ax.set_xlim(-20, x0 / 1e6 + 120)
    ax.set_ylim(0, 7.4)
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.set_title(title, loc="left", fontsize=8.5, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def fig1():
    fig, axs = plt.subplots(2, 1, figsize=(7.2, 4.6), constrained_layout=True)
    manhattan(axs[0], load("g09_S1_FTall_A.csv.gz"),
              "(a) Frost tolerance across all 12 environments (FT_all), n = 144",
              [("Fr-A2 region", "5A", 514.68e6, 0), ("1B", "1B", 653.8e6, 0), ("4B", "4B", 28.71e6, 8), ("4A", "4A", 631.7e6, -8)])
    manhattan(axs[1], load("g09_RU_FTm12_A.csv.gz"),
              "(b) Frost tolerance at −12 °C (FT_m12), n = 91",
              [("Fr-A2 region", "5A", 514.68e6, 0), ("4A", "4A", 642.5e6, 0), ("3A", "3A", 0.45e6, 0)])
    fig.savefig(OUT / "Figure1.jpg", dpi=600, pil_kwargs={"quality": 100, "subsampling": 0, "optimize": True})


LOCI = [("Fr-A2 region, 5A 475–566 Mb", "5A", 475e6, 566.2e6),
        ("4A 618–650 Mb", "4A", 617.9e6, 650.3e6),
        ("3A 0–14 Mb", "3A", 0, 14.3e6),
        ("4B 28.2–30.2 Mb", "4B", 28.2e6, 30.3e6),
        ("1B 650–661 Mb", "1B", 650e6, 661e6)]
REG = [("FT_all_ang", "All\n(n = 144)"), ("FT_m12_ang", "−12 °C\n(n = 91)"), ("FT_m15_ang", "−15 °C\n(n = 76)"),
       ("FT_m16_ang", "−16 °C\n(n = 75)"), ("FT_m17_ang", "−17 °C\n(n = 70)")]


def regime_matrix():
    g = load("g04_s1_GEMMA_A_ang.csv.gz")
    M = np.zeros((len(LOCI), len(REG)))
    for i, (_, c, a, b) in enumerate(LOCI):
        for j, (t, _) in enumerate(REG):
            s = g[(g.trait == t) & (g.chr == c) & (g.pos >= a) & (g.pos <= b)]
            M[i, j] = -np.log10(s.p.min())
    return M


def fig2():
    M = regime_matrix()
    pd.DataFrame(M, index=[l[0] for l in LOCI], columns=[r[0] for r in REG]).round(2).to_csv(OUT / "Figure2_values.csv")
    fig, ax = plt.subplots(figsize=(5.2, 2.6), constrained_layout=True)
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("blue", ["#f3f6fa", "#c9d8ea", "#7ea3cc", "#2f5d8f", "#173a63"])
    im = ax.imshow(M, cmap=cmap, vmin=0, vmax=6.5, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7.5,
                    color="white" if v > 3.6 else INK, fontweight="bold" if v >= MEFF else "normal")
            ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False, ec="white", lw=2))
    ax.set_xticks(range(len(REG)))
    ax.set_xticklabels([r[1] for r in REG])
    ax.set_yticks(range(len(LOCI)))
    ax.set_yticklabels([l[0] for l in LOCI])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label(r"best $-\log_{10}(p)$ in interval", fontsize=7)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=7, length=0)
    fig.savefig(OUT / "Figure2.jpg", dpi=600, pil_kwargs={"quality": 100, "subsampling": 0, "optimize": True})


if __name__ == "__main__":
    fig1()
    fig2()
    print("ok")
