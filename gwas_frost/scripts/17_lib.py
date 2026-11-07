# -*- coding: utf-8 -*-
"""Раздел 17. Запуск скриптов раздела 16 (и 14) без изменения метода: исходный текст скрипта исполняется с заменой путей
(s16 -> s17, genes16 -> genes17, tables/16_* -> tables/17_* и т. п.). Каждая замена обязана найтись в тексте (иначе ошибка),
поэтому расхождение версий скрипта и замен не пройдёт молча. Список замен задаётся в вызывающем скрипте 17_*.
variants17() — модуль 16_06_variants.py с путями раздела 17 (для скриптов, которые его импортируют: 16_05b, 16_07b)."""
import os, sys, types

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
UNC = "/work/scratch/"
WK = UNC if os.name == "nt" else "/work/scratch/"      # запуск из Windows или из WSL
S16, S17 = WK + "s16/", WK + "s17/"
LOAD_V = 'spec = importlib.util.spec_from_file_location("v", os.path.join(HERE, "16_06_variants.py")); v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)'
VAR_SUBS = [('S16 = "/work/scratch/s16/"', f'S16 = "{S17}"'),
            ('ALN = os.path.join(OUT, "genes16", "aln")', 'ALN = os.path.join(OUT, "genes17", "aln")'),
            ('os.path.join(OUT, "tables/16_gene_copies_by_genome.csv")', 'os.path.join(OUT, "tables/17_gene_copies_by_genome.csv")'),
            ('os.path.join(OUT, "tables/16_variants_functional.csv")', 'os.path.join(OUT, "tables/17_variants_functional.csv")'),
            ('os.path.join(OUT, "tables/16_haplotype_groups.csv")', 'os.path.join(OUT, "tables/17_haplotype_groups.csv")')]


def source(script, subs=(), cut=None):
    src = open(os.path.join(HERE, script), encoding="utf-8").read()
    if cut: src = src[:src.index(cut)]
    for a, b in subs:
        if a not in src: raise ValueError(f"{script}: не найдено «{a}»")
        src = src.replace(a, b)
    if WK != UNC: src = src.replace(UNC, WK)          # из WSL — пути WSL
    return src


def run16(script, subs=(), argv=None, cut=None, name="__main__", inject=None):
    """исполнить скрипт с заменами; argv — аргументы командной строки; cut — исполнить текст только до этой строки (библиотека);
    inject — глобальные имена, подставляемые до исполнения (например, модуль вариантов с путями 17)"""
    src = source(script, subs, cut)
    old = sys.argv
    sys.argv = [os.path.join(HERE, script)] + list(argv or [])
    ns = {"__name__": name, "__file__": os.path.join(HERE, script), **(inject or {})}
    try:
        exec(compile(src, script, "exec"), ns)
    finally:
        sys.argv = old
    return ns


def variants17():
    return types.SimpleNamespace(**run16("16_06_variants.py", VAR_SUBS, name="v"))
