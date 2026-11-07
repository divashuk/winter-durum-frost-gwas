# -*- coding: utf-8 -*-
"""Этап 2.0. Списки образцов потоков (только генотипированные И фенотипированные) + таблица переименования хромосом.
Файлы пишутся с LF (их читают bcftools/plink в WSL)."""
import pandas as pd, os
T = r"/work/project/gwas_frost/tables"
b = pd.read_csv(os.path.join(T, "01_blues.csv"))
idm = pd.read_csv(os.path.join(T, "01_id_map.csv"))
ok = set(idm[~idm.role.str.startswith("исключён")].TP_MO)
s1 = sorted(set(b[(b.trait == "FT_all")].TP_MO) & ok)
s2 = sorted((set(b[b.trait.isin(["FT_all", "ITFT_m12"])].TP_MO)) & ok)

def w(name, text):
    with open(os.path.join(T, name), "w", newline="\n") as f:
        f.write(text)

w("02_samples_stream1.txt", "\n".join(s1) + "\n")
w("02_samples_stream2.txt", "\n".join(s2) + "\n")
chrs = ["1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B", "5A", "5B", "6A", "6B", "7A", "7B"]
w("02_chr_rename.txt", "".join(f"CM1182{60 + i}.1\t{c}\n" for i, c in enumerate(chrs)))
print("stream1", len(s1), "stream2", len(s2))
