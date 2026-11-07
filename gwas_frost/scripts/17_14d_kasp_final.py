# -*- coding: utf-8 -*-
"""Раздел 17, шаг 7е. Окончательный выбор KASP — scripts/16_07e_kasp_final.py без изменения метода (посадки праймеров, посторонние
продукты, ожидаемые генотипы Langdon, Svevo, родителей F2) -> tables/17_kasp_primers.csv, tables/17_kasp_specificity.csv.
Строки доминантной ПЦР на делецию CBF-B (раздел 16) здесь не повторяются — они остаются в tables/16_kasp_primers.csv и входят в часть D."""
import os, sys, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
lib = importlib.import_module("17_lib")
src = lib.source("16_07e_kasp_final.py")
pcr = src[src.index("# доминантная ПЦР по делеции CBF-B12/B14/B15"):src.index("P = pd.DataFrame(out)")]
lib.run16("16_07e_kasp_final.py", [
    ('S16 = "/work/scratch/s16/"', f'S16 = "{lib.S17}"'),
    ('os.path.join(OUT, "tables/16_kasp_primers.csv")', 'os.path.join(OUT, "tables/17_kasp_primers.csv")'),
    ('os.path.join(OUT, "tables/16_kasp_specificity.csv")', 'os.path.join(OUT, "tables/17_kasp_specificity.csv")'),
    (pcr, ""),
    # тип столбца target: смесь bool Python и numpy давала object, и «~» не было логическим отрицанием (в 16 не проявилось)
    ('hx = hx.assign(target=is_t.values if len(hx) else [])', 'hx = hx.assign(target=(is_t.values if len(hx) else np.array([])).astype(bool))')])
