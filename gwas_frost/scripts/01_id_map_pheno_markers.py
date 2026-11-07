# -*- coding: utf-8 -*-
"""Этап 1. Таблица соответствия ID, длинная таблица фенотипов, маркеры-ковариаты.

Источники (только чтение):
  durum_frost/00_raw/абоминация_но_по_тз.xlsx       — ключ «Внешний шифр» N_TP_MO ↔ «Номер образца (ТП)», сырые вызовы, итал. фенотипы
  Oz_tv_sic/Sequncing runs database (1).xlsx         — TP_MO_xxx ↔ название (выгрузка tables/seqdb_TP_MO_463.csv)
  durum_frost/08_deposit/data/genotypes_161x11.csv   — финальные KASP 161 × 11 (= Табл. S1)
  durum_frost/08_deposit/data/phenotype_long.csv     — финальные фенотипы статьи (по средам)
  durum_frost/01_data/markers/genotypes_clean_v2.csv — KASP итальянской коллекции (первая панель)
  durum_frost/00_raw/ТП. Италия по годам.13.06.23.xlsx — итал. промораживание −12 °C, сезон 2019–20
  durum_frost/14_submission_final_2026-09-15/03_supplementary/Supplementary_Material.docx — Табл. S4 (происхождение)
Выход (gwas_frost/tables):
  01_id_map.csv, 01_manual_check.csv, 01_not_matched.csv,
  01_pheno_long.csv, 01_markers_012.csv, 01_marker_coding.csv
"""
import os, re, sys, gzip
import numpy as np, pandas as pd, docx

sys.stdout.reconfigure(encoding="utf-8")
DF = r"/work/phenotypes"
PJ = r"/work/project"
T = os.path.join(PJ, "gwas_frost", "tables")
ab = pd.read_excel(os.path.join(DF, "00_raw", "абоминация_но_по_тз.xlsx"))
sq = pd.read_csv(os.path.join(T, "seqdb_TP_MO_463.csv"))
g161 = pd.read_csv(os.path.join(DF, "08_deposit", "data", "genotypes_161x11.csv"))
ph = pd.read_csv(os.path.join(DF, "08_deposit", "data", "phenotype_long.csv"))
gv2 = pd.read_csv(os.path.join(DF, "01_data", "markers", "genotypes_clean_v2.csv"))

with gzip.open(os.path.join(PJ, "FINAL-181124.vcf.gz"), "rt") as f:
    for l in f:
        if l.startswith("#CHROM"):
            vcf = set(l.rstrip().split("\t")[9:]); break

# ---------------- Табл. S4 (происхождение) ----------------
d = docx.Document(os.path.join(DF, "14_submission_final_2026-09-15", "03_supplementary", "Supplementary_Material.docx"))
t4 = [t for t in d.tables if len(t.columns) >= 3 and t.rows[0].cells[2].text.strip() == "Originator"][0]
rows = [[c.text.strip() for c in r.cells] for r in t4.rows]
s4 = pd.DataFrame(rows[1:], columns=["sample_id", "name_s4", "originator", "country", "year_GSI", "origin_group"])
s4["sample_id"] = s4.sample_id.astype(int)

# ---------------- 1. Таблица соответствия ----------------
ab["TP_MO"] = ab["Внешний шифр"].str.extract(r"^(\d+)_TP_MO")[0].astype(int).map(lambda n: f"TP_MO_{n:03d}")
idm = ab[["TP_MO", "Номер образца (ТП)", "Название"]].rename(columns={"Номер образца (ТП)": "sample_id_TP", "Название": "name_pheno"})
idm = idm.merge(sq.rename(columns={"ID образца": "TP_MO", "Образец": "name_seqdb"})[["TP_MO", "name_seqdb"]], on="TP_MO", how="left")
idm["collection"] = np.where(ab["п/п"] <= 163, "RU_winter", "italy")
idm["in_vcf"] = idm.TP_MO.isin(vcf)
art = g161[["sample_id", "name_ru", "name_en"]].rename(columns={"name_ru": "name_article_ru", "name_en": "name_article_en"})
idm["in_article"] = (idm.collection == "RU_winter") & idm.sample_id_TP.isin(g161.sample_id)
idm["_key"] = np.where(idm.in_article, idm.sample_id_TP, -1)
idm = idm.merge(art.rename(columns={"sample_id": "_key"}), on="_key", how="left")
idm = idm.merge(s4[["sample_id", "originator", "country", "year_GSI", "origin_group"]].rename(columns={"sample_id": "_key"}),
                on="_key", how="left").drop(columns="_key")
idm.loc[idm.collection == "italy", "origin_group"] = "italian collection"
# происхождение для графиков: Россия / Италия / прочее
def orig3(r):
    if r.collection == "italy":
        return "Италия (итал. коллекция)"
    if isinstance(r.country, str) and r.country == "Russia":
        return "Россия"
    if isinstance(r.country, str) and r.country == "Italy":
        return "Италия"
    return "прочее" if isinstance(r.country, str) else "нет данных"
idm["origin3"] = idm.apply(orig3, axis=1)
# тип развития: статья — все озимые; итальянские поимённо неизвестны
vrn = gv2.set_index(gv2.collection.astype(str) + "_" + gv2.sample_id.astype(str))["VrnA1_9K0001"]
idm["VrnA1_9K0001"] = [vrn.get(("winter_" if c == "RU_winter" else "italy_") + str(s)) for c, s in zip(idm.collection, idm.sample_id_TP)]
idm["growth_habit"] = np.where(idm.collection == "RU_winter", "озимый (статья: VRN-A1, VRN-B1 озимые)", "неизвестен поимённо (панель Habib 2026: преимущественно яровые)")
# совпадение названий
norm = lambda s: re.sub(r"[\s\-\.\(\)/,]", "", str(s).lower())
idm["name_match_seqdb"] = idm.name_pheno.map(norm) == idm.name_seqdb.map(norm)
idm["name_match_article"] = np.where(idm.in_article, idm.name_pheno.map(norm) == idm.name_article_ru.map(norm), np.nan)
# итог: роль в анализе
def role(r):
    if not r.in_vcf:
        return "исключён: нет в VCF"
    if r.collection == "RU_winter" and not r.in_article:
        return "исключён: не входит в 161 статьи (объединён с «Синьора»)"
    return "поток1+поток2" if r.collection == "RU_winter" else "поток2"
idm["role"] = idm.apply(role, axis=1)
idm.to_csv(os.path.join(T, "01_id_map.csv"), index=False, encoding="utf-8-sig")

# несопоставленные / исключённые на этапе ID
idm[idm.role.str.startswith("исключён")][["TP_MO", "sample_id_TP", "name_pheno", "collection", "role"]].to_csv(
    os.path.join(T, "01_not_matched.csv"), index=False, encoding="utf-8-sig")

# ---------------- ручная проверка: возможные дубликаты названий среди всех 463 TP_MO ----------------
TR = dict(zip("абвгдеёжзийклмнопрстуфхцчшщъыьэюя",
              ["a","b","v","g","d","e","e","zh","z","i","i","k","l","m","n","o","p","r","s","t","u","f","h","c","ch","sh","sch","","y","","e","yu","ya"]))
def key(s):
    s = str(s).lower()
    s = re.sub(r"\(.*?\)", " ", s)
    s = "".join(TR.get(ch, ch) for ch in s)
    s = s.replace("h", "").replace("j", "i").replace("w", "v")
    s = re.sub(r"^(ama|ama)\s*", "", s)
    return re.sub(r"[^a-z0-9]", "", s)
allq = sq.rename(columns={"ID образца": "TP_MO", "Образец": "name"})[["TP_MO", "name"]].copy()
allq["k"] = allq.name.map(key)
paren = allq.name.str.extract(r"\((.*?)\)")[0].map(lambda x: key(x) if isinstance(x, str) else None)
allq["k_paren"] = paren
mc = []
for k, grp in allq.groupby("k"):
    if len(grp) > 1 and k:
        mc.append(dict(reason="одинаковое имя после нормализации/транслита", TP_MO=";".join(grp.TP_MO), names=";".join(grp.name)))
for _, r in allq[allq.k_paren.notna()].iterrows():
    hit = allq[(allq.k == r.k_paren) & (allq.TP_MO != r.TP_MO)]
    if len(hit):
        mc.append(dict(reason="название в скобках совпадает с другим образцом", TP_MO=";".join([r.TP_MO] + hit.TP_MO.tolist()),
                       names=";".join([r["name"]] + hit.name.tolist())))
# подстрока-код линии (напр. «15908ds1122-07» и «AMA 15908ds1122-07»)
codes = allq[allq.k.str.len() >= 7]
for _, r in codes.iterrows():
    hit = allq[(allq.k != r.k) & allq.k.str.contains(re.escape(r.k)) & (allq.TP_MO != r.TP_MO)]
    if len(hit):
        mc.append(dict(reason="код линии входит в другое название", TP_MO=";".join([r.TP_MO] + hit.TP_MO.tolist()),
                       names=";".join([r["name"]] + hit.name.tolist())))
mc = pd.DataFrame(mc).drop_duplicates("TP_MO")
mc["in_vcf"] = mc.TP_MO.map(lambda s: ";".join(str(x in vcf) for x in s.split(";")))
mc["used_in_gwas"] = mc.TP_MO.map(lambda s: ";".join(str(x in set(idm[~idm.role.str.startswith("исключён")].TP_MO)) for x in s.split(";")))
mc.to_csv(os.path.join(T, "01_manual_check.csv"), index=False, encoding="utf-8-sig")

# ---------------- 2. Фенотипы (длинный формат) ----------------
fr = ph[ph.trait_en.str.startswith("frost")].copy()
fr["regime"] = fr.regime.replace({"15_12": "-12"})
fr = fr.merge(idm[idm.collection == "RU_winter"][["TP_MO", "sample_id_TP"]], left_on="sample_id", right_on="sample_id_TP")
fr = fr.assign(collection="RU_winter", env=lambda x: "RU_" + x.year.astype(str) + "_" + x.regime.astype(str),
               source="08_deposit/data/phenotype_long.csv", use="да")[["TP_MO", "collection", "env", "year", "regime", "value", "source", "use"]]
it_ab = ab[ab["п/п"] > 163].copy()
it_ab["TP_MO"] = it_ab["TP_MO"]
c21 = "Промораживание, % живых растений при -12_италия_2021_занятой пар_15_12"
c22 = "Морозостойкость, %_италия_2022_занятой пар_15_12"
it21 = pd.DataFrame(dict(TP_MO=it_ab.TP_MO, value=pd.to_numeric(it_ab[c21], errors="coerce"), year=2021,
                         source="абоминация: " + c21, use="да"))
it22 = pd.DataFrame(dict(TP_MO=it_ab.TP_MO, value=pd.to_numeric(it_ab[c22], errors="coerce"), year=2022,
                         source="абоминация: " + c22, use="нет: источник не найден, r≈0.03–0.10 с сезонами 2020/2021"))
x = pd.read_excel(os.path.join(DF, "00_raw", "ТП. Италия по годам.13.06.23.xlsx"), sheet_name="2019-20", header=None)
s1 = x.iloc[6:, [3, 4, 7]].copy(); s1.columns = ["tp", "name", "value"]
s1["tp"] = pd.to_numeric(s1.tp, errors="coerce"); s1["value"] = pd.to_numeric(s1.value, errors="coerce")
s1 = s1.dropna(subset=["tp"])
chk = s1.merge(it_ab[["Номер образца (ТП)", "Название", "TP_MO"]], left_on="tp", right_on="Номер образца (ТП)")
okname = chk.name.map(norm).str[:5] == chk["Название"].map(norm).str[:5]
print(f"Италия 2019–20: строк с ТП из коллекции {len(chk)}, название совпадает {okname.sum()}")
it20 = pd.DataFrame(dict(TP_MO=chk.TP_MO, value=chk.value, year=2020,
                         source="ТП. Италия по годам.13.06.23.xlsx, лист 2019-20", use=np.where(okname, "да", "нет: название не совпало")))
itp = pd.concat([it20, it21, it22]).dropna(subset=["value"])
itp = itp.assign(collection="italy", regime="-12", env=lambda x: "IT_" + x.year.astype(str) + "_-12")
pl = pd.concat([fr, itp[fr.columns]], ignore_index=True)
pl.to_csv(os.path.join(T, "01_pheno_long.csv"), index=False, encoding="utf-8-sig")
print(pl.groupby(["collection", "env", "use"]).value.agg(["size", "mean", lambda v: (v == 0).mean()]).round(2).to_string())

# ---------------- 3. Маркеры 0/1/2 ----------------
# Кодировка: 2 = гомозигота по указанному классу (ALLELE_2), 0 = гомозигота по другому классу, 1 = гетерозигота;
# heterogeneous / нет вызова -> NA. Для самоопылителя классы KASP практически всегда гомозиготные.
CODE = {  # marker: (класс, кодируемый как 2, описание)
    "S2269949": ("FAM", "Fr-A2 CNV-тег; FAM (A) — мажорный в РФ (благоприятный по статье)"),
    "S1077313": ("HEX", "Fr-A2 CNV-тег; HEX (A) — мажорный в РФ"),
    "S1862541": ("HEX", "Fr-A2 CNV-тег; HEX (C) — мажорный в РФ"),
    "S1298957": ("HEX", "Fr-A2 CNV-тег; HEX — мажорный в РФ"),
    "S1051014": ("FAM", "Fr-A2 не-CNV; FAM = C"),
    "Vrn1_new": ("2147", "VRN-A1 экзон 4; 2 = vrn-A1b (длинная яровизация)"),
    "Mor534": ("G", "5B_517276534; 2 = G (благоприятный в РФ)"),
    "Mor587": ("T", "1B_41099587; 2 = T (благоприятный в РФ)"),
    "Mor866": ("A", "7B_598228866; 2 = A"),
}
def code(v, cls):
    s = str(v).strip()
    if s in ("nan", "") or "heterogen" in s.lower():
        return np.nan
    if "heterozyg" in s.lower():
        return 1.0
    return 2.0 if s.upper().startswith(cls.upper()) else 0.0
ru = idm[idm.collection == "RU_winter"][["TP_MO", "sample_id_TP"]].merge(g161, left_on="sample_id_TP", right_on="sample_id")
itg = idm[idm.collection == "italy"][["TP_MO", "sample_id_TP"]].merge(gv2[gv2.collection == "italy"], left_on="sample_id_TP", right_on="sample_id")
out = []
for df, coll in [(ru, "RU_winter"), (itg, "italy")]:
    m = pd.DataFrame({"TP_MO": df.TP_MO, "collection": coll})
    for mk, (cls, _) in CODE.items():
        m[mk] = df[mk].map(lambda v: code(v, cls)) if mk in df.columns else np.nan
    out.append(m)
mk = pd.concat(out, ignore_index=True)
mk.to_csv(os.path.join(T, "01_markers_012.csv"), index=False, encoding="utf-8-sig")
pd.DataFrame([(k, v[0], v[1]) for k, v in CODE.items()], columns=["marker", "class_coded_2", "note"]).to_csv(
    os.path.join(T, "01_marker_coding.csv"), index=False, encoding="utf-8-sig")
print(mk.groupby("collection")[list(CODE)].agg(lambda s: f"{(s==2).sum()}/{(s==0).sum()}/{(s==1).sum()}/NA{s.isna().sum()}").T.to_string())

print("\nID map:", idm.role.value_counts().to_dict())
print("Названия совпали (фенотип vs реестр):", int(idm.name_match_seqdb.sum()), "из", len(idm))
print("Ручная проверка, групп:", len(mc))
print(mc.to_string())
