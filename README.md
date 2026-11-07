# Winter durum wheat frost-tolerance analysis

Analysis source for **Association Mapping of Frost Tolerance in a Winter Durum Wheat Gene Pool with a Nearly Fixed Fr-A2 Haplotype**.

Version: **1.0.0**. Release date: **2026-11-07**.

This repository preserves the analysis scripts associated with the manuscript. It includes phenotype preparation, genotype quality control and imputation, population structure, association and conditional analyses, candidate-gene comparisons, primer design, protein analyses, and manuscript figure generation.

## Data availability

This is a **code-only release**. Individual-level phenotypes and genotypes, association-result tables, reference genome files, model weights, and generated analysis outputs are not included. The underlying phenotypic and genotypic data and additional analysis outputs are available from the corresponding author, Mikhail Divashuk (divashuk@gmail.com), upon reasonable request. Summary results are provided in the article and its supplementary tables.

## Contents and scope

| Location | Content |
| --- | --- |
| `gwas_frost/scripts/00_*`–`07_*` | Environment, phenotypes and BLUEs, QC, structure, GWAS and initial interpretation |
| `gwas_frost/scripts/09_*`–`14_*` | Origin-adjusted and conditional analyses, LD, locus interpretation and candidate genes |
| `gwas_frost/scripts/16_*`, `17_*` | Genome comparisons, sequence variants, KASP design and protein analyses |
| `figures/make_figures.py` | Final manuscript Figures 1 and 2 as 600 dpi JPEG files |
| `environment/` | Software versions and R session information recorded during the original analyses |
| `configure_paths.py` | Create a separate working copy with installation-specific paths |
| `SOURCE_MANIFEST.json` | Source and release SHA-256 checksums for the 112 analysis scripts |

The manuscript uses the Russian winter durum panel (`stream1`, `S1` and `RU` in the original scripts). Some shared scripts retain branches for other panels or parental materials used during the wider project. Those branches are retained to preserve source dependencies; their inclusion does not mean their results are reported in this manuscript. Inspect the script arguments and inputs before running a stage. The files are an analysis-source archive, not a single-command workflow.

## Environment and paths

The original workflow ran under Ubuntu/WSL with R, Python and command-line bioinformatics tools. Recorded software and R package versions are in `environment/`; they describe the original analysis environment, not a newly resolved dependency lockfile. Sequence and protein stages also require their stated external tools and model weights, which are not redistributed here.

Only installation-specific paths and a compute-host alias were normalized in the archived analysis scripts. Statistical formulas, thresholds and analytical parameters were retained. The figure script now accepts explicit input and output directories. `SOURCE_MANIFEST.json` records the original and normalized source hashes.

The default layout is `/work/project` (project), `/work/scratch` (large intermediate files), `/work/phenotypes` (source phenotype project), `/work/remote` (remote assembly work), and `/opt/micromamba` (environments). To use different directories, make a separate configured copy, for example:

```bash
python configure_paths.py --destination /srv/durum-code \
  --project /srv/durum-project --scratch /srv/durum-scratch \
  --phenotypes /srv/durum-phenotypes --remote /srv/durum-remote \
  --micromamba /opt/micromamba --compute-host compute-host
```

Copy the configured `gwas_frost/scripts` directory into the project layout and supply the corresponding inputs. Run scripts explicitly in the intended environment; do not launch all shell scripts indiscriminately. Server helpers copy files and launch remote jobs. Reference download helpers access public resources. Historical source-file names and table schemas are retained because downstream scripts depend on them.

## Inputs and analysis stages

1. Phenotypes and identifier/marker tables: `01_id_map_pheno_markers.py`, then `01_blue_h2.R`. Prepared inputs include `01_pheno_long.csv`, `01_id_map.csv`, and marker tables under `gwas_frost/tables/`. The source data layout is documented in the preparation script.
2. Genotypes: `02_make_sample_lists.py` and `02_qc_impute.sh stream1`. The QC script expects the project VCF `FINAL-181124.vcf.gz`, sample lists and chromosome renaming table; it uses bcftools, PLINK and Beagle.
3. Structure and associations: `03_structure.sh`, `03_structure.R`, `04_gwas.R`, `09_gwas.R` and their run helpers. The scripts expect prepared phenotype, covariate, genotype and kinship files. Consult each script for its supported arguments and output naming.
4. Interpretation: the `05_*`, `09_*`, `10_*` and `14_*` scripts consume association results, LD/marker information and gene annotations. Sequence comparisons and primer/protein analyses are staged in `16_*` and `17_*`; their scripts specify the reference and intermediate files needed at each stage.

Full statistical and biological methods are provided in the manuscript and Supplementary File S1. Independent execution requires the requested input data and external reference resources; this public archive does not contain enough data to rerun the GWAS by itself.

## Reproduce the final figures

Install the plotting dependencies and run:

```bash
python -m pip install -r figures/requirements.txt
python figures/make_figures.py --data-dir /srv/association-results --output-dir /srv/manuscript-figures
```

The input directory must contain the following uncompressed CSV files:

| File | Required columns | Use |
| --- | --- | --- |
| `g09_S1_FTall_A.csv` | `chr`, `pos`, `p` | Origin-adjusted FT_all association results |
| `g09_RU_FTm12_A.csv` | `chr`, `pos`, `p` | Origin-adjusted FT_m12 association results |
| `g04_s1_GEMMA_A_ang.csv` | `trait`, `chr`, `pos`, `p` | Regime-specific associations without origin adjustment |

Chromosome labels are `1A`–`7B`, positions are base pairs, and `p` contains positive p-values. The third file uses `FT_all_ang`, `FT_m12_ang`, `FT_m15_ang`, `FT_m16_ang`, and `FT_m17_ang`. Provide the complete association results, including the chromosomes and intervals plotted by the manuscript. Outputs are `Figure1.jpg`, `Figure2.jpg` and `Figure2_values.csv`. This command plots existing associations; it does not estimate association p-values.

## Citation and authors

Use `CITATION.cff` for the software citation. The release identifies the manuscript authors in their approved order: Mikhail Divashuk, Aleksey Yanovsky, Varvara Gorskaia, Ludmila Bespalova, Aleksey Ermolaev, Kristina Petrova and Gennady Karlov. Aleksey Ermolaev contributed bioinformatics software. Detailed author contributions are provided in the manuscript.

The research was supported by the Russian Science Foundation, grant 24-16-00274.

## License

The code is released under the MIT License; see `LICENSE`. This license covers the released code and documentation. It does not grant access to, or license, the underlying research data or third-party reference resources and model weights.

## Release validation

Python syntax (49 files) and shell syntax (38 files) were checked. The figure-generation utility was run against the private local association results: both JPEGs were produced at 600 dpi, their layout was compared with the submission figures, and the numerical matrix for Figure 2 matched. A separately configured copy passed Python syntax checks. The package was scanned for credentials and original machine-specific paths. R analyses, GWAS, imputation and protein models were not rerun during release preparation. The checks do not establish full end-to-end reproducibility without the input data and original software environment.
