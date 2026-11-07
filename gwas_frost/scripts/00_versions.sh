#!/usr/bin/env bash
# Версии программ окружения gwas_frost -> gwas_frost/software_versions.txt
source "$(dirname "$0")/00_env.sh"
{
  echo "# Версии программ (сгенерировано $(date -Iseconds))"
  echo "OS: $(lsb_release -ds 2>/dev/null || uname -a)"
  echo "bcftools: $(bcftools --version | head -1)"
  echo "htslib/tabix: $(tabix --version | head -1)"
  echo "plink2: $(plink2 --version)"
  echo "plink1.9: $(plink --version)"
  echo "beagle: $(ls $ENVBIN/../share | grep -i beagle | head -1) ; $(beagle 2>&1 | head -1)"
  echo "gemma: $(gemma -license 2>&1 | grep -m1 GEMMA)"
  echo "java: $(java -version 2>&1 | head -1)"
  echo "R: $(R --version | head -1)"
  echo "vcftools: $(vcftools --version 2>&1 | head -1)"
} > "$OUT/software_versions.txt" 2>&1
cat "$OUT/software_versions.txt"
{
  echo "GAPIT: $(grep -m1 'GAPIT.Version=' /work/scratch/soft/gapit_functions.txt) (source() из https://zzlab.net/GAPIT/gapit_functions.txt, загружен 2026-10-25; пакет GitHub не загружается из-за зашитого пути автора)"
  echo "minimap2: $(minimap2 --version 2>/dev/null || echo 'не установлен')"
} >> "$OUT/software_versions.txt"
