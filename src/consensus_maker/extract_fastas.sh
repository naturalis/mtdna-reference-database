#!/usr/bin/env bash
SAMPLE="$1"
DEPTH=5
DATA=/Users/rutger/Dropbox/documents/projects/dropbox-projects/parchment/data
SNPS=${DATA}/snps

/usr/local/bin/python3.9 extract_fastas.py \
  --snp-vcf ${DATA}/ChrMT-Run8-TAUIND-public.vcf.gz \
  --bed ${SNPS}/${SAMPLE}_regionfile.bedgraph \
  --parchment-vcf ${SNPS}/${SAMPLE}_vcf.vcf \
  --breeds ../../doc/Breed_Country_Region.tsv \
  --ref ${DATA}/ref.fa \
  --depth ${DEPTH} > ${SAMPLE}.fa
