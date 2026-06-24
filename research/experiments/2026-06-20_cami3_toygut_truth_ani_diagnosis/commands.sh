#!/usr/bin/env bash
set -euo pipefail

ROOT=/mnt/new3T/minco_cami3_toygut_20260620

# Extract the sample0 read-to-source mapping from the read archive.
tar -xOzf "$ROOT/sample_0_reads.tar.gz" sample_0_reads/reads_mapping.tsv.gz \
  > "$ROOT/sample_0_reads_mapping.tsv.gz"

# Download source-genome metadata and source genomes.
curl --retry 5 --retry-all-errors --retry-delay 5 -fL \
  https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/gsa_pooled_mapping.tsv.gz \
  -o "$ROOT/gsa_pooled_mapping.tsv.gz"

curl --retry 5 --retry-all-errors --retry-delay 5 -C - -fL \
  https://s3.bi.denbi.de/swift/v1/cami/cami3_toydata/human-gut-toy/source_genomes.tar.gz \
  -o "$ROOT/source_genomes.tar.gz"

# Extract only source genomes for the first ANI-discrepant species panel.
mkdir -p "$ROOT/source_targets"
tar -xzf "$ROOT/source_genomes.tar.gz" -C "$ROOT/source_targets" --wildcards \
  'source_genomes/Faecalibacterium*' \
  'source_genomes/Phascolarctobacterium_succinatutens*' \
  'source_genomes/Blautia_obeum*' \
  'source_genomes/Bacteroides_thetaiotaomicron*'

# Source-ID counts for the diagnostic taxa; read mapping has one row per read mate.
gzip -dc "$ROOT/sample_0_reads_mapping.tsv.gz" |
  awk -F'\t' 'NR>1 {cnt[$2]++; tax[$2]=$3} END {for (g in cnt) if(tax[g]==853 || tax[g]==626940 || tax[g]==40520 || tax[g]==818) print tax[g],g,cnt[g]}' |
  sort -k1,1n -k3,3nr

# Match ASV source IDs to extracted source FASTAs by contig accession.
rg -n 'NZ_QVFB01000006\.1|NZ_CYXN01000028\.1|NZ_JAQDKL010000046\.1|NZ_CP097284\.1|NZ_QSHL01000003\.1|NZ_JBDQBI010000015\.1|NZ_JAQNVQ010000013\.1|NZ_JADMXB010000013\.1' \
  "$ROOT/source_targets/source_genomes"

# fastANI truth checks against the references selected by both minco and Sylph.
fastANI -q "$ROOT/source_targets/source_genomes/Faecalibacterium_prausnitzii_AM37_13AC_strain_genomic.fna" \
  -r /mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_003324185.1_ASM332418v1_genomic.fna.gz \
  -o /tmp/cami3_truth_ASV507.1_vs_GCF003324185.fastani
fastANI -q "$ROOT/source_targets/source_genomes/Faecalibacterium_prausnitzii_2789STDY5834970_genomic.fna" \
  -r /mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_003324185.1_ASM332418v1_genomic.fna.gz \
  -o /tmp/cami3_truth_ASV507.0_vs_GCF003324185.fastani
fastANI -q "$ROOT/source_targets/source_genomes/Faecalibacterium_prausnitzii_AM100_B16A_strain_genomic.fna" \
  -r /mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_003324185.1_ASM332418v1_genomic.fna.gz \
  -o /tmp/cami3_truth_ASV507.3_vs_GCF003324185.fastani
fastANI -q "$ROOT/source_targets/source_genomes/Phascolarctobacterium_succinatutens_NB2A_12_FMU_strain_genomic.fna" \
  -r /mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_000188175.1_ASM18817v1_genomic.fna.gz \
  -o /tmp/cami3_truth_ASV579.0_vs_GCF000188175.fastani
fastANI -q "$ROOT/source_targets/source_genomes/Blautia_obeum_AM37_4AC_strain_genomic.fna" \
  -r /mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_025147765.1_ASM2514776v1_genomic.fna.gz \
  -o /tmp/cami3_truth_ASV406.4_vs_GCF025147765.fastani
fastANI -q "$ROOT/source_targets/source_genomes/Bacteroides_thetaiotaomicron_E8_m1001262Bd2_200113_strain_genomic.fna" \
  -r /mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_000011065.1_ASM1106v1_genomic.fna.gz \
  -o /tmp/cami3_truth_ASV153.8_vs_GCF000011065.fastani
fastANI -q "$ROOT/source_targets/source_genomes/Bacteroides_thetaiotaomicron_1001095st1_D10_1001095H_141210_strain_genomic.fna" \
  -r /mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_000011065.1_ASM1106v1_genomic.fna.gz \
  -o /tmp/cami3_truth_ASV102.10_vs_GCF000011065.fastani
fastANI -q "$ROOT/source_targets/source_genomes/Bacteroides_thetaiotaomicron_D52t1_170925_B7_strain_genomic.fna" \
  -r /mnt/new3T/gtdbr220/GTDBr226_genomes/GCF_000011065.1_ASM1106v1_genomic.fna.gz \
  -o /tmp/cami3_truth_ASV99.3_vs_GCF000011065.fastani

# Feature extraction for the same selected references from minco unfiltered output.
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[$i]=i; print "ref\tminco_ani\tXnY_ctx\tbreadth\tmean_depth\thit_mean_depth\tdepth_cv\tzip_af\treads_with_match\ttotal_reads"} /GCF_003324185\.1/ || /GCF_000188175\.1/ || /GCF_025147765\.1/ || /GCF_000011065\.1/ {print $15"\t"$3"\t"$7"\t"$(h["Ref_breadth"])"\t"$(h["Ref_mean_depth"])"\t"$(h["Ref_hit_mean_depth"])"\t"$(h["Ref_depth_cv"])"\t"$(h["Ref_zip_af"])"\t"$(h["Reads_with_ctx_match"])"\t"$(h["Total_reads"])}' \
  "$ROOT/minco_s1000_gtdbonly_sample0_unfiltered.tsv"
