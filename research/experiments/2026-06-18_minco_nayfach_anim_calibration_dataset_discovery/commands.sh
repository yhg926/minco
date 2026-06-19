#!/usr/bin/env bash
set -euo pipefail

# Search local repo for prior ANIm references.
rg -n -i "nayfach|anim|ani[m_ -]?ground|groundtruth|ground truth" .

# Search mounted data for Nayfach/ANIm artifacts.
find /mnt -maxdepth 5 -iname '*nayfach*' -o -iname '*anim*' -o -iname '*ground*truth*'

# Inspect the primary Nayfach ANIm label table.
head -3 /mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/Nayfach52k.kssd3_codenpatternT10_vs_ANIm.tsv
wc -l /mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/Nayfach52k.kssd3_codenpatternT10_vs_ANIm.tsv
awk 'BEGIN{FS="\t"} NR==1{print "first_row_cols", NF} {n++; truth=$10; raw=$9; xny=$3; if(NR==1||truth<mintruth)mintruth=truth; if(NR==1||truth>maxtruth)maxtruth=truth; if(NR==1||raw<minraw)minraw=raw; if(NR==1||raw>maxraw)maxraw=raw; if(NR==1||xny<minxny)minxny=xny; if(NR==1||xny>maxxny)maxxny=xny} END{print "rows", n; print "truth_min", mintruth; print "truth_max", maxtruth; print "raw_min", minraw; print "raw_max", maxraw; print "xny_min", minxny; print "xny_max", maxxny}' /mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/Nayfach52k.kssd3_codenpatternT10_vs_ANIm.tsv

# Verify alternate T11 table and Nayfach source files.
head -3 /mnt/new3T/skani_data/Nayfach_data/Nayfach52k.kssd3_codenpatternT11_vs_ANIm
wc -l /mnt/new3T/skani_data/Nayfach_data/Nayfach52k.kssd3_codenpatternT11_vs_ANIm
find /mnt/new3T/skani_data/Nayfach_data/fna -maxdepth 1 -type f -name '*.fna' | wc -l
head -3 /mnt/new3T/skani_data/Nayfach_data/genome_metadata.tsv

# Inspect prior model artifacts as references.
sed -n '1,220p' /mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/optimization_summary.md
sed -n '1,220p' /mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/final_optimization_conclusion.md
head -20 /mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/t10_model_parameters_raw_scale.tsv
head -20 /mnt/new3T/gtdbr220/eval_runs/kssd3_ani_model_optimization_20260609/model_parameters_raw_scale.tsv
