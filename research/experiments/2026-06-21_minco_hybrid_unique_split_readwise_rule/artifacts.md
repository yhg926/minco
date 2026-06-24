# Artifacts

Record paths, sizes, and whether files are temporary or should be preserved.

| Kind | Path | Description | Preserve? |
| --- | --- | --- | --- |
| input | `/mnt/new3T/minco_cami3_toygut_20260620/sample_0_anonymous_reads.fq.gz` | CAMI3 toy human-gut sample0 reads | yes |
| input | `/mnt/new3T/minco_cami3_toygut_20260620/taxonomic_profiles/taxonomic_profile_0.txt` | CAMI3 toy sample0 gold profile | yes |
| input | `/tmp/cami_marine_sample0_reads.fq.gz` | CAMI2 marine sample0 reads cache | maybe |
| input | `/tmp/gs_marine_short.profile` | CAMI2 marine gold profile cache | maybe |
| input | `/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno` | GTDB-only S1000 minco reference sketch | yes |
| input | `/tmp/gtdb232_refseqvirus_s1000.minco.cami_taxmap.full_lineage.tsv` | Species taxid/lineage map used for scoring | maybe |
| input | `/mnt/new3T/minco_cami3_toygut_20260620/minco_s1000_gtdbonly_sample0_unfiltered.tsv` | Existing toy unique ZIP unfiltered candidates | yes |
| input | `/tmp/minco_toy_s1000_gtdbonly_split_zip_unfiltered.tsv` | Existing toy split ZIP unfiltered candidates | maybe |
| input | `/mnt/new3T/minco_cami3_toygut_20260620/sylph_sample0/profile.tsv` | Toy Sylph profile | yes |
| input | `/tmp/sylph_marine_sample0/profile.tsv` | Marine sample0 Sylph profile | maybe |
| output | `/tmp/minco_hybrid_20260621/marine_s1000_gtdb_split_zip_sample0_unfiltered.tsv` | Newly generated marine split ZIP unfiltered candidates, 2:11 wall, 4.12 GB RSS | maybe |
| output | `/tmp/minco_hybrid_20260621/marine_s1000_gtdb_unique_zip_sample0_unfiltered.tsv` | Newly generated marine unique ZIP unfiltered candidates, 1:53 wall, 3.68 GB RSS | maybe |
| output | `/tmp/minco_hybrid_20260621/toy_gtdb_bacteria.joined_features.tsv` | Toy joined best-row species features | maybe |
| output | `/tmp/minco_hybrid_20260621/marine_gtdb_species.joined_features.tsv` | Marine joined best-row species features | maybe |
| output | `/tmp/minco_hybrid_20260621/hybrid_grid_all.tsv` | Full hybrid threshold grid results | maybe |
| output | `/tmp/minco_hybrid_20260621/toy_bacteria_gold_taxids.tsv` | Toy bacterial gold species taxids used by local scorer | maybe |
| output | `/tmp/minco_hybrid_20260621/marine_gold_taxids.tsv` | Marine gold species taxids used by local scorer | maybe |
| output | `/tmp/minco_hybrid_20260621/classifier_leave_one_dataset_out.tsv` | Classifier generalization output before repo copy | maybe |
| output | `research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/summary.tsv` | Main comparison table copied into repo | yes |
| output | `research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/hybrid_best_params.tsv` | Joint-best hybrid parameters | yes |
| output | `research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/hybrid_grid_best50.tsv` | Best 50 joint hybrid parameter rows | yes |
| output | `research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/classifier_leave_one_dataset_out.tsv` | Leave-one-dataset-out classifier comparison copied into repo | yes |
| output | `research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/toy_bacteria_gold_taxids.tsv` | Toy bacterial gold species taxids copied into repo | yes |
| output | `research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/marine_gold_taxids.tsv` | Marine gold species taxids copied into repo | yes |
| script | `research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/scripts/hybrid_unique_split_score.py` | Reusable hybrid scoring/grid script | yes |
| script | `research/experiments/2026-06-21_minco_hybrid_unique_split_readwise_rule/scripts/train_hybrid_classifier.py` | Leave-one-dataset-out classifier test script | yes |
