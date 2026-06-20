# CAMI II Clinical Pathogen Detection Minco Benchmark

Date: 2026-06-20

## Question

Can minco analyze the CAMI II Clinical pathogen detection challenge metagenome, and what signal does it report against the current S1000 viral and appended GTDB+viral reference sketches?

## Dataset

CAMI dataset: CAMI II Clinical pathogen detection, `patmgCAMI2.tar.gz`.

Remote URL: `https://frl.publisso.de/data/frl:6425521/patmgCAMI2.tar.gz`

Local tarball: `/mnt/new3T/gtdbr220/test/patmgCAMI2.tar.gz`

Extracted reads:

- `/tmp/patmg_CAMI2_reads_20260620/patmg_CAMI2/patmg_CAMI2_short_read_R1.fastq.gz`
- `/tmp/patmg_CAMI2_reads_20260620/patmg_CAMI2/patmg_CAMI2_short_read_R2.fastq.gz`

Metadata reports one real blood metagenome, MiSeq, 1x100 bp, total size 0.69 Gbp. The paired FASTQ files each contain 13,778,280 lines, or 3,444,570 reads per mate; combined minco input was 6,889,140 reads.

## References

Appended S1000 DB:

`/mnt/new3T/gtdbr220/gtdb232/GTDBr232_minco_T_S1000_anno_20260619/sketch_T_S1000_anno_plus_refseq_virus_20260620`

Viral S1000 DB:

`/mnt/new3T/IMG_VR_data/refseq_viro/refseq_virus_minco_S1000_anno_20260619`

## Commands

Commands are recorded in `commands.sh`.

Main parameters:

`--qraw - --query-density ref --abundance-est depth -p16 --density-block-ctx 100`

Filtered scan:

`-f0.01 -n0.9 -m0`

Permissive viral scan:

`-f0 -n0 -t1 -m0`

## Results

Filtered appended DB scan:

- Runtime: 1:32.82
- CPU: 402%
- Max RSS: 10.54 GB
- Output rows: 81 plus header
- Total reads: 6,889,140
- Unique query contexts: 261,493,794
- Density blocks: 3,444,564

Filtered viral-only scan:

- Runtime: 1:16.14
- CPU: 463%
- Max RSS: 7.51 GB
- Output rows: 44 plus header

The strongest filtered signal was Lambda-like phage:

- `NC_049953.1 Escherichia phage Lambda_ev099`: ANI 0.996799, XnY 769,694, RefAF 0.015, Reads 1,539,356, normalized abundance 0.4986.
- `NC_049949.1 Escherichia phage Lambda_ev207`: ANI 0.996783, XnY 769,827, RefAF 0.016, Reads 1,539,466, normalized abundance 0.4987.

The permissive viral-only scan found a clinical pathogen candidate that was filtered out by `-n0.9`:

- `NC_005301.3 Crimean-Congo hemorrhagic fever virus segment L`: ANI 0.892146, XnY 90, RefAF 0.041, Reads 122, unique ref hit 41, breadth 0.041, normalized abundance 0.000007.
- `NC_078225.1 Crimean-Congo hemorrhagic fever virus 2 strain AP92 segment L`: ANI 0.852468, XnY 35, RefAF 0.015, Reads 48.

## Interpretation

For this CAMI clinical pathogen dataset, a strict `ANI >= 0.9` filter emphasizes Lambda-like phage and misses the CCHF candidate. A permissive viral scan recovers Crimean-Congo hemorrhagic fever virus segment L, but with low read support and low abundance. This suggests minco can surface the expected pathogen only if pathogen-detection mode uses a lower ANI threshold or a special low-abundance viral reporting rule.

## Caveats

- This is one CAMI sample only.
- The combined appended DB uses universal density 1.0 due to small viral references, increasing query context load.
- Current readwise default calls mark all rows as `weak`; pathogen-detection reporting needs different thresholds from abundance/major-call reporting.
- The CAMI challenge truth was not encoded in the downloaded metadata, only the challenge description link was present.

