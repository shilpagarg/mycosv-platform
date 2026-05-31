**Table S2.** Pangenome-stress simulation with donor-derived off-reference sequence. This displayed subset includes 12 ecological scenarios, 3 queries/scenario, 5 contigs/query, 100 kbp total/query, and 1% background divergence as Table S1, but sets `--off-ref-contigs 3 --pangenome-stress` so that 60% of each query's embedded events are sequence absent from the conspecific single reference but represented in a non-conspecific donor reference elsewhere in the pangenome. MycoSV is scored on the same two axes as Table S1 after excluding `NOVEL_WEAK` calls. Stress-mode `svim-asm` and PGGB comparator columns report single-reference F1; PGGB decomposed VCF records were scored as SVs when the inferred REF/ALT allele-length difference was >=30 bp.

| Scenario | Truth | SR P | SR R | **SR F1** | PG P | PG R | **PG F1** | svim-asm F1 | PGGB F1 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| core | 15 | 1.000 | 0.400 | **0.571** | 1.000 | 0.733 | **0.846** | 0.571 | 0.571 |
| hgt_receiver | 15 | 0.750 | 0.400 | **0.522** | 0.882 | 1.000 | **0.938** | 0.000 | 0.000 |
| cross_phylum_hgt_stress | 15 | 1.000 | 0.400 | **0.571** | 1.000 | 0.667 | **0.800** | 0.000 | 0.000 |
| giant_amf | 15 | 0.857 | 0.400 | **0.545** | 0.917 | 0.733 | **0.815** | 0.400 | 0.276 |
| ectomycorrhizal | 15 | 1.000 | 0.400 | **0.571** | 1.000 | 0.800 | **0.889** | 0.222 | 0.207 |
| soil_endophyte | 15 | 1.000 | 0.400 | **0.571** | 1.000 | 0.667 | **0.800** | 0.500 | 0.414 |
| necrotrophic | 15 | 0.857 | 0.400 | **0.545** | 0.900 | 0.600 | **0.720** | 0.571 | 0.387 |
| pathogenic | 15 | 1.000 | 0.400 | **0.571** | 1.000 | 0.867 | **0.929** | 0.500 | 0.148 |
| rust_smut_te_heavy | 15 | 1.000 | 0.400 | **0.571** | 1.000 | 0.733 | **0.846** | 0.500 | 0.188 |
| two_speed_pathogen_extreme | 15 | 1.000 | 0.400 | **0.571** | 1.000 | 0.733 | **0.846** | 0.333 | 0.071 |
| smut | 15 | 0.857 | 0.400 | **0.545** | 0.917 | 0.733 | **0.815** | 0.500 | 0.571 |
| strict_endophyte | 15 | 0.429 | 0.200 | **0.273** | 0.765 | 0.867 | **0.812** | 0.571 | 0.571 |
