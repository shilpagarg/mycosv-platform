# MycoSV manuscript — key take-home messages

Headline framing for the cover-letter, abstract, and Discussion. Each message
is anchored in a specific figure or table in this folder so a reviewer (or
co-author drafting prose) can verify the number in one click.

---

## TL;DR

> MycoSV is the first pangenome-graph-routed structural variant caller for
> fungi. On a panel of phylogenetically diverse fungal genomes it
> (i) recovers **15× more SV loci per sample** than the leading
> single-reference comparators, (ii) supports those calls at a **0.81
> independent-read-validation rate** (vs 0.51–0.67 for comparators on the
> same evidence), (iii) uniquely captures **off-reference novel sequence
> and cross-clade lineage labels** that single-reference approaches cannot
> represent by construction, and (iv) makes **76% of its loci pangenome-only
> (not recoverable on any held-out single reference)**.

---

## 1. Pangenome-only loci are the dominant signal

> Across the 5-sample completed panel, **29,400 of 38,588 MycoSV loci (76.2%)
> are pangenome-only** — they cannot be represented in single-reference
> coordinate space. Single-reference comparators (AnchorWave, svim-asm) do
> not target this category by construction; minigraph captures only the
> subset traversable in its bubble graph.

- **Source:** [tables/table_method_comparison_pangenome_vs_single_ref.md](tables/table_method_comparison_pangenome_vs_single_ref.md), "Pangenome-only loci" column
- **Visual:** [figures/fig1b_panel_pangenome_lift.png](figures/fig1b_panel_pangenome_lift.png) (per-sample stacked bars; green = pangenome-only read-supported)

## 2. MycoSV is not noisy — it has a higher quality bar than comparators

> Per-sample read-validation rate (fraction of each caller's truth set
> independently supported by raw-read split alignments in the same query):
> **MycoSV 0.81, minigraph 0.67, AnchorWave 0.53, svim-asm 0.51.** MycoSV
> trades neither quantity nor quality — both win.

- **Source:** [tables/table_method_comparison_pangenome_vs_single_ref.md](tables/table_method_comparison_pangenome_vs_single_ref.md), "Read-validation rate" column
- **Visual:** [figures/fig1a_panel_read_validation_rate.png](figures/fig1a_panel_read_validation_rate.png) (per-sample grouped bars, Wilson 95% CIs)
- **Anchor:** the validation BAMs are built from raw FASTQs independently of any caller; comparator and MycoSV truth sets are scored against the same per-sample alignment.

## 3. Off-reference novel sequence is a category only MycoSV represents

> MycoSV reports **5,194 OFF_REF insertions** across the panel — insertions
> whose alternate allele has no homologous reference anchor. Single-reference
> callers (AnchorWave, svim-asm) and minigraph's bubble extraction emit
> zero such calls by design. This category is where new fungal accessory
> genome and HGT cargo lives.

- **Source:** [tables/table_method_comparison_pangenome_vs_single_ref.md](tables/table_method_comparison_pangenome_vs_single_ref.md), "OFF_REF novel sequence" column
- **Detail:** [tables/table1_pangenome_vs_single_reference.md](tables/table1_pangenome_vs_single_reference.md), per-sample HGT/Starship and TE/RIP columns

## 4. Even on single-reference scale, MycoSV wins

> Even if you restrict to MycoSV's single-reference-equivalent subset
> (9,188 calls across the panel), that's still **~2× more than the best
> comparator** (svim-asm 4,031; minigraph 4,174; AnchorWave 4,918). The
> pangenome advantage is additive, not a tradeoff: MycoSV doesn't sacrifice
> single-reference recall to gain pangenome-only loci.

- **Source:** [tables/table_method_comparison_pangenome_vs_single_ref.md](tables/table_method_comparison_pangenome_vs_single_ref.md), "Single-reference calls" column

## 5. Cross-clade lineage and element-class labels are unique to MycoSV

> Every MycoSV call carries a routed-clade rank tag (`CLADE` / `CLADE_RANK`
> at the species / family / class / phylum level) and an element-class tag
> (`EC`: HGT / Starship / RIP / TE_LTR / TE_TIR / TE_LINE / TE_SINE / REPEAT
> / NONE). None of the three comparators provide either. These labels are
> what power Figure 2 cross-guild biology stratification and Table 1 per-sample
> biology breakdowns. RIP labels are gated to Pezizomycotina classes per
> standard fungal genome-defense biology.

- **Source:** [tables/table_method_comparison_pangenome_vs_single_ref.md](tables/table_method_comparison_pangenome_vs_single_ref.md), "Cross-clade lineage / element-class labels" column

---

## Figure 1 (this folder) — manuscript validation panels

| Panel | Claim | File |
|---|---|---|
| 1A | MycoSV's per-sample read-validation rate exceeds every comparator on the same independent raw-read evidence | [figures/fig1a_panel_read_validation_rate.svg](figures/fig1a_panel_read_validation_rate.svg) (vector); .png also present |
| 1B | 60–95% of MycoSV's per-sample loci are pangenome-only and 100% of those are read-supported on this panel | [figures/fig1b_panel_pangenome_lift.svg](figures/fig1b_panel_pangenome_lift.svg) (vector); .png also present |

## Figure 2 (pending 165-sample run) — biology at scale

Both panels build automatically once the 165-sample shards finish. Panel A
re-uses the cross-guild enrichment volcano on the 4-way guild split (AMF /
Filamentous / Basidio / Yeast); Panel B compares circos SV landscapes side-
by-side for representative AMF + Filamentous genomes.

## Tables

| Table | Granularity | Audience | File |
|---|---|---|---|
| Headline method-comparison | Method-level (rows = callers, columns = aggregate panel metrics) | Reviewer's first read — "why pangenome?" | [tables/table_method_comparison_pangenome_vs_single_ref.md](tables/table_method_comparison_pangenome_vs_single_ref.md) |
| Per-sample detail | Sample-level (rows = each genome, columns = recovery + validation + biology breakdown) | Supplementary / methods reader who wants every sample visible | [tables/table1_pangenome_vs_single_reference.md](tables/table1_pangenome_vs_single_reference.md) |

## Data files

Raw panel-fold aggregates underlying the figures/tables, for downstream
reuse without re-running the pipeline:

- [data/panel_read_validation_rate.tsv](data/panel_read_validation_rate.tsv) — per-(query, source) yes / total / rate + Wilson 95% CI
- [data/panel_pangenome_lift.tsv](data/panel_pangenome_lift.tsv) — per-query raw / dedup / single-ref / pangenome-only / read-supported counts

---

## Numbers a reviewer will quote

| Claim | Number | Source |
|---|---|---|
| Pangenome lift (panel) | **76.2%** | Method-comparison table, MycoSV row, "Pangenome-only loci" |
| MycoSV read-validation rate | **0.814** | Method-comparison table, MycoSV row |
| Best comparator read-validation rate | 0.670 (minigraph) | Method-comparison table, minigraph row |
| MycoSV median calls per sample | **6,264** | Method-comparison table, MycoSV row |
| Best comparator median calls per sample | 429 (AnchorWave) | Method-comparison table, AnchorWave row |
| MycoSV total panel loci | **38,588** | Method-comparison table, MycoSV row |
| MycoSV OFF_REF novel sequence | **5,194** | Method-comparison table, MycoSV row |
| Comparator OFF_REF novel sequence | 0 (not supported) | Method-comparison table, all comparator rows |

All numbers are reproducible by re-running:

```
.venv/bin/python plot_mycosv_pangenome_calls.py \
  --panel-dir experiments/million_real/full_fungal_assembly_20260522_070100/assembly/by_query \
  --outdir manuscript/_regen
```

Folder structure is identical to this one; diff the regenerated outputs to
verify reproducibility.

---

## Caveats to flag in the Discussion

- **Panel completeness:** the headline numbers come from 5 of 14 shards that
  completed the full MycoSV pipeline. The other 9 timed out and have a
  `MYCOSV_FAILED.txt` marker. Both tables disclose this in their captions so
  the panel size is never silently misrepresented.
- **165-sample scale-up:** the 165-sample panel is staged (prepared/) but
  not yet computed. Figure 2 (cross-guild biology) and an updated headline
  table on the larger panel will be regenerated automatically once those
  shards finish.
- **Comparator-truth quality flag:** the per-query
  `exact_benchmark_summary.tsv` carries a `comparator_truth_quality` column
  (`concordant` / `discordant_overcalled` / `discordant_undercalled`) that
  flags cases where comparators disagree by >5× (e.g., minigraph
  over-segmenting bubbles on the Trichoderma genome). Aggregate numbers in
  the headline table treat all comparator truth equally; per-sample numbers
  in Table 1 allow stratification.
- **Single-reference F1:** intentionally not in the headline table. Pangenome
  calls scored against single-reference truth in single-reference coordinates
  carry a known contig-mismatch artifact (54% of failure rows in
  `match_failures.tsv` per query) that we addressed with a query-projected
  scoring fallback in the pipeline. The read-validation rate column is the
  yardstick this table presents.
