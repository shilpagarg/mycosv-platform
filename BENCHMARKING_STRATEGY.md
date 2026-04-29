# Real-data SV benchmarking: ground-truth strategy

This document explains why "minigraph as truth" gives misleading numbers,
how the current pipeline avoids that with multi-comparator consensus, and
which bias-free alternatives are available for raw-data validation.

---

## 1 — Why minigraph alone is not a clean ground truth

`minigraph -cxggs` builds a pangenome graph and `gfatools bubble` extracts
"bubbles" (regions of structural divergence). On the
[20260428_182921 compact_yeast assembly](experiments/real_data/20260428_182921/compact_yeast/benchmark_assembly/exact_benchmark_summary.tsv)
benchmark this caller emits **8 SVs** on `GCA_030569935_1` while
cactus/svim_asm/anchorwave emit 110–170. The reasons are well-known:

1. **Bubble extraction is conservative.** `gfatools bubble` only emits
   isolated bubbles — long divergent stretches collapse into a single
   bubble or are skipped entirely. Truth callers that work on pairwise
   alignments (`svim_asm`, `anchorwave`) recover those events as
   individual INS/DEL/INV.
2. **Default `minigraph -xggs` uses asm5 chains** (≤5% divergence). For
   intra-genus comparisons (e.g. *Nakaseomyces glabratus* haplotypes,
   *Saccharomyces* hybrids) regions exceeding that drop out of the graph
   entirely.
3. **No nested-event recovery.** If a 50 kb inversion contains a 3 kb
   deletion, minigraph reports one bubble; a per-base aligner reports
   both.
4. **Bias against small SVs near assembly contig ends.** The bubble
   extraction trims terminal bubbles; svim_asm doesn't.

The empirical signal: across all four assembly truth labels in the
20260428 run, minigraph consistently emits 5–20 % of the SV count of the
other three callers. Treating its output as the single ground truth makes
mycosv (or any tool) look like it has 0 % recall against minigraph even
when 7–10 % of mycosv's predictions match the consensus of the other
three callers. Use minigraph as **one signal**, not as the truth.

---

## 2 — What the pipeline now reports (after the fixes)

[exact_benchmark_summary.tsv](experiments/real_data/20260428_182921/compact_yeast/benchmark_assembly/exact_benchmark_summary.tsv) now carries three views per `(query, mode)`:

| `coordinate_space`        | what it scores |
|----------------------------|----------------|
| `reference`                | mycosv ref-coord calls *filtered* to `benchmark_ref_fasta` contigs vs each truth caller |
| `reference_any_clade`      | **same truth, but mycosv predictions NOT filtered** — exposes how many mycosv calls were correctly anchored but on a sibling clade |
| `reference` w/ `truth_label=consensus_2of_N` | mycosv vs the high-confidence consensus (an SV is in the consensus iff it's compatible across position + length + type tolerance with calls from at least 2 of the N comparators) |

The consensus row is the **primary number**; the per-tool rows let you
see which comparator drives the consensus and how badly any single tool
disagrees. `benchmark_summary.json` additionally reports
`mycosv_calls.off_ref_dropped` (calls with no REFPOS, un-matchable in
ref space) and `mycosv_calls.misrouted_to_sibling_clade` (calls whose
REFCONTIG did not pass the `benchmark_ref_fasta` filter), so the
prediction count is no longer opaque.

Implementation: see
[run_real_fungal_benchmark.py:build_consensus_truth](run_real_fungal_benchmark.py)
and the per-query metric loop just below it.

---

## 3 — Truly bias-free truth, ranked by how much it costs

For each option below, "raw data" means going back to the FASTQ and not
trusting any single caller's annotation. They compose: you can apply 3
on top of 1+2.

### 3.1 — Multi-comparator consensus *(implemented, free, recommended primary)*

Truth = SVs supported by ≥ 2 of {cactus, svim_asm, anchorwave, pggb}.

- Pros: zero new infrastructure, removes single-tool bias, runs on the
  comparator outputs already produced by `benchmark_<mode>/comparators/`.
- Cons: still bounded by what **all** assembly comparators can find
  (e.g. complex recombination is undercalled by every tool).
- Knob: `min_support` in `build_consensus_truth()`. 2-of-4 is sensitive,
  3-of-4 is highest-confidence.
- **What the existing TSV now uses.** Treat `truth_label=consensus_2of_4`
  as the headline metric; treat per-tool rows as diagnostic.

### 3.2 — Read-level support filter *(implemented as a post-filter)*

For every SV in any truth set, require ≥ K split-read alignments from raw
long reads of the same sample. SVs with no read support are removed from
both truth and predictions before scoring. This anchors the metric in the
raw FASTQ instead of in any caller.

- Pros: filters caller artefacts (paftools INS at MNV sites, cactus
  duplications from heterozygosity collapses) without re-implementing
  variant calling.
- Cons: needs raw long reads for the same sample as the assembly; for
  NCBI assembly queries the matching SRA run usually exists but not
  always.
- See [validate_truth_with_reads.py](validate_truth_with_reads.py)
  (added by these fixes).

### 3.3 — De-novo HiFi assembly as truth *(workflow, not yet wired in)*

If raw HiFi reads exist for the sample:

```
hifiasm -o asm -t 32 reads.hifi.fa.gz
# asm.bp.p_ctg.fa is the diploid primary assembly, ≥ Q40 per base
svim-asm haploid out_dir asm.bp.p_ctg.fa ref.fa > truth.vcf
```

That `truth.vcf` is the cleanest single-sample truth available because
the assembly itself was reconstructed from the raw reads, and the SV
caller only sees the resulting consensus contigs (no graph collapse).
For fungi this typically gives Q50+ contigs and resolves nested events.
Plug the resulting VCF into the pipeline via `--other-vcf
hifi_truth=truth.vcf` and the per-query loop already picks it up.

### 3.4 — Synthetic spike-in *(implemented as standalone benchmark)*

For a given reference, programmatically inject N known SVs (INS/DEL/DUP/
INV at random positions of random length sampled from a realistic
distribution), regenerate the "query" FASTA, and run the full pipeline.
Truth is exact by construction.

- Pros: 100 % perfect ground truth; tests positional accuracy directly.
- Cons: misses real-data effects (heterozygosity, transposons,
  centromeric satellite regions) — pair with §3.1–3.3 for a complete
  picture.
- See [synthetic_sv_benchmark.py](synthetic_sv_benchmark.py)
  (added by these fixes).

### 3.5 — Trio / family validation *(workflow, not implemented)*

For trios (parent×parent→F1) call SVs in all three samples. Mendelian
inconsistency is a strong "false positive" signal. For fungi this works
on dikaryotic crosses and on lab progeny (e.g. *Saccharomyces*
crosses, *Zymoseptoria tritici* mating populations). Not in scope here —
flagged for the read-mode work.

---

## 4 — Recommended report layout

When generating the visualization report and biological findings:

1. **Headline P/R/F1**: use the `consensus_2of_N` row from
   `exact_benchmark_summary.tsv`. Drop minigraph from the headline.
2. **Per-tool rows**: keep, but render as a small-multiples panel under
   the headline so the user sees disagreement.
3. **Routing diagnostics**: render `mycosv_calls.off_ref_dropped`
   (un-matchable novel-sequence events) and
   `mycosv_calls.misrouted_to_sibling_clade` so a 70 / 111 / 233 split
   is visible at a glance.
4. **Synthetic-spike-in scorecard**: include as a "calibration" panel —
   the pipeline ought to clear ≥ 90 % recall on synthetic data; if it
   doesn't, real-data numbers below that ceiling reflect real-data
   difficulty rather than caller bias.
5. **Biological findings**: only include SV candidates that pass
   *both* the consensus filter and the read-level support filter
   (§3.2). Anything below `consensus_2of_2` and read-support ≥ 3 should
   be marked as "low-confidence" rather than featured.

---

## 5 — Concrete next steps

- [x] Multi-comparator consensus, any-clade row, OFF_REF /
      misrouted diagnostics (this commit).
- [x] [retry_real_panels.sh](retry_real_panels.sh) to rerun under
      `srun --mem=32G` and bypass the 12 GiB cap.
- [x] [validate_truth_with_reads.py](validate_truth_with_reads.py) for
      §3.2.
- [x] [synthetic_sv_benchmark.py](synthetic_sv_benchmark.py) for §3.4.
- [ ] Wire the consensus row into
      [sv_visualization_report.py](sv_visualization_report.py) as the
      headline metric.
- [ ] Add HiFi-assembly truth ingestion when the source FASTQ is
      identified for a query (§3.3).
