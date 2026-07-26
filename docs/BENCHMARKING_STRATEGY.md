# Real-data SV benchmarking for fungi

This benchmark separates two different questions:

1. **Comparator comparison**: how MycoSV behaves relative to existing SV
   tools on the same fungal inputs.
2. **Independent validation**: whether MycoSV calls are supported by raw
   FASTQ/read evidence from the same sample.

For fungi, comparator outputs are not treated as ground truth. Generic SV
callers can miss or misrepresent repeat-rich, accessory-chromosome,
heterokaryotic/dikaryotic, TE-rich, AMF-scale, and highly fragmented fungal
genomes. They are useful baselines, but their shared blind spots should not
define biological truth.

---

## 1. Comparator Outputs Are Baselines

Tools such as minigraph, SVIM-asm, AnchorWave, pggb, cactus, Sniffles, cuteSV,
Delly, and Manta are reported as comparator baselines. They answer:

- Does MycoSV recover calls also found by other methods?
- Where does MycoSV disagree with each comparator?
- Which SV types or fungal architectures drive disagreement?

`exact_benchmark_summary.tsv` keeps the historical `truth_label` column for
compatibility, but now also writes `validation_basis`:

| `validation_basis` | Meaning |
|--------------------|---------|
| `raw_read_validated` | Tool-agnostic candidate union validated directly against raw reads/FASTQ. This is the preferred independent validation basis. |
| `comparator_agreement_read_supported` | Comparator consensus or comparator candidate set after a raw-read support filter. Useful, but still comparator-seeded. |
| `comparator_agreement` | Calls supported by multiple comparator tools. Use as baseline agreement, not ground truth. |
| `comparator_baseline` | MycoSV scored against one comparator output. Diagnostic only. |
| `no_independent_validation` | MycoSV ran, but no comparator/read validation basis was available. Metrics are undefined. |

---

## 2. Primary Fungal Validation

The primary real-data validation target is `truth_label=read_level_union` with
`validation_basis=raw_read_validated`.

This path pools candidate loci from MycoSV and available comparators, then
keeps only loci supported by raw reads. When `force_external=True`, MycoSV's
own internal support cannot self-validate a call. This makes the row a direct
read-evidence question:

> Given all proposed loci, which ones are independently supported by raw
> FASTQ/read alignments, and how well does each method recover them?

For long reads, validation uses split/clipped alignment evidence around
breakpoints. For short reads, it uses paired/split/clipped/depth-style support
where available. For assembly-mode queries, matching raw reads are still the
best validation source; without them, assembly-only evidence is reported but
should not be presented as independently validated biology.

---

## 3. Comparator Agreement

Comparator consensus rows such as `consensus_2of_N` remain useful, but only as
baseline agreement. They are not the headline biological truth for fungi.

Use them to answer:

- Does MycoSV agree with multiple existing tools?
- Are disagreements concentrated in INS/DEL/DUP/INV/TRA/OFF_REF?
- Is a result sensitive to one comparator? Check `loo_consensus_summary.tsv`
  and `loo_consensus_variance.json`.

Do not describe `consensus_2of_N` as ground truth in reports. Prefer
phrases such as “comparator agreement,” “baseline agreement,” or
“multi-tool support.”

---

## 4. Reporting Rules

Recommended ordering for real fungal reports:

1. **Independent validation**: use `read_level_union` rows when present.
2. **Read-filtered comparator agreement**: use `*_read_supported` consensus
   rows when `read_level_union` is unavailable.
3. **Comparator agreement**: use `consensus_2of_N` as baseline comparison
   only.
4. **Single-comparator rows**: show as diagnostics, not as truth.
5. **MycoSV-only runs**: report call burden, evidence tiers, TE/biology
   annotations, and raw-read support if available; do not report P/R/F1 as
   accuracy metrics when validation is absent.

Biology tables should feature calls with independent read support first:

- `strong`: comparator-supported and raw-read supported
- `moderate`: comparator-supported or raw-read supported
- `intrinsic_only`: supported only by MycoSV internal evidence
- `weak`: exploratory

For fungal novelty, `mycosv_only_read_supported` is a meaningful category:
the call is not found by current comparators but is supported by raw reads.

---

## 5. Remaining Gaps

- Add explicit matching of assembly queries to same-sample FASTQ/SRA accessions
  wherever possible.
- Add HiFi reassembly-derived validation when suitable reads exist.
- Keep synthetic spike-in benchmarks as calibration, because they provide exact
  truth by construction but do not replace real fungal raw-read validation.
- Make every visualization label use “validation basis” or “baseline
  agreement” instead of “ground truth” for real-data comparator rows.
