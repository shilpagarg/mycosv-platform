# MycoSV-platform v0.1.0 — first tagged release

## Self-contained Starship / HGT annotation
- Element **boundary** resolution + flanking repeat detection (TSD / DR / TIR)
- Starship **captain gene** — tyrosine recombinase (DUF3435) via six-frame ORF
  translation + catalytic-tetrad motif scan
- **Donor lineage** direction and **transfer age** (amelioration: RECENT /
  AMELIORATING / ANCIENT / HOST_NATIVE) from GC / GC3 / dinucleotide signature

All alignment-free — no external tools, HMMs, or reference databases.

## Output
- New VCF INFO tags: `EBND`, `FLANK`, `CAPTAIN`, `CAPTAIN_SC`, `DONOR`, `XFERAGE`
- New GFA S-line tags: `EB`, `FL`, `CP`, `DN`, `XA`

## Tests
- Annotator unit tests (boundary / captain / donor) and VCF header conformance
