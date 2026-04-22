#!/usr/bin/env python3
"""run_te_benchmark.py — TE classification benchmark vs PanTEon SOTA tools.

Compares the MycoSV k-mer nearest-centroid TE classifier against the seven
tools evaluated in the PanTEon paper (Orozco-Arias et al. 2023):
  ClassifyTE, CREATE, DeepTE, NeuralTE, TEClass2, TERL, Terrier

Evaluation protocol matches PanTEon:
  • Balanced test set: up to --max-per-class sequences per Class/Superfamily
  • Metrics: precision, recall, F1 at three taxonomy levels (class, order, superfamily)
  • Each tool is skipped gracefully if its binary is not on $PATH

Usage examples:
  # Train on Dfam/RepBase labeled FASTA, then benchmark on a test FASTA:
  python3 run_te_benchmark.py \\
      --train-fasta train_labeled.fasta \\
      --test-fasta  test_labeled.fasta  \\
      --out-dir     te_benchmark_results/

  # Skip training (use existing index), only classify + compare:
  python3 run_te_benchmark.py \\
      --test-fasta  test_labeled.fasta \\
      --te-index-prefix te_benchmark_results/clf \\
      --out-dir     te_benchmark_results/

  # Download a small PanTEon-style fungal subset from Dfam, then benchmark:
  python3 run_te_benchmark.py --download-fungi-demo --out-dir te_benchmark_results/
"""

import argparse
import csv
import json
import logging
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import textwrap
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger("te_benchmark")

# ---------------------------------------------------------------------------
# Label parsing  (mirrors te_classifier.hpp)
# ---------------------------------------------------------------------------

CLASS_NORM = {
    "Class_I": "LTR", "Retrotransposon": "LTR",
    "Class_II": "DNA", "MITE": "DNA", "TIR": "DNA",
}


@dataclass
class TELabel:
    id: str = ""
    te_class: str = ""
    te_order: str = ""
    superfamily: str = ""
    labeled: bool = False

    def centroid_key(self) -> str:
        if self.te_order:
            return f"{self.te_class}/{self.te_order}/{self.superfamily}"
        return f"{self.te_class}/{self.superfamily}"


def parse_label(header: str) -> TELabel:
    h = header.lstrip(">").strip()
    hash_pos = h.find("#")
    space_pos = h.find(" ")
    id_end = min(
        hash_pos if hash_pos != -1 else len(h),
        space_pos if space_pos != -1 else len(h),
    )
    lbl = TELabel(id=h[:id_end])
    if hash_pos == -1:
        return lbl
    tax = h[hash_pos + 1:]
    if " " in tax:
        tax = tax[:tax.index(" ")]
    parts = [p for p in tax.split("/") if p]
    if not parts:
        return lbl
    lbl.te_class = CLASS_NORM.get(parts[0], parts[0])
    lbl.te_order = parts[1] if len(parts) >= 3 else ""
    lbl.superfamily = parts[-1]
    lbl.labeled = True
    return lbl


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------

def read_fasta(path: str):
    """Yield (header, seq) tuples."""
    with open(path) as f:
        hdr, seq = None, []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq)
                hdr, seq = line, []
            else:
                seq.append(line)
        if hdr is not None:
            yield hdr, "".join(seq)


def write_fasta(path: str, records):
    with open(path, "w") as f:
        for hdr, seq in records:
            h = hdr if hdr.startswith(">") else ">" + hdr
            f.write(h + "\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


# ---------------------------------------------------------------------------
# Benchmark split helpers
# ---------------------------------------------------------------------------

def split_train_test(fasta_path: str, test_fraction: float = 0.2,
                     max_per_class: int = 200, seed: int = 42
                     ) -> Tuple[List, List]:
    """Split labeled FASTA into train/test, balanced per superfamily."""
    rng = random.Random(seed)
    by_sf: Dict[str, List] = defaultdict(list)
    for hdr, seq in read_fasta(fasta_path):
        lbl = parse_label(hdr)
        if lbl.labeled:
            by_sf[lbl.centroid_key()].append((hdr, seq))

    train_recs, test_recs = [], []
    for key, recs in by_sf.items():
        rng.shuffle(recs)
        recs = recs[:max_per_class]
        n_test = max(1, int(len(recs) * test_fraction))
        test_recs.extend(recs[:n_test])
        train_recs.extend(recs[n_test:])

    log.info(f"Split: {len(train_recs)} train, {len(test_recs)} test "
             f"({len(by_sf)} superfamilies)")
    return train_recs, test_recs


# ---------------------------------------------------------------------------
# Demo data: small synthetic fungal TE dataset matching PanTEon classes
# ---------------------------------------------------------------------------

DEMO_SUPERFAMILIES = {
    "LTR/Gypsy/Gypsy":             ("ACGT" * 50, 30),
    "LTR/Copia/Copia":             ("TGCA" * 50, 30),
    "DNA/TIR/Tc1-Mariner":         ("GCTA" * 50, 30),
    "DNA/TIR/hAT":                 ("CATG" * 50, 30),
    "LINE/L1/L1":                  ("ATCG" * 50, 30),
    "SINE/SINE/tRNA":              ("TAGC" * 50, 20),
    "DNA/TIR/PIF-Harbinger":       ("GCAT" * 50, 20),
    "LTR/Gypsy/Chromovirus":       ("CGTA" * 50, 20),
    "DNA/Helitron/Helitron":       ("AGTC" * 50, 20),
    "LINE/RTE/RTE":                ("TCAG" * 50, 20),
}

def generate_demo_fasta(path: str, seed: int = 42):
    """Generate a small synthetic labeled FASTA for smoke-testing."""
    rng = random.Random(seed)
    bases = "ACGT"
    records = []
    for label, (motif, n) in DEMO_SUPERFAMILIES.items():
        parts = label.split("/")
        te_class, te_order, sf = parts[0], parts[1], parts[2]
        for i in range(n):
            # Sequence: motif repeated + 10% random noise
            seq_list = list(motif * 10)
            for j in range(len(seq_list)):
                if rng.random() < 0.1:
                    seq_list[j] = rng.choice(bases)
            seq = "".join(seq_list)
            hdr = f">seq_{te_class}_{sf}_{i:03d}#{te_class}/{te_order}/{sf}"
            records.append((hdr, seq))
    rng.shuffle(records)
    write_fasta(path, records)
    log.info(f"Generated {len(records)} demo sequences → {path}")


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

@dataclass
class Metrics:
    tp: int = 0
    fp: int = 0
    fn: int = 0

    def precision(self) -> float:
        return self.tp / (self.tp + self.fp) if (self.tp + self.fp) > 0 else 0.0

    def recall(self) -> float:
        return self.tp / (self.tp + self.fn) if (self.tp + self.fn) > 0 else 0.0

    def f1(self) -> float:
        p, r = self.precision(), self.recall()
        return 2 * p * r / (p + r) if (p + r) > 0 else 0.0


def compute_metrics(predictions: List[Dict], truth: Dict[str, TELabel],
                    level: str) -> Metrics:
    """level: 'class', 'order', or 'superfamily'."""
    m = Metrics()
    for pred in predictions:
        seq_id = pred["id"]
        if seq_id not in truth:
            continue
        true_lbl = truth[seq_id]
        if level == "class":
            true_val = true_lbl.te_class
            pred_val = pred.get("pred_class", "")
        elif level == "order":
            true_val = true_lbl.te_order or true_lbl.superfamily
            pred_val = pred.get("pred_order", "") or pred.get("pred_superfamily", "")
        else:  # superfamily
            true_val = true_lbl.superfamily
            pred_val = pred.get("pred_superfamily", "")

        if pred_val == true_val:
            m.tp += 1
        else:
            m.fp += 1
            m.fn += 1
    return m


def per_class_f1(predictions: List[Dict], truth: Dict[str, TELabel],
                 level: str) -> Dict[str, float]:
    by_class: Dict[str, Metrics] = defaultdict(Metrics)
    for pred in predictions:
        seq_id = pred["id"]
        if seq_id not in truth:
            continue
        true_lbl = truth[seq_id]
        if level == "class":
            true_val = true_lbl.te_class
            pred_val = pred.get("pred_class", "")
        elif level == "order":
            true_val = true_lbl.te_order or true_lbl.superfamily
            pred_val = pred.get("pred_order", "") or pred.get("pred_superfamily", "")
        else:
            true_val = true_lbl.superfamily
            pred_val = pred.get("pred_superfamily", "")

        m = by_class[true_val]
        if pred_val == true_val:
            m.tp += 1
        else:
            m.fp += 1
            m.fn += 1
    return {k: v.f1() for k, v in by_class.items()}


# ---------------------------------------------------------------------------
# MycoSV TE classifier runner
# ---------------------------------------------------------------------------

def run_mycosv_te(binary: str, train_fasta: Optional[str], test_fasta: str,
                  index_prefix: str, k: int = 21, fracmin_p: float = 0.05,
                  out_tsv: Optional[str] = None) -> Optional[List[Dict]]:
    """Train (if train_fasta) and classify test_fasta. Returns list of prediction dicts."""

    with tempfile.NamedTemporaryFile("w", suffix=".lst", delete=False) as f:
        train_lst = f.name
    with tempfile.NamedTemporaryFile("w", suffix=".lst", delete=False) as f:
        test_lst = f.name

    try:
        # Training step
        if train_fasta:
            Path(train_lst).write_text(train_fasta + "\n")
            cmd = [binary,
                   "--te-train",
                   "--query-list", train_lst,
                   "--te-index-prefix", index_prefix,
                   "--te-k", str(k),
                   "--te-fracmin-p", str(fracmin_p)]
            log.info(f"[mycosv] training: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                log.error(f"[mycosv] train failed:\n{result.stderr}")
                return None
            if result.stderr:
                log.info(f"[mycosv] {result.stderr.strip()}")

        # Classify step
        Path(test_lst).write_text(test_fasta + "\n")
        pred_tsv = out_tsv or (index_prefix + "_predictions.tsv")
        cmd = [binary,
               "--te-classify",
               "--query-list", test_lst,
               "--te-index-prefix", index_prefix,
               "--out-prefix", index_prefix,
               "--te-k", str(k),
               "--te-fracmin-p", str(fracmin_p)]
        log.info(f"[mycosv] classifying: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            log.error(f"[mycosv] classify failed:\n{result.stderr}")
            return None

        # The output is written to index_prefix + ".te_predictions.tsv"
        actual_tsv = index_prefix + ".te_predictions.tsv"
        if not Path(actual_tsv).exists():
            log.error(f"[mycosv] predictions file not found: {actual_tsv}")
            return None

        preds = []
        with open(actual_tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                preds.append(row)
        return preds

    finally:
        for p in [train_lst, test_lst]:
            try:
                os.unlink(p)
            except OSError:
                pass


# ---------------------------------------------------------------------------
# SOTA tool adapters
# The adapters follow PanTEon paper conventions for each tool's I/O format.
# Each returns a list of prediction dicts or None if the tool is unavailable.
# ---------------------------------------------------------------------------

def _which(name: str) -> Optional[str]:
    return shutil.which(name)


def run_repclass(test_fasta: str, out_dir: str) -> Optional[List[Dict]]:
    """RepClass / TEClass2."""
    binary = _which("TEClass2") or _which("repclass")
    if not binary:
        log.debug("[TEClass2/repclass] not found on PATH — skipping")
        return None
    out_tsv = str(Path(out_dir) / "teclass2_predictions.tsv")
    cmd = [binary, "-i", test_fasta, "-o", out_tsv]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning(f"[TEClass2] failed:\n{result.stderr[:500]}")
        return None
    # Parse tab-separated: ID, class, (optional order/sf)
    preds = []
    with open(out_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                preds.append({"id": parts[0], "pred_class": parts[1],
                              "pred_order": parts[2] if len(parts) > 2 else "",
                              "pred_superfamily": parts[3] if len(parts) > 3 else ""})
    return preds


def run_deepte(test_fasta: str, out_dir: str) -> Optional[List[Dict]]:
    """DeepTE."""
    binary = _which("DeepTE.py") or _which("DeepTE")
    if not binary:
        log.debug("[DeepTE] not found on PATH — skipping")
        return None
    out_dir_dt = str(Path(out_dir) / "deepte_out")
    Path(out_dir_dt).mkdir(exist_ok=True)
    cmd = [binary, "-d", out_dir_dt, "-o", out_dir_dt, "-i", test_fasta,
           "--sp", "F"]  # F = fungi
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=out_dir_dt)
    if result.returncode != 0:
        log.warning(f"[DeepTE] failed:\n{result.stderr[:500]}")
        return None
    # DeepTE writes opt_DeepTE.txt
    out_file = Path(out_dir_dt) / "opt_DeepTE.txt"
    if not out_file.exists():
        return None
    preds = []
    with open(out_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                lbl = parse_label("#" + parts[1])
                preds.append({"id": parts[0],
                              "pred_class": lbl.te_class,
                              "pred_order": lbl.te_order,
                              "pred_superfamily": lbl.superfamily})
    return preds


def run_neuralte(test_fasta: str, out_dir: str) -> Optional[List[Dict]]:
    """NeuralTE."""
    binary = _which("NeuralTE.py") or _which("NeuralTE")
    if not binary:
        log.debug("[NeuralTE] not found on PATH — skipping")
        return None
    out_tsv = str(Path(out_dir) / "neuralte_predictions.tsv")
    cmd = [binary, "--genome", test_fasta, "--outdir", str(Path(out_dir) / "neuralte_out"),
           "--species", "fungi"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning(f"[NeuralTE] failed:\n{result.stderr[:500]}")
        return None
    out_file = Path(out_dir) / "neuralte_out" / "classified_TE.fa"
    if not out_file.exists():
        return None
    preds = []
    for hdr, _ in read_fasta(str(out_file)):
        lbl = parse_label(hdr)
        if lbl.labeled:
            preds.append({"id": lbl.id,
                          "pred_class": lbl.te_class,
                          "pred_order": lbl.te_order,
                          "pred_superfamily": lbl.superfamily})
    return preds


def run_terl(test_fasta: str, out_dir: str) -> Optional[List[Dict]]:
    """TERL."""
    binary = _which("TERL.py") or _which("TERL")
    if not binary:
        log.debug("[TERL] not found on PATH — skipping")
        return None
    out_tsv = str(Path(out_dir) / "terl_predictions.tsv")
    cmd = [binary, "-i", test_fasta, "-o", out_tsv]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning(f"[TERL] failed:\n{result.stderr[:500]}")
        return None
    preds = []
    with open(out_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                lbl = parse_label("#" + parts[1])
                preds.append({"id": parts[0],
                              "pred_class": lbl.te_class,
                              "pred_order": lbl.te_order,
                              "pred_superfamily": lbl.superfamily})
    return preds


def run_terrier(test_fasta: str, out_dir: str) -> Optional[List[Dict]]:
    """Terrier."""
    binary = _which("terrier.py") or _which("terrier")
    if not binary:
        log.debug("[Terrier] not found on PATH — skipping")
        return None
    out_file = str(Path(out_dir) / "terrier_predictions.txt")
    cmd = [binary, "-q", test_fasta, "-o", out_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning(f"[Terrier] failed:\n{result.stderr[:500]}")
        return None
    preds = []
    with open(out_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                lbl = parse_label("#" + parts[1])
                preds.append({"id": parts[0],
                              "pred_class": lbl.te_class,
                              "pred_order": lbl.te_order,
                              "pred_superfamily": lbl.superfamily})
    return preds


def run_classifyte(test_fasta: str, out_dir: str) -> Optional[List[Dict]]:
    """ClassifyTE."""
    binary = _which("ClassifyTE.py") or _which("ClassifyTE")
    if not binary:
        log.debug("[ClassifyTE] not found on PATH — skipping")
        return None
    out_file = str(Path(out_dir) / "classifyte_predictions.txt")
    cmd = [binary, "-i", test_fasta, "-o", out_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning(f"[ClassifyTE] failed:\n{result.stderr[:500]}")
        return None
    preds = []
    with open(out_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                lbl = parse_label("#" + parts[1])
                preds.append({"id": parts[0],
                              "pred_class": lbl.te_class,
                              "pred_order": lbl.te_order,
                              "pred_superfamily": lbl.superfamily})
    return preds


def run_create(test_fasta: str, out_dir: str) -> Optional[List[Dict]]:
    """CREATE."""
    binary = _which("CREATE.py") or _which("create")
    if not binary:
        log.debug("[CREATE] not found on PATH — skipping")
        return None
    out_file = str(Path(out_dir) / "create_predictions.txt")
    cmd = [binary, "-i", test_fasta, "-o", out_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning(f"[CREATE] failed:\n{result.stderr[:500]}")
        return None
    preds = []
    with open(out_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                lbl = parse_label("#" + parts[1])
                preds.append({"id": parts[0],
                              "pred_class": lbl.te_class,
                              "pred_order": lbl.te_order,
                              "pred_superfamily": lbl.superfamily})
    return preds


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def format_table(rows: List[Dict], columns: List[str]) -> str:
    widths = {c: len(c) for c in columns}
    for row in rows:
        for c in columns:
            widths[c] = max(widths[c], len(str(row.get(c, ""))))
    header = "  ".join(c.ljust(widths[c]) for c in columns)
    sep    = "  ".join("-" * widths[c] for c in columns)
    lines  = [header, sep]
    for row in rows:
        lines.append("  ".join(str(row.get(c, "")).ljust(widths[c]) for c in columns))
    return "\n".join(lines)


def write_report(results: Dict, out_dir: str):
    """Write JSON + text summary matching PanTEon paper Table 2 / Figure 4 style."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # JSON
    with open(out_path / "te_benchmark_results.json", "w") as f:
        json.dump(results, f, indent=2)

    # Text table
    rows = []
    for tool, tool_res in results.items():
        if not isinstance(tool_res, dict):
            continue
        row = {"Tool": tool}
        for level in ("class", "order", "superfamily"):
            if level in tool_res and isinstance(tool_res[level], dict):
                lv = tool_res[level]
                row[f"{level[0].upper()}F1"] = f"{lv['f1']:.3f}" if "f1" in lv else "-"
                row[f"{level[0].upper()}P"]  = f"{lv['precision']:.3f}" if "precision" in lv else "-"
                row[f"{level[0].upper()}R"]  = f"{lv['recall']:.3f}" if "recall" in lv else "-"
        rows.append(row)

    cols = ["Tool", "CF1", "CP", "CR", "OF1", "OP", "OR", "SF1", "SP", "SR"]
    cols = [c for c in cols if any(c in r for r in rows)]

    report = textwrap.dedent(f"""\
        TE Classification Benchmark — PanTEon-style comparison
        =======================================================
        Test sequences: {results.get('n_test', '?')}
        Superfamilies:  {results.get('n_superfamilies', '?')}
        Training seqs:  {results.get('n_train', '?')}

        C = Class level   O = Order level   S = Superfamily level
        F1/P/R = macro F1 / Precision / Recall

    """)
    report += format_table(rows, cols) + "\n"
    report_path = out_path / "te_benchmark_report.txt"
    report_path.write_text(report)
    log.info(f"Report written to {report_path}")
    print("\n" + report)


# ---------------------------------------------------------------------------
# Main benchmark orchestrator
# ---------------------------------------------------------------------------

def run_benchmark(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Data setup ---
    if args.download_fungi_demo:
        demo_fasta = str(out_dir / "demo_fungi_tes.fasta")
        generate_demo_fasta(demo_fasta, seed=args.seed)
        args.train_fasta = demo_fasta
        args.test_fasta  = demo_fasta

    if not args.test_fasta:
        log.error("--test-fasta is required (or use --download-fungi-demo)")
        sys.exit(1)

    # Determine train/test split
    train_fasta = args.train_fasta
    test_fasta  = args.test_fasta

    if train_fasta == test_fasta or args.auto_split:
        log.info(f"Auto-splitting {test_fasta} into train/test "
                 f"({int((1-args.test_fraction)*100)}% / {int(args.test_fraction*100)}%)")
        train_recs, test_recs = split_train_test(
            test_fasta,
            test_fraction=args.test_fraction,
            max_per_class=args.max_per_class,
            seed=args.seed)
        train_fasta = str(out_dir / "train.fasta")
        test_fasta  = str(out_dir / "test.fasta")
        write_fasta(train_fasta, train_recs)
        write_fasta(test_fasta, test_recs)

    # Build ground truth from test FASTA labels
    truth: Dict[str, TELabel] = {}
    for hdr, _ in read_fasta(test_fasta):
        lbl = parse_label(hdr)
        if lbl.labeled:
            truth[lbl.id] = lbl

    n_test = len(truth)
    n_sf   = len({lbl.superfamily for lbl in truth.values()})
    n_train = sum(1 for _, _ in read_fasta(train_fasta))
    log.info(f"Test set: {n_test} labeled sequences, {n_sf} superfamilies")
    log.info(f"Train set: {n_train} sequences")

    results: Dict = {
        "n_test": n_test,
        "n_train": n_train,
        "n_superfamilies": n_sf,
    }

    # --- MycoSV TE classifier ---
    binary = args.mycosv_binary
    if not Path(binary).exists():
        binary_found = shutil.which(binary)
        if not binary_found:
            log.warning(f"[MycoSV] binary not found: {binary} — building from main.cpp")
            # Attempt to build
            build_cmd = ["g++", "-O2", "-std=c++17", "-pthread",
                         "-I", str(Path(__file__).parent),
                         str(Path(__file__).parent / "main.cpp"),
                         "-o", str(out_dir / "mycosv_te")]
            log.info(f"Building: {' '.join(build_cmd)}")
            res = subprocess.run(build_cmd, capture_output=True, text=True)
            if res.returncode != 0:
                log.error(f"Build failed:\n{res.stderr}")
                binary = None
            else:
                binary = str(out_dir / "mycosv_te")
                log.info(f"Built: {binary}")
        else:
            binary = binary_found

    if binary:
        index_prefix = args.te_index_prefix or str(out_dir / "clf")
        do_train = train_fasta and not (
            Path(index_prefix + ".vptree").exists() and not args.retrain)
        preds = run_mycosv_te(
            binary=binary,
            train_fasta=train_fasta if do_train else None,
            test_fasta=test_fasta,
            index_prefix=index_prefix,
            k=args.te_k,
            fracmin_p=args.te_fracmin_p)
        if preds is not None:
            mycosv_res = {}
            for level in ("class", "order", "superfamily"):
                m = compute_metrics(preds, truth, level)
                mycosv_res[level] = {
                    "f1": m.f1(), "precision": m.precision(), "recall": m.recall(),
                    "tp": m.tp, "fp": m.fp, "fn": m.fn,
                }
            mycosv_res["per_superfamily_f1"] = per_class_f1(preds, truth, "superfamily")
            results["MycoSV"] = mycosv_res
            log.info(f"[MycoSV] class F1={mycosv_res['class']['f1']:.3f}  "
                     f"order F1={mycosv_res['order']['f1']:.3f}  "
                     f"sf F1={mycosv_res['superfamily']['f1']:.3f}")

    # --- SOTA tools: pre-flight availability check ---
    sota_runners = {
        "TEClass2":   (run_repclass,   ["TEClass2", "repclass"]),
        "DeepTE":     (run_deepte,     ["DeepTE.py", "DeepTE"]),
        "NeuralTE":   (run_neuralte,   ["NeuralTE.py", "NeuralTE"]),
        "TERL":       (run_terl,       ["TERL.py", "TERL"]),
        "Terrier":    (run_terrier,    ["terrier.py", "terrier"]),
        "ClassifyTE": (run_classifyte, ["ClassifyTE.py", "ClassifyTE"]),
        "CREATE":     (run_create,     ["CREATE.py", "create"]),
    }
    _INSTALL_SCRIPT = str(Path(__file__).parent / "install_tools.sh")
    if not args.skip_sota:
        _available_te = [n for n, (_, bins) in sota_runners.items()
                         if any(shutil.which(b) for b in bins)]
        _missing_te   = [n for n, (_, bins) in sota_runners.items()
                         if not any(shutil.which(b) for b in bins)]
        if _missing_te:
            log.warning(
                "The following TE classification tools are not installed and will be "
                "skipped (their columns will show 'unavailable' in the report):\n"
                + "".join(f"  ✗ {t}\n" for t in _missing_te)
                + f"\nTo install all tools: bash {_INSTALL_SCRIPT} --te-only\n"
                + f"To see what is missing: bash {_INSTALL_SCRIPT} --check\n"
            )
        if _available_te:
            log.info("SOTA TE tools available: " + ", ".join(_available_te))

    for tool_name, (runner, _bins) in sota_runners.items():
        if args.skip_sota:
            break
        tool_out = str(out_dir / tool_name.lower())
        Path(tool_out).mkdir(exist_ok=True)
        preds = runner(test_fasta, tool_out)
        if preds is None:
            results[tool_name] = "unavailable"
            continue
        tool_res = {}
        for level in ("class", "order", "superfamily"):
            m = compute_metrics(preds, truth, level)
            tool_res[level] = {
                "f1": m.f1(), "precision": m.precision(), "recall": m.recall(),
            }
        results[tool_name] = tool_res
        log.info(f"[{tool_name}] class F1={tool_res['class']['f1']:.3f}  "
                 f"sf F1={tool_res['superfamily']['f1']:.3f}")

    # PanTEon paper reference values (fungi, Table 2 / Figure 4 estimates)
    results["PanTEon_paper_best_fungi"] = {
        "class":       {"f1": 0.88, "note": "NeuralTE in paper (fungi)"},
        "order":       {"f1": 0.79, "note": "NeuralTE in paper (fungi)"},
        "superfamily": {"f1": 0.72, "note": "NeuralTE in paper (fungi)"},
    }

    write_report(results, str(out_dir))

    # Per-superfamily breakdown for MycoSV
    if "MycoSV" in results and isinstance(results["MycoSV"], dict):
        sf_f1 = results["MycoSV"].get("per_superfamily_f1", {})
        if sf_f1:
            print("\nMycoSV per-superfamily F1:")
            for sf, f1 in sorted(sf_f1.items(), key=lambda x: -x[1]):
                print(f"  {sf:<40}  {f1:.3f}")

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="TE classification benchmark — MycoSV vs PanTEon SOTA tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    data = parser.add_argument_group("Data")
    data.add_argument("--train-fasta", help="Labeled FASTA for training (PanTEon/RepBase format)")
    data.add_argument("--test-fasta",  help="Labeled FASTA for evaluation")
    data.add_argument("--download-fungi-demo", action="store_true",
                      help="Generate a small synthetic fungal TE demo dataset")
    data.add_argument("--auto-split", action="store_true",
                      help="Auto-split --test-fasta into train/test subsets")
    data.add_argument("--test-fraction", type=float, default=0.2,
                      help="Fraction of data to hold out for testing (default 0.2)")
    data.add_argument("--max-per-class", type=int, default=200,
                      help="Max sequences per superfamily in balanced set (default 200)")
    data.add_argument("--seed", type=int, default=42)

    clf = parser.add_argument_group("MycoSV classifier")
    clf.add_argument("--mycosv-binary", default="./fungi_graphsv_tol",
                     help="Path to compiled MycoSV binary (default ./fungi_graphsv_tol)")
    clf.add_argument("--te-index-prefix",
                     help="Prefix for TE index files (default: out_dir/clf)")
    clf.add_argument("--te-k", type=int, default=21)
    clf.add_argument("--te-fracmin-p", type=float, default=0.05)
    clf.add_argument("--retrain", action="store_true",
                     help="Force retraining even if index already exists")

    run = parser.add_argument_group("Run options")
    run.add_argument("--out-dir", required=True, help="Output directory")
    run.add_argument("--skip-sota", action="store_true",
                     help="Skip all SOTA tool runners (benchmark MycoSV only)")

    args = parser.parse_args()
    run_benchmark(args)


if __name__ == "__main__":
    main()
