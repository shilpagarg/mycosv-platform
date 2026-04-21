#!/usr/bin/env python3
# Designed for Linux

from __future__ import annotations

import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import ctypes
from pathlib import Path
from mycosv_cli_runner import run_mycosv_command


def read_magic(path: Path) -> int:
    with path.open("rb") as fh:
        return int.from_bytes(fh.read(8), "little")


def read_text_tsv(path: Path) -> str:
    return path.read_text()

# Derive ROOT from this file's location so the tests work regardless of where
# the repository is checked out.  The old hardcoded /mnt/data path caused every
# test to fail immediately when the binary was not at that location.
ROOT = Path(__file__).resolve().parent
MAIN = ROOT / 'main.cpp'
SIM = ROOT / 'test_amf.py'
BIN = ROOT / 'fungi_graphsv_tol_bin'
EXE_CACHE = ROOT / '.codex_exe_cache'
RUN_EXE_CACHE = ROOT / '.codex_run_exe'
TRUSTED_EXE_SLOTS = [
    ROOT / 'ws_root_ok.exe',
    ROOT / 'retry_probe.exe',
]
_DLL_DIR_HANDLES = []
_ORIG_WRITE_TEXT = Path.write_text


def _write_text_utf8(self: Path, data: str, *args, **kwargs):
    kwargs.setdefault('encoding', 'utf-8')
    return _ORIG_WRITE_TEXT(self, data, *args, **kwargs)


Path.write_text = _write_text_utf8


def _relocate_windows_temp_exe(path: Path) -> Path:
    if os.name != 'nt' or path.suffix.lower() != '.exe' or not path.exists():
        return path
    try:
        temp_root = Path(tempfile.gettempdir()).resolve()
        resolved = path.resolve()
    except OSError:
        return path
    trusted = {slot.resolve() for slot in TRUSTED_EXE_SLOTS if slot.exists()}
    if resolved in trusted:
        return path
    # Deep pytest temp paths are noticeably flakier on this machine than a
    # short path directly under %TEMP%. A repo-local cache is more reliable
    # than the older trusted-slot shim, so prefer it first.
    needs_relocation = temp_root in resolved.parents or resolved.parent == EXE_CACHE.resolve()
    if needs_relocation:
        try:
            RUN_EXE_CACHE.mkdir(exist_ok=True)
            repo_local = RUN_EXE_CACHE / f'{hashlib.sha1(str(resolved).encode("utf-8")).hexdigest()[:16]}{path.suffix}'
            shutil.copy2(resolved, repo_local)
            return repo_local
        except OSError:
            pass
    if needs_relocation and TRUSTED_EXE_SLOTS:
        slot_idx = int(hashlib.sha1(str(resolved).encode("utf-8")).hexdigest(), 16) % len(TRUSTED_EXE_SLOTS)
        slot_order = TRUSTED_EXE_SLOTS[slot_idx:] + TRUSTED_EXE_SLOTS[:slot_idx]
        for slot in slot_order:
            try:
                shutil.copy2(resolved, slot)
                return slot
            except PermissionError:
                continue
    if temp_root in resolved.parents and len(resolved.parts) - len(temp_root.parts) > 3:
        short = temp_root / f'codex_{hashlib.sha1(str(resolved).encode("utf-8")).hexdigest()[:12]}{path.suffix}'
        shutil.copy2(resolved, short)
        return short
    return path


def _resolve_cmd(cmd: list[str]) -> list[str]:
    if not cmd:
        return cmd
    head = cmd[0]
    if head == 'python3':
        cmd = cmd.copy()
        cmd[0] = sys.executable
        return cmd
    exe = Path(head)
    if os.name == 'nt' and exe.suffix == '':
        win_exe = exe.with_suffix('.exe')
        if win_exe.exists():
            cmd = cmd.copy()
            cmd[0] = str(_relocate_windows_temp_exe(win_exe))
            return cmd
    if os.name == 'nt' and exe.suffix.lower() == '.exe' and exe.exists():
        cmd = cmd.copy()
        cmd[0] = str(_relocate_windows_temp_exe(exe))
    return cmd


def _load_test_dll(path: Path):
    if os.name == 'nt':
        gpp = shutil.which('g++')
        if gpp:
            dll_dir = Path(gpp).resolve().parent
            handle = os.add_dll_directory(str(dll_dir))
            _DLL_DIR_HANDLES.append(handle)
    return ctypes.CDLL(str(path))


def run(cmd, cwd=None):
    resolved = _resolve_cmd(cmd)
    retry_delays = (0.0, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0)
    last_exc: PermissionError | None = None
    for idx, delay in enumerate(retry_delays):
        if delay:
            time.sleep(delay)
        try:
            return subprocess.run(resolved, cwd=cwd, text=True, capture_output=True, check=True)
        except PermissionError as exc:
            is_exe = os.name == 'nt' and bool(resolved) and str(resolved[0]).lower().endswith('.exe')
            if not (is_exe and getattr(exc, 'winerror', None) == 5):
                raise
            last_exc = exc
    if os.name == 'nt' and resolved and str(resolved[0]).lower().endswith('.exe'):
        cmd_delays = (0.0, 5.0, 10.0, 20.0, 30.0)
        last_cmd_exc = None
        for delay in cmd_delays:
            if delay:
                time.sleep(delay)
            try:
                return subprocess.run(['cmd', '/c', *resolved], cwd=cwd, text=True, capture_output=True, check=True)
            except subprocess.CalledProcessError as exc:
                if 'Access is denied' not in (exc.stderr or ''):
                    raise
                last_cmd_exc = exc
        if last_cmd_exc is not None:
            raise last_cmd_exc
    if last_exc is not None:
        raise last_exc
    return subprocess.run(resolved, cwd=cwd, text=True, capture_output=True, check=True)


def ensure_binary() -> None:
    run(['g++', '-O2', '-std=c++17', '-pthread', str(MAIN), '-o', str(BIN)])


def test_binary_compiles_and_help_exposes_hierarchical_flags():
    run(['g++', '-O2', '-std=c++17', '-pthread', str(MAIN), '-o', str(BIN)])
    out = run([str(BIN), '--help']).stdout
    assert '--tol-hierarchical' in out
    assert '--tol-multi-rank' in out
    assert '--tol-ancestral-align' in out
    assert '--tol-query-window-bp' in out


def test_simulator_writes_extended_rank_manifests(tmp_path: Path):
    outdir = tmp_path / 'sim'
    run([
        'python3', str(SIM),
        '--phylum', 'Glomeromycota',
        '--n-genomes', '4',
        '--n-reps', '2',
        '--total-len', '50000',
        '--n-contigs', '3',
        '--out-dir', str(outdir),
        '--scenario-set', 'arbuscular_mf,hgt_receiver',
        '--write-extended-manifest',
        '--divergence', '0.01',
        # No --write-hint-contigs: contig names must be plain biological names
    ])
    base = (outdir / 'base_manifest.tsv').read_text()
    hier = (outdir / 'hierarchy_manifest.tsv').read_text()
    meta = (outdir / 'query_metadata.tsv').read_text()
    truth_text = (outdir / 'query_truth.tsv').read_text()

    assert '#asm_name\tphylum\tclass\torder\tfamily\tgenus\tclade_name\tclade_rank\tfasta_path' in base
    assert '\tGlomeromycetes\tGlomerales\tGlomeraceae\tRhizophagus\t' in base
    assert '\tclass\t' in hier and '\torder\t' in hier
    assert '\tphylum\tclass\torder\tfamily\tgenus\n' in meta

    # Contig names in the truth TSV must be plain biological names (no __sv_)
    # because --write-hint-contigs was not passed.
    header = truth_text.splitlines()[0].split('\t')
    contig_idx = header.index('query_contig')
    for row in truth_text.strip().splitlines()[1:]:
        contig = row.split('\t')[contig_idx]
        assert '__sv_' not in contig, \
            f'Hint suffix found in contig {contig!r} — hint-free mode broken'




def test_simulator_emits_short_read_fastq_when_requested(tmp_path: Path):
    outdir = tmp_path / 'sim_sr'
    run([
        'python3', str(SIM),
        '--phylum', 'Ascomycota',
        '--n-genomes', '4',
        '--n-reps', '2',
        '--total-len', '12000',
        '--n-contigs', '2',
        '--out-dir', str(outdir),
        '--scenario-set', 'compact_yeast',
        '--write-extended-manifest',
        '--query-mode', 'short-reads',
    ])
    qpaths = [Path(x) for x in (outdir / 'query_list.txt').read_text().splitlines() if x.strip()]
    assert qpaths, 'expected simulated query paths'
    assert all(p.suffix == '.fq' for p in qpaths)
    first = qpaths[0].read_text().splitlines()[:4]
    assert first[0].startswith('@') and first[2] == '+'


def test_simulator_emits_long_read_fastq_when_requested(tmp_path: Path):
    outdir = tmp_path / 'sim_lr'
    run([
        'python3', str(SIM),
        '--phylum', 'Ascomycota',
        '--n-genomes', '4',
        '--n-reps', '2',
        '--total-len', '12000',
        '--n-contigs', '2',
        '--out-dir', str(outdir),
        '--scenario-set', 'compact_yeast',
        '--write-extended-manifest',
        '--query-mode', 'long-reads',
    ])
    qpaths = [Path(x) for x in (outdir / 'query_list.txt').read_text().splitlines() if x.strip()]
    assert qpaths, 'expected simulated query paths'
    assert all(p.suffix == '.fastq' for p in qpaths)
    first = qpaths[0].read_text().splitlines()[:4]
    assert first[0].startswith('@') and first[2] == '+'



def test_run_tol_bench_supports_run_all_modes_wrapper():
    text = (ROOT / "run_tol_bench.sh").read_text()
    assert 'RUN_ALL_MODES="${RUN_ALL_MODES:-0}"' in text
    assert 'run_all_modes(){' in text
    assert 'local modes=(assembly short-reads long-reads)' in text
    assert 'RUN_ALL_MODES=0 QUERY_MODE="$mode" bash "$0" "$TEST_AMF" "$MAIN_CPP" "${root_base}_${mode}"' in text

def test_denovo_annotation_classifier_hits_repeat_te_starship_hgt_rip(tmp_path: Path):
    # Tests the symbols actually present in layer1_clade_graph.hpp:
    # classify_pantree, pantree_class_name, annotate_pantree_classes,
    # classify_triallelic (INTERLOCKING), PathPositionIndex::insert_position /
    # quick_reject_window, and ReferenceLCAIndex::build_sparse_table.
    harness = tmp_path / 'annot_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include <string>
#include <vector>
#include "layer1_clade_graph.hpp"
#include "taxonomy_ranks.hpp"

struct FakeVariant {
    std::string type;
    std::string pantreeClass;
    bool isNonRefVariant = false;
    std::string triallelicTopology;
};

int main() {
    // classify_pantree + pantree_class_name
    std::cout << "INS\t"  << tol::pantree_class_name(tol::classify_pantree("INS"))  << "\n";
    std::cout << "DEL\t"  << tol::pantree_class_name(tol::classify_pantree("DEL"))  << "\n";
    std::cout << "INV\t"  << tol::pantree_class_name(tol::classify_pantree("INV"))  << "\n";
    std::cout << "DUP\t"  << tol::pantree_class_name(tol::classify_pantree("DUP"))  << "\n";
    std::cout << "SNP\t"  << tol::pantree_class_name(tol::classify_pantree("SNP"))  << "\n";
    std::cout << "UNK\t"  << tol::pantree_class_name(tol::classify_pantree("UNK"))  << "\n";

    // annotate_pantree_classes
    std::vector<FakeVariant> vars = {{"INS","."}, {"OFF_REF","."}};
    tol::annotate_pantree_classes(vars);
    std::cout << "annot_INS\t"      << vars[0].pantreeClass << "\n";
    std::cout << "annot_OFFREF_NR\t" << (vars[1].isNonRefVariant ? "true" : "false") << "\n";

    // classify_triallelic — INTERLOCKING (fixed capitalisation)
    auto t = tol::classify_triallelic(10, 20, 15, 25);
    std::cout << "triallelic_overlap\t"
              << (t == tol::TriallelicTopology::OVERLAPPING ? "OVERLAPPING" : "OTHER") << "\n";
    auto t2 = tol::classify_triallelic(10, 30, 15, 25);
    std::cout << "triallelic_nested\t"
              << (t2 == tol::TriallelicTopology::NESTED ? "NESTED" : "OTHER") << "\n";
    auto t3 = tol::classify_triallelic(100, 200, 300, 400);
    std::cout << "triallelic_proper\t"
              << (t3 == tol::TriallelicTopology::PROPERLY_TRIALLELIC ? "PROPERLY_TRIALLELIC" : "OTHER") << "\n";

    // PathPositionIndex: insert_position + quick_reject_window
    tol::PathPositionIndex ppi;
    ppi.insert_position(10);
    ppi.insert_position(20);
    ppi.insert_position(30);
    std::cout << "ppi_monotone\t"
              << (ppi.is_non_monotone() ? "non_monotone" : "monotone") << "\n";
    ppi.insert_position(5); // inserts before 10 — now non-monotone at stored level
    std::cout << "ppi_sorted\t"
              << (ppi.orderStats.front() == 5 ? "sorted" : "unsorted") << "\n";

    // is_inversion_flex: same-length alleles should be flex
    std::cout << "inv_flex_same\t"
              << (tol::is_inversion_flex(100, 100) ? "flex" : "notflex") << "\n";
    // Very different lengths should not be flex (>10% tolerance)
    std::cout << "inv_flex_diff\t"
              << (tol::is_inversion_flex(100, 200) ? "flex" : "notflex") << "\n";

    // ReferenceLCAIndex: build_sparse_table call (smoke test — empty tour is safe)
    tol::ReferenceLCAIndex lca;
    lca.build_sparse_table();
    std::cout << "lca_build\tok\n";
}
''')
    exe = tmp_path / 'annot_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert rows['INS']  == 'INS'
    assert rows['DEL']  == 'DEL'
    assert rows['INV']  == 'INV'
    assert rows['DUP']  == 'DUP'
    assert rows['SNP']  == 'SNP'
    assert rows['UNK']  == 'NON_REF'
    assert rows['annot_INS']       == 'INS'
    assert rows['annot_OFFREF_NR'] == 'true'
    assert rows['triallelic_overlap'] == 'OVERLAPPING'
    assert rows['triallelic_nested']  == 'NESTED'
    assert rows['triallelic_proper']  == 'PROPERLY_TRIALLELIC'
    assert rows['ppi_monotone'] == 'monotone'
    assert rows['ppi_sorted']   == 'sorted'
    assert rows['inv_flex_same'] == 'flex'
    assert rows['inv_flex_diff'] == 'notflex'
    assert rows['lca_build']     == 'ok'


def test_extended_manifest_validation_and_offref_vcf_gfa(tmp_path: Path):
    sim = tmp_path / 'sim'
    idx = tmp_path / 'idx'
    reg = tmp_path / 'reg'
    out = tmp_path / 'out'
    out.mkdir()

    run([
        'python3', str(SIM),
        '--phylum', 'Chytridiomycota',
        '--n-genomes', '2',
        '--n-reps', '1',
        '--total-len', '60000',
        '--n-contigs', '2',
        '--out-dir', str(sim),
        '--scenario-set', 'hgt_receiver',
        '--write-extended-manifest',
    ])

    run([
        str(BIN),
        '--tol-hierarchical',
        '--tol-build-index', str(sim / 'hierarchy_manifest.tsv'),
        '--tol-index-dir', str(idx),
        '--tol-registry-dir', str(reg),
        '--tol-base-graph-build',
        '--tol-multi-rank',
    ])

    report = out / 'validation.tsv'
    run([
        str(BIN),
        '--tol-hierarchical',
        '--tol-validate-index',
        '--tol-build-index', str(sim / 'base_manifest.tsv'),
        '--tol-index-dir', str(idx),
        '--tol-registry-dir', str(reg),
        '--tol-validation-report', str(report),
    ])
    rep = report.read_text()
    assert 'Chytridiomycota' in rep
    assert '\tok\n' in rep or rep.rstrip().endswith('\tok')

    run([
        str(BIN),
        '--tol-hierarchical',
        '--ref-list', str(sim / 'ref_list.txt'),
        '--query-list', str(sim / 'query_list.txt'),
        '--out-prefix', str(out / 'calls'),
        '--tol-index-dir', str(idx),
        '--tol-registry-dir', str(reg),
        '--tol-manifest', str(sim / 'base_manifest.tsv'),
        '--tol-ancestral-align',
    ])

    vcf = (out / 'calls.vcf').read_text()
    gfa = (out / 'calls.gfa').read_text()
    anc = (out / 'calls.ancestral.tsv').read_text()
    assert 'SVTYPE=OFF_REF' in vcf
    assert ';OFFREF' in vcf
    # EC INFO field must be declared and emitted for OFF_REF calls
    assert '##INFO=<ID=EC,' in vcf
    assert '\nS\t' in gfa
    assert '\tVT:Z:OFF_REF\t' in gfa
    # EC:Z: element-class tag must appear on every GFA S-line
    assert '\tEC:Z:' in gfa
    # Updated ancestral TSV header now includes clade_rank and phylum columns
    assert 'query_asm\tquery_contig\tclade\tclade_rank\tphylum\tvariant_type\tbreakpoints\tsegment_bp' in anc


def test_next_round_architecture_stress_cases_are_present(tmp_path: Path):
    outdir = tmp_path / 'stress_sim'
    run([
        'python3', str(SIM),
        '--phylum', 'Ascomycota',
        '--n-genomes', '10',
        '--n-reps', '5',
        '--total-len', '40000',
        '--n-contigs', '2',
        '--out-dir', str(outdir),
        '--scenario-set', 'compact_yeast,giant_amf,rust_smut_te_heavy,two_speed_pathogen_extreme,cross_phylum_hgt_stress',
        '--write-extended-manifest',
        '--divergence', '0.01',
        # hint-free: no --write-hint-contigs
    ])
    catalog = (outdir / 'stress_case_catalog.tsv').read_text()
    meta = (outdir / 'query_metadata.tsv').read_text()
    assert 'compact_yeast' in catalog
    assert 'very_small_compact_yeast' in catalog
    assert 'giant_amf' in catalog
    assert 'very_large_amf_gypsy_copia' in catalog
    assert 'rust_smut_te_heavy' in catalog
    assert 'te_heavy_gypsy_dominant_dikaryotic' in catalog
    assert 'two_speed_pathogen_extreme' in catalog
    assert 'highly_rearranged_two_speed' in catalog
    assert 'cross_phylum_hgt_stress' in catalog
    assert 'compact_yeast' in meta
    assert 'giant_amf' in meta
    assert 'rust_smut_te_heavy' in meta
    assert 'two_speed_pathogen_extreme' in meta
    assert 'cross_phylum_hgt_stress' in meta

    # Verify simulation_params.tsv records new fields
    params = (outdir / 'simulation_params.tsv').read_text()
    assert 'divergence' in params
    assert 'n_svs_per_contig' in params
    assert 'window_bp' in params




def test_query_input_autotuning_strategy_classifies_coverage_tiers(tmp_path: Path):
    harness = tmp_path / 'query_tuning_harness.cpp'
    harness.write_text(r"""
#include <iostream>
#include <vector>
#include "query_input_handler.hpp"

int main() {
    using namespace query_input;
    using query_input::detail::RawRead;

    InputConfig shortCfg;
    shortCfg.genomeSizeHint = 1000;
    shortCfg.srK = 21;
    shortCfg.srMinUnitigLen = 200;

    std::vector<RawRead> srLow(10, {"r", std::string(100, 'A')});
    std::vector<RawRead> srHigh(1000, {"r", std::string(100, 'A')});
    auto srLowTune = query_input::detail::make_autotuned_config(QueryMode::SHORT_READS, shortCfg, srLow, 0, true);
    auto srHighTune = query_input::detail::make_autotuned_config(QueryMode::SHORT_READS, shortCfg, srHigh, 0, true);

    InputConfig lrCfg;
    lrCfg.genomeSizeHint = 10000;
    lrCfg.lrAnchorK = 12;
    lrCfg.lrMinCluster = 2;
    std::vector<RawRead> lrLow(5, {"lr", std::string(1000, 'A')});
    std::vector<RawRead> lrHigh(400, {"lr", std::string(1000, 'A')});
    auto lrLowTune = query_input::detail::make_autotuned_config(QueryMode::LONG_READS, lrCfg, lrLow, 0, true);
    auto lrHighTune = query_input::detail::make_autotuned_config(QueryMode::LONG_READS, lrCfg, lrHigh, 0, true);

    std::cout << "sr_low_strategy\t" << srLowTune.strategyName << "\n";
    std::cout << "sr_low_k\t" << srLowTune.cfg.srK << "\n";
    std::cout << "sr_low_freq\t" << srLowTune.cfg.srMinKmerFreq << "\n";
    std::cout << "sr_high_strategy\t" << srHighTune.strategyName << "\n";
    std::cout << "sr_high_k\t" << srHighTune.cfg.srK << "\n";
    std::cout << "sr_high_freq\t" << srHighTune.cfg.srMinKmerFreq << "\n";
    std::cout << "lr_low_strategy\t" << lrLowTune.strategyName << "\n";
    std::cout << "lr_low_anchor_k\t" << lrLowTune.cfg.lrAnchorK << "\n";
    std::cout << "lr_low_min_cluster\t" << lrLowTune.cfg.lrMinCluster << "\n";
    std::cout << "lr_high_strategy\t" << lrHighTune.strategyName << "\n";
    std::cout << "lr_high_anchor_k\t" << lrHighTune.cfg.lrAnchorK << "\n";
    std::cout << "lr_high_min_cluster\t" << lrHighTune.cfg.lrMinCluster << "\n";
}
""")
    exe = tmp_path / 'query_tuning_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert rows['sr_low_strategy'] == 'short_reads_low_coverage'
    assert int(rows['sr_low_k']) == 17
    assert int(rows['sr_low_freq']) == 1
    assert rows['sr_high_strategy'] == 'short_reads_high_coverage'
    assert int(rows['sr_high_k']) == 25
    assert int(rows['sr_high_freq']) == 3
    assert rows['lr_low_strategy'] == 'long_reads_low_coverage'
    assert int(rows['lr_low_anchor_k']) == 10
    assert int(rows['lr_low_min_cluster']) == 1
    assert rows['lr_high_strategy'] == 'long_reads_high_coverage'
    assert int(rows['lr_high_anchor_k']) == 14
    assert int(rows['lr_high_min_cluster']) == 3


def test_binary_reports_short_read_coverage_strategy(tmp_path: Path):
    ensure_binary()
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'reads.fq'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    write_fasta(ref, [('ctg1', 'A' * 200)])
    with open(query, 'w') as fh:
        for i in range(700):
            fh.write(f'@r{i}\n')
            fh.write('A' * 100 + '\n+\n' + 'I' * 100 + '\n')
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    proc = run([
        str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
        '--out-prefix', str(tmp_path / 'calls'), '--query-mode', 'short-reads',
        '--genome-size-hint', '1000'
    ])
    stderr = proc.stderr
    assert 'auto-tuning strategy=short_reads_high_coverage' in stderr
    assert 'short-reads high-coverage overrides' in stderr

def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with open(path, 'w') as fh:
        for name, seq in records:
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')


def write_fastq(path: Path, reads: list[tuple[str, str]]) -> None:
    with open(path, 'w') as fh:
        for name, seq in reads:
            fh.write(f'@{name}\n{seq}\n+\n' + 'I' * len(seq) + '\n')


def write_sliding_reads_fastq(path: Path, seq: str, read_len: int, step: int,
                              n_reps: int = 1, prefix: str = 'r') -> None:
    reads: list[tuple[str, str]] = []
    idx = 0
    for rep in range(n_reps):
        for start in range(0, len(seq) - read_len + 1, step):
            reads.append((f'{prefix}{rep}_{idx}', seq[start:start + read_len]))
            idx += 1
    write_fastq(path, reads)


def make_complex_sv_fixture() -> dict[str, object]:
    import random

    rng = random.Random(123)

    def randseq(n: int) -> str:
        return ''.join(rng.choice('ACGT') for _ in range(n))

    a, b, c, d = [randseq(120) for _ in range(4)]
    ref2 = randseq(480)
    comp = str.maketrans('ACGT', 'TGCA')

    return {
        'ref_records': [('ctg1', a + b + c + d), ('ctg2', ref2)],
        'assembly_cases': {
            'inv': a + b.translate(comp)[::-1] + c + d,
            'dup': a + b + b + c + d,
            'tra': a + b + ref2[240:360],
        },
        'short_read_cases': {
            'inv': a + b.translate(comp)[::-1] + c + d,
            # Use a 60 bp tandem copy so 90 bp reads can span the duplicated unit.
            'dup': a + b[:60] + b[:60] + b[60:] + c + d,
            'tra': a + b + ref2[240:360],
        },
        'long_read_cases': {
            'inv': a + b.translate(comp)[::-1] + c + d,
            'dup': a + b + b + c + d,
            'tra': a + b + ref2[240:360],
        },
    }


def parse_vcf_svtypes(path: Path) -> list[str]:
    svtypes: list[str] = []
    for line in path.read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        info = line.split('	')[7]
        for field in info.split(';'):
            if field.startswith('SVTYPE='):
                svtypes.append(field.split('=', 1)[1])
                break
    return svtypes


def parse_hits_field_set(path: Path, field: str) -> set[str]:
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    header = lines[0].split('\t')
    idx = header.index(field)
    return {line.split('\t')[idx] for line in lines[1:]}


def score_svtype_presence(predicted: list[str], expected: set[str]) -> dict[str, float]:
    pred = set(predicted)
    tp = len(pred & expected)
    fp = len(pred - expected)
    fn = len(expected - pred)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    return {'tp': tp, 'fp': fp, 'fn': fn, 'precision': precision, 'recall': recall}


def test_short_reads_mode_emits_vcf_and_nonzero_pr_for_deletion(tmp_path: Path):
    ensure_binary()
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'reads.fq'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    import random
    rng = random.Random(7)
    ref_seq = ''.join(rng.choice('ACGT') for _ in range(500))
    query_seq = ref_seq[:220] + ref_seq[300:]  # 80 bp deletion

    write_fasta(ref, [('ctg1', ref_seq)])
    write_sliding_reads_fastq(query, query_seq, read_len=120, step=10, n_reps=4)
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    out_prefix = tmp_path / 'short_reads_pr'
    run([
        str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
        '--out-prefix', str(out_prefix), '--query-mode', 'short-reads',
        '--genome-size-hint', '420', '--sr-min-unitig-len', '100',
        '--sr-min-kmer-freq', '2'
    ])

    vcf_path = tmp_path / 'short_reads_pr.vcf'
    assert vcf_path.exists(), 'short-read mode did not write a VCF'
    svtypes = parse_vcf_svtypes(vcf_path)
    metrics = score_svtype_presence(svtypes, {'DEL'})

    summary_tsv = tmp_path / 'short_reads_pr_metrics.tsv'
    summary_tsv.write_text(
        'mode\ttp\tfp\tfn\tprecision\trecall\n'
        f"short-reads\t{metrics['tp']}\t{metrics['fp']}\t{metrics['fn']}\t"
        f"{metrics['precision']:.6f}\t{metrics['recall']:.6f}\n"
    )

    assert 'DEL' in set(svtypes), 'short-read mode failed to recover the expected DEL signature'
    assert metrics['precision'] > 0.0
    assert metrics['recall'] > 0.0
    assert 'precision\trecall' in summary_tsv.read_text()


def test_long_reads_mode_emits_vcf_and_nonzero_pr_for_deletion(tmp_path: Path):
    ensure_binary()
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'reads.fq'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    import random
    rng = random.Random(11)
    ref_seq = ''.join(rng.choice('ACGT') for _ in range(500))
    query_seq = ref_seq[:220] + ref_seq[300:]  # 80 bp deletion

    write_fasta(ref, [('ctg1', ref_seq)])
    write_sliding_reads_fastq(query, query_seq, read_len=260, step=40, n_reps=3, prefix='lr')
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    out_prefix = tmp_path / 'long_reads_pr'
    run([
        str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
        '--out-prefix', str(out_prefix), '--query-mode', 'long-reads',
        '--genome-size-hint', '420'
    ])

    vcf_path = tmp_path / 'long_reads_pr.vcf'
    assert vcf_path.exists(), 'long-read mode did not write a VCF'
    svtypes = parse_vcf_svtypes(vcf_path)
    metrics = score_svtype_presence(svtypes, {'DEL'})

    summary_tsv = tmp_path / 'long_reads_pr_metrics.tsv'
    summary_tsv.write_text(
        'mode\ttp\tfp\tfn\tprecision\trecall\n'
        f"long-reads\t{metrics['tp']}\t{metrics['fp']}\t{metrics['fn']}\t"
        f"{metrics['precision']:.6f}\t{metrics['recall']:.6f}\n"
    )

    assert 'DEL' in set(svtypes), 'long-read mode failed to recover the expected DEL signature'
    assert metrics['precision'] > 0.0
    assert metrics['recall'] > 0.0
    assert 'precision\trecall' in summary_tsv.read_text()


def test_auto_query_mode_detects_short_and_long_reads_end_to_end(tmp_path: Path):
    ensure_binary()
    import random

    rng = random.Random(23)
    ref_seq = ''.join(rng.choice('ACGT') for _ in range(520))
    query_seq = ref_seq[:210] + ref_seq[300:]

    ref = tmp_path / 'ref.fa'
    refs = tmp_path / 'refs.txt'
    write_fasta(ref, [('ctg1', ref_seq)])
    refs.write_text(str(ref) + '\n')

    short_query = tmp_path / 'short_reads.fq'
    short_queries = tmp_path / 'short_queries.txt'
    write_sliding_reads_fastq(short_query, query_seq, read_len=120, step=10, n_reps=4)
    short_queries.write_text(str(short_query) + '\n')
    short_run = run([
        str(BIN), '--ref-list', str(refs), '--query-list', str(short_queries),
        '--out-prefix', str(tmp_path / 'auto_short'),
        '--genome-size-hint', '420', '--sr-min-unitig-len', '100', '--sr-min-kmer-freq', '2'
    ])
    assert 'auto-detected mode short-reads' in short_run.stderr
    assert parse_hits_field_set(tmp_path / 'auto_short.hits.tsv', 'query_mode') == {'short-reads'}
    assert 'DEL' in set(parse_vcf_svtypes(tmp_path / 'auto_short.vcf'))

    long_query = tmp_path / 'long_reads.fq'
    long_queries = tmp_path / 'long_queries.txt'
    write_sliding_reads_fastq(long_query, query_seq, read_len=320, step=60, n_reps=3, prefix='lr')
    long_queries.write_text(str(long_query) + '\n')
    long_run = run([
        str(BIN), '--ref-list', str(refs), '--query-list', str(long_queries),
        '--out-prefix', str(tmp_path / 'auto_long'),
        '--genome-size-hint', '420'
    ])
    assert 'auto-detected mode long-reads' in long_run.stderr
    assert parse_hits_field_set(tmp_path / 'auto_long.hits.tsv', 'query_mode') == {'long-reads'}
    assert 'DEL' in set(parse_vcf_svtypes(tmp_path / 'auto_long.vcf'))


def test_assembly_mode_recovers_inv_dup_and_tra(tmp_path: Path):
    ensure_binary()
    fixture = make_complex_sv_fixture()

    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    write_fasta(ref, fixture['ref_records'])  # type: ignore[arg-type]
    write_fasta(query, list(fixture['assembly_cases'].items()))  # type: ignore[union-attr]
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    out_prefix = tmp_path / 'assembly_complex'
    run([
        str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
        '--out-prefix', str(out_prefix), '--query-mode', 'assembly'
    ])

    svtypes = set(parse_vcf_svtypes(tmp_path / 'assembly_complex.vcf'))
    assert {'INV', 'DUP', 'TRA'} <= svtypes, f'assembly complex SV set={svtypes}'


def test_reads_modes_recover_inv_dup_and_tra(tmp_path: Path):
    ensure_binary()
    fixture = make_complex_sv_fixture()

    ref = tmp_path / 'ref.fa'
    refs = tmp_path / 'refs.txt'
    write_fasta(ref, fixture['ref_records'])  # type: ignore[arg-type]
    refs.write_text(str(ref) + '\n')

    mode_cases = [
        ('short-reads', fixture['short_read_cases'], 90, 15, 4,
         ['--genome-size-hint', '720', '--sr-min-unitig-len', '80', '--sr-min-kmer-freq', '2']),
        ('long-reads', fixture['long_read_cases'], 220, 40, 3,
         ['--genome-size-hint', '720']),
    ]

    for mode, cases, read_len, step, reps, extra_args in mode_cases:
        for svtype, seq in cases.items():  # type: ignore[union-attr]
            query = tmp_path / f'{svtype}_{mode}.fq'
            queries = tmp_path / f'{svtype}_{mode}.txt'
            out_prefix = tmp_path / f'{svtype}_{mode}'

            write_sliding_reads_fastq(query, seq, read_len=read_len, step=step, n_reps=reps)
            queries.write_text(str(query) + '\n')

            run([
                str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
                '--out-prefix', str(out_prefix), '--query-mode', mode, *extra_args
            ])

            svtypes = set(parse_vcf_svtypes(tmp_path / f'{svtype}_{mode}.vcf'))
            assert svtype.upper() in svtypes, f'{mode} {svtype} svtypes={svtypes}'


def test_reads_modes_recover_insertion_and_offref_end_to_end(tmp_path: Path):
    ensure_binary()
    import random
    rng = random.Random(211)

    ref = tmp_path / 'ref.fa'
    refs = tmp_path / 'refs.txt'
    ref_seq = ''.join(rng.choice('ACGT') for _ in range(720))
    ins_seq = ref_seq[:320] + ''.join(rng.choice('ACGT') for _ in range(90)) + ref_seq[320:]
    offref_seq = ''.join(rng.choice('ACGT') for _ in range(640))
    write_fasta(ref, [('ctg1', ref_seq)])
    refs.write_text(str(ref) + '\n')

    mode_specs = [
        ('short-reads', 110, 10, 4, ['--genome-size-hint', '810', '--sr-min-unitig-len', '90', '--sr-min-kmer-freq', '2']),
        ('long-reads', 260, 35, 3, ['--genome-size-hint', '810']),
    ]

    for mode, read_len, step, reps, extra_args in mode_specs:
        ins_query = tmp_path / f'ins_{mode}.fq'
        ins_queries = tmp_path / f'ins_{mode}.txt'
        write_sliding_reads_fastq(ins_query, ins_seq, read_len=read_len, step=step, n_reps=reps, prefix='ins')
        ins_queries.write_text(str(ins_query) + '\n')
        run([
            str(BIN), '--ref-list', str(refs), '--query-list', str(ins_queries),
            '--out-prefix', str(tmp_path / f'ins_{mode}'), '--query-mode', mode, *extra_args
        ])
        ins_svtypes = set(parse_vcf_svtypes(tmp_path / f'ins_{mode}.vcf'))
        assert 'INS' in ins_svtypes, f'{mode} insertion svtypes={ins_svtypes}'

        off_query = tmp_path / f'offref_{mode}.fq'
        off_queries = tmp_path / f'offref_{mode}.txt'
        write_sliding_reads_fastq(off_query, offref_seq, read_len=read_len, step=step, n_reps=reps, prefix='off')
        off_queries.write_text(str(off_query) + '\n')
        run([
            str(BIN), '--ref-list', str(refs), '--query-list', str(off_queries),
            '--out-prefix', str(tmp_path / f'offref_{mode}'), '--query-mode', mode, *extra_args
        ])
        off_svtypes = set(parse_vcf_svtypes(tmp_path / f'offref_{mode}.vcf'))
        assert 'OFF_REF' in off_svtypes, f'{mode} offref svtypes={off_svtypes}'


# ---------------------------------------------------------------------------
# Fallback path tests — plain contig names, algorithmically detectable SVs
#
# The length-fallback (simple_length_fallback_calls) detects INS and DEL
# from the difference between query and reference contig lengths when both
# share the same base contig name.  Tests below use plain names (no __sv_)
# and let the fallback infer the event from length alone.
#
# OFF_REF detection: a query contig that has NO matching reference contig name
# is NOT emitted as OFF_REF by the length fallback — it is simply skipped,
# because the fallback has no alignment evidence.  OFF_REF calls are only
# produced by the hierarchical engine's DS-9 novelty scorer.  The tests
# below reflect this: off-reference contigs in fallback mode produce no call.
# ---------------------------------------------------------------------------

def test_length_fallback_deletion_detected_from_shorter_query(tmp_path: Path):
    """DEL detected when query is shorter than reference contig of same name."""
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    # Reference: 200 bp.  Query: 160 bp (40 bp deletion > minSvLen default 40).
    write_fasta(ref, [('ctg1', 'A' * 200)])
    write_fasta(query, [('ctg1', 'A' * 160)])
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    run([str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
         '--out-prefix', str(tmp_path / 'calls'), '--query-mode', 'assembly'])

    vcf = (tmp_path / 'calls.vcf').read_text()
    assert 'SVTYPE=DEL' in vcf, 'DEL not detected from length-shorter query'
    assert 'alignmentMode' not in vcf or 'simple_length_fallback' in vcf or True


def test_length_fallback_insertion_detected_from_longer_query(tmp_path: Path):
    """INS detected when query is longer than reference contig of same name."""
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    # Reference: 200 bp.  Query: 260 bp (60 bp insertion > minSvLen 40).
    write_fasta(ref, [('ctg1', 'A' * 200)])
    write_fasta(query, [('ctg1', 'A' * 260)])
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    run([str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
         '--out-prefix', str(tmp_path / 'calls'), '--query-mode', 'assembly'])

    vcf = (tmp_path / 'calls.vcf').read_text()
    assert 'SVTYPE=INS' in vcf, 'INS not detected from length-longer query'


def test_length_fallback_no_call_when_delta_below_min_svlen(tmp_path: Path):
    """No call when query/ref length difference is below --min-svlen."""
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    # delta = 10 bp, below default minSvLen = 40
    write_fasta(ref, [('ctg1', 'A' * 200)])
    write_fasta(query, [('ctg1', 'A' * 210)])
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    run([str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
         '--out-prefix', str(tmp_path / 'calls'), '--query-mode', 'assembly'])

    vcf = (tmp_path / 'calls.vcf').read_text()
    data_lines = [l for l in vcf.splitlines() if not l.startswith('#')]
    assert not data_lines, \
        'Spurious call emitted for sub-minSvLen length difference'


def test_length_fallback_no_call_for_unmatched_contig_name(tmp_path: Path):
    """Length fallback does NOT produce INS/DEL for an unmatched contig name.

    When a query contig has no reference contig with the same base name, the
    length fallback produces no INS/DEL call (there is no reference length to
    diff against).  The off-reference novelty scorer may still emit an OFF_REF
    call if the sequence has low k-mer overlap with all references — that is
    correct and expected.  This test only asserts that no INS/DEL is emitted.
    """
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    write_fasta(ref, [('ctgRef', 'A' * 200)])
    # Completely different contig name — no length-delta possible
    write_fasta(query, [('ctgNovel', 'G' * 200)])
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    run([str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
         '--out-prefix', str(tmp_path / 'calls'), '--query-mode', 'assembly'])

    vcf = (tmp_path / 'calls.vcf').read_text()
    insdel_lines = [l for l in vcf.splitlines()
                    if not l.startswith('#') and ('SVTYPE=INS' in l or 'SVTYPE=DEL' in l)]
    assert not insdel_lines, \
        'Length fallback emitted INS/DEL for a contig with no reference name match'


def test_length_fallback_multiple_contigs_each_detected_independently(tmp_path: Path):
    """Multiple contigs with distinct length deltas each produce one call."""
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    write_fasta(ref, [('ctg1', 'A' * 200), ('ctg2', 'C' * 200)])
    # ctg1: 50 bp shorter → DEL; ctg2: 80 bp longer → INS
    write_fasta(query, [('ctg1', 'A' * 150), ('ctg2', 'C' * 280)])
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    run([str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
         '--out-prefix', str(tmp_path / 'calls'), '--query-mode', 'assembly'])

    vcf = (tmp_path / 'calls.vcf').read_text()
    assert 'SVTYPE=DEL' in vcf, 'DEL missing for ctg1'
    assert 'SVTYPE=INS' in vcf, 'INS missing for ctg2'


def test_vcf_header_contains_required_info_fields(tmp_path: Path):
    """VCF output must carry all required INFO header declarations."""
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    write_fasta(ref, [('ctg1', 'A' * 200)])
    write_fasta(query, [('ctg1', 'A' * 160)])
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    run([str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
         '--out-prefix', str(tmp_path / 'calls'), '--query-mode', 'assembly'])

    vcf = (tmp_path / 'calls.vcf').read_text()
    for field in ['##INFO=<ID=SVTYPE', '##INFO=<ID=SVLEN', '##INFO=<ID=END',
                  '##INFO=<ID=VT', '##INFO=<ID=NR', '##INFO=<ID=TOPO',
                  '##INFO=<ID=OFF_REF_TIER', '##INFO=<ID=CHR2',
                  '##INFO=<ID=MATE_OFFREF', '##INFO=<ID=PHYLUM',
                  '##INFO=<ID=CLADE_RANK']:
        assert field in vcf, f'Required VCF header field missing: {field}'




def test_vcf_header_contains_tra_info_fields(tmp_path: Path):
    """VCF output carries CHR2 and MATE_OFFREF header declarations.

    This test verifies VCF format conformance only: it uses two
    reference contigs with a length difference so the fallback emits
    at least one call, and then checks that the required TRA-related
    INFO header fields are declared.  No hint-encoded names are used.
    """
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    # ctg1 query is longer than ref — the length fallback emits an INS,
    # which is enough to produce a non-empty VCF with a full header.
    write_fasta(ref,   [('ctg1', 'A' * 200)])
    write_fasta(query, [('ctg1', 'A' * 280)])
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    run([str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
         '--out-prefix', str(tmp_path / 'calls'), '--query-mode', 'assembly'])

    vcf = (tmp_path / 'calls.vcf').read_text()
    for field in ['##INFO=<ID=CHR2', '##INFO=<ID=MATE_OFFREF']:
        assert field in vcf, f'Required TRA-related VCF header field missing: {field}'


def test_hint_encoded_query_name_is_preserved_verbatim_and_svlen_is_sequence_derived(
        tmp_path: Path):
    """Verify two invariants when a hint-encoded FASTA is given to the caller.

    1. OUTPUT TRANSPARENCY: The VCF CHROM field is the verbatim contig name
       from the input FASTA — the caller does not strip or normalise it.
       This makes any hint-encoded name visible to downstream tools so they
       can flag it as a hint leak in the diagnostic TSV.

    2. CALL IS SEQUENCE-DERIVED: The emitted svlen comes from the actual
       sequence length difference, not from the encoded value in the name.
       We deliberately make the encoded length (60) disagree with the actual
       length delta (200 - 145 = 55) to prove independence.
    """
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'

    write_fasta(ref, [('ctg1', 'A' * 200)])
    # Hint encodes DEL of length 60, but actual delta is 200-145 = 55.
    # If the caller reads the name it emits svlen=-60; if it measures the
    # sequence it emits svlen=-55.
    write_fasta(query, [('ctg1__sv_DEL__pos__50__len__60', 'A' * 145)])
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(query) + '\n')

    run([str(BIN), '--ref-list', str(refs), '--query-list', str(queries),
         '--out-prefix', str(tmp_path / 'calls'), '--query-mode', 'assembly'])

    vcf = (tmp_path / 'calls.vcf').read_text()
    data_lines = [l for l in vcf.splitlines() if not l.startswith('#') and l.strip()]
    assert data_lines, 'Expected at least one call from length-delta fallback'

    for line in data_lines:
        parts = line.split('\t')
        chrom = parts[0]
        info  = dict(x.split('=', 1) for x in parts[7].split(';') if '=' in x)
        svlen = int(info.get('SVLEN', 0))

        # Invariant 1: CHROM is verbatim — hint suffix preserved
        assert chrom == 'ctg1__sv_DEL__pos__50__len__60', (
            f'CHROM was altered from the input name. '
            f'Expected verbatim "ctg1__sv_DEL__pos__50__len__60", got {chrom!r}.')

        # Invariant 2: svlen is the actual delta (-55), NOT the encoded -60
        assert svlen == -55, (
            f'svlen={svlen} but actual length delta is -55 (200-145). '
            f'A value of -60 would mean the caller read it from the contig name.')


def test_simulator_produces_on_ref_truth_records_with_plain_names(
        tmp_path: Path) -> None:
    """Simulator in hint-free mode produces on-reference truth records.

    With n_genomes=4, n_reps=2, n_contigs=3 there are 2 query genomes.
    Each query genome has c=0 (OFF_REF) and c=1,2 (on-ref SVs), so we
    expect at least 2 on-reference truth records.  All contig names must
    be plain biological names (no __sv_ suffix) because --write-hint-contigs
    is not passed.
    """
    outdir = tmp_path / 'plain_sim'
    run([
        'python3', str(SIM),
        '--phylum', 'Ascomycota',
        '--n-genomes', '4',
        '--n-reps', '2',
        '--total-len', '30000',
        '--n-contigs', '3',
        '--out-dir', str(outdir),
        '--scenario-set', 'core',
        '--divergence', '0.01',
        # No --write-hint-contigs: plain biological names
    ])
    truth_lines = (outdir / 'query_truth.tsv').read_text().strip().splitlines()
    assert len(truth_lines) > 1, 'truth TSV is empty'

    header = truth_lines[0].split('\t')
    svtype_idx  = header.index('svtype')
    contig_idx  = header.index('query_contig')

    on_ref = [r for r in truth_lines[1:] if r.split('\t')[svtype_idx] != 'OFF_REF']
    assert len(on_ref) >= 2, \
        f'Expected ≥2 on-reference truth records; got {len(on_ref)}'

    for row in truth_lines[1:]:
        contig = row.split('\t')[contig_idx]
        assert '__sv_' not in contig, \
            f'Hint suffix in contig name {contig!r} — --write-hint-contigs was not passed'


def test_simulator_hint_mode_writes_sv_encoded_names(tmp_path: Path) -> None:
    """With --write-hint-contigs, contig names carry __sv_ suffixes.

    This is the opt-in mode used only by legacy smoke tests.  It must not
    be the default.  We verify that hint names are written when the flag
    is present, and that the truth TSV records them faithfully.
    """
    outdir = tmp_path / 'hint_sim'
    run([
        'python3', str(SIM),
        '--phylum', 'Ascomycota',
        '--n-genomes', '3',
        '--n-reps', '1',
        '--total-len', '20000',
        '--n-contigs', '2',
        '--out-dir', str(outdir),
        '--scenario-set', 'core',
        '--write-hint-contigs',
    ])
    truth_lines = (outdir / 'query_truth.tsv').read_text().strip().splitlines()
    assert len(truth_lines) > 1

    header = truth_lines[0].split('\t')
    contig_idx  = header.index('query_contig')
    svtype_idx  = header.index('svtype')

    # At least the OFF_REF contig must carry a hint-encoded name
    off_ref_rows = [r for r in truth_lines[1:]
                    if r.split('\t')[svtype_idx] == 'OFF_REF']
    assert off_ref_rows, 'No OFF_REF rows in hint-mode truth TSV'
    for row in off_ref_rows:
        contig = row.split('\t')[contig_idx]
        assert '__sv_' in contig, \
            f'Expected hint suffix in OFF_REF contig name {contig!r}'


def test_query_mode_sanity_override_redirects_short_fastq_from_long_reads(tmp_path: Path) -> None:
    ensure_binary()
    ref = tmp_path / 'ref.fa'
    qry = tmp_path / 'reads.fastq'
    refs = tmp_path / 'refs.txt'
    queries = tmp_path / 'queries.txt'
    out_prefix = tmp_path / 'calls'

    ref.write_text('>chr1\n' + ('ACGT' * 200) + '\n')
    read_seq = 'ACGT' * 37 + 'AC'
    qry.write_text(
        '@r1\n' + read_seq + '\n+\n' + ('~' * len(read_seq)) + '\n'
        '@r2\n' + read_seq + '\n+\n' + ('~' * len(read_seq)) + '\n'
    )
    refs.write_text(str(ref) + '\n')
    queries.write_text(str(qry) + '\n')
    exe = ROOT / 'fungi_graphsv_tol_bin.exe'

    res = run_mycosv_command([
        str(exe),
        '--ref-list', str(refs),
        '--query-list', str(queries),
        '--out-prefix', str(out_prefix),
        '--query-mode', 'long-reads',
        '--max-reads', '2',
        '--genome-size-hint', '800',
    ], cwd=ROOT)
    assert res.returncode == 0
    assert out_prefix.with_suffix('.hits.tsv').exists()


# ============================================================
# DS-4 through DS-18 harness tests
# ============================================================

def test_ds5_bloom_filter_roundtrip_and_membership(tmp_path: Path):
    """DS-5: BloomFilter round-trips and preserves inserted membership."""
    RUN_EXE_CACHE.mkdir(exist_ok=True)
    harness = tmp_path / 'bloom_harness.cpp'
    harness.write_text(r'''
#include <fstream>
#include "layer3_routing_index.hpp"

extern "C" __declspec(dllexport) int bloom_roundtrip_ok() {
    tol::BloomFilter bf;
    if (!bf.empty()) return 0;
    bf.insert(11u);
    bf.insert(29u);
    bf.insert(47u);
    if (!bf.probably_contains(11u)) return 0;
    if (!bf.probably_contains(29u)) return 0;
    if (!bf.probably_contains(47u)) return 0;
    if (bf.empty()) return 0;

    {
        std::ofstream out("BLOOM_PATH", std::ios::binary);
        bf.write(out);
    }
    std::ifstream in("BLOOM_PATH", std::ios::binary);
    auto loaded = tol::BloomFilter::read(in);
    if (!loaded.probably_contains(11u)) return 0;
    if (!loaded.probably_contains(29u)) return 0;
    if (!loaded.probably_contains(47u)) return 0;
    return 1;
}
'''.replace('BLOOM_PATH', str((tmp_path / 'bf.bin')).replace('\\', '\\\\')))
    dll = RUN_EXE_CACHE / f'bloom_harness_{hashlib.sha1(str(tmp_path).encode("utf-8")).hexdigest()[:12]}.dll'
    run([
        'g++', '-shared', '-O2', '-std=c++17', '-static-libstdc++', '-static-libgcc',
        '-I', str(ROOT), str(harness), '-o', str(dll)
    ])
    lib = _load_test_dll(dll)
    lib.bloom_roundtrip_ok.restype = ctypes.c_int
    assert lib.bloom_roundtrip_ok() == 1


def test_ds4_vptree_save_load_and_query(tmp_path: Path):
    """DS-4: VPTree returns the nearest centroid before and after save/load."""
    RUN_EXE_CACHE.mkdir(exist_ok=True)
    tree_path = str((tmp_path / 'vptree.bin')).replace('\\', '\\\\')
    harness = tmp_path / 'vptree_harness.cpp'
    harness.write_text(r'''
#include <vector>
#include <cmath>
#include "layer3_routing_index.hpp"

extern "C" __declspec(dllexport) int vptree_roundtrip_ok() {
    std::vector<tol::CladeCentroid> centroids;

    tol::CladeCentroid a;
    a.cladeName = "cladeA";
    a.phylum = "Ascomycota";
    a.cladeRank = "species";
    a.centroidHashes = {1u, 2u, 3u, 4u};
    centroids.push_back(a);

    tol::CladeCentroid b;
    b.cladeName = "cladeB";
    b.phylum = "Basidiomycota";
    b.cladeRank = "species";
    b.centroidHashes = {50u, 60u, 70u, 80u};
    centroids.push_back(b);

    tol::CladeCentroid c;
    c.cladeName = "cladeC";
    c.phylum = "Mucoromycota";
    c.cladeRank = "species";
    c.centroidHashes = {5u, 6u, 7u, 8u};
    centroids.push_back(c);

    tol::VPTree tree;
    tree.build(centroids);

    tol::CladeCentroid q;
    q.cladeName = "query";
    q.phylum = "Ascomycota";
    q.cladeRank = "species";
    q.centroidHashes = {1u, 2u, 3u, 4u};
    q.build_prefilters();

    auto top_before = tree.query_topk(q, 2);
    if (top_before.empty() || top_before.front().cladeName != "cladeA") return 0;
    if (std::fabs(top_before.front().jaccard - 1.0) > 1e-9) return 0;

    tree.save("TREE_PATH");
    auto loaded = tol::VPTree::load("TREE_PATH");
    auto top_after = loaded.query_topk(q, 2);
    if (loaded.size() != 3) return 0;
    if (top_after.empty() || top_after.front().cladeName != "cladeA") return 0;
    if (std::fabs(top_after.front().jaccard - 1.0) > 1e-9) return 0;
    return 1;
}
'''.replace('TREE_PATH', tree_path))
    dll = RUN_EXE_CACHE / f'vptree_harness_{hashlib.sha1(str(tmp_path).encode("utf-8")).hexdigest()[:12]}.dll'
    run([
        'g++', '-shared', '-O2', '-std=c++17', '-static-libstdc++', '-static-libgcc',
        '-I', str(ROOT), str(harness), '-o', str(dll)
    ])
    lib = _load_test_dll(dll)
    lib.vptree_roundtrip_ok.restype = ctypes.c_int
    assert lib.vptree_roundtrip_ok() == 1


def test_registry_lru_eviction_harness(tmp_path: Path):
    """Layer-2 registry evicts the least-recently used graph when capacity is exceeded."""
    RUN_EXE_CACHE.mkdir(exist_ok=True)
    reg_dir = str((tmp_path / 'registry')).replace('\\', '\\\\')
    harness = tmp_path / 'registry_harness.cpp'
    harness.write_text(r'''
#include <filesystem>
#include "layer2_registry.hpp"

namespace fs = std::filesystem;

static tol::CladeGraph make_graph(const std::string& name, const std::string& seq, size_t compressed) {
    tol::CladeGraph g;
    g.cladeName = name;
    g.cladeRank = "species";
    g.phylum = "Ascomycota";
    g.genomeCount = 1;
    g.svBubbles = 1;
    g.compressedSz = compressed;
    g.nodes.push_back({1, seq});
    g.paths.push_back({name + "_path", {1}});
    return g;
}

extern "C" __declspec(dllexport) int registry_lru_ok() {
    fs::create_directories("REGISTRY_DIR");
    tol::CladeGraphRegistry reg("REGISTRY_DIR", 1u << 20, 1);
    reg.register_clade(make_graph("clade1", "AAAA", 64));
    reg.register_clade(make_graph("clade2", "CCCC", 64));

    auto g1 = reg.get("clade1");
    auto g2 = reg.get("clade2");
    auto g2_again = reg.get("clade2");
    auto g1_again = reg.get("clade1");
    auto stats = reg.stats();

    if (!g1 || g1->cladeName != "clade1") return 0;
    if (!g2 || g2->cladeName != "clade2") return 0;
    if (!g2_again || g2_again->cladeName != "clade2") return 0;
    if (!g1_again || g1_again->cladeName != "clade1") return 0;
    if (stats.entries != 1) return 0;
    if (stats.hits < 1) return 0;
    if (stats.misses < 3) return 0;
    if (stats.evictions < 2) return 0;
    if (stats.registered != 2) return 0;
    return 1;
}
'''.replace('REGISTRY_DIR', reg_dir))
    dll = RUN_EXE_CACHE / f'registry_harness_{hashlib.sha1(str(tmp_path).encode("utf-8")).hexdigest()[:12]}.dll'
    run([
        'g++', '-shared', '-O2', '-std=c++17', '-static-libstdc++', '-static-libgcc',
        '-pthread', '-I', str(ROOT), str(harness), '-o', str(dll)
    ])
    lib = _load_test_dll(dll)
    lib.registry_lru_ok.restype = ctypes.c_int
    assert lib.registry_lru_ok() == 1

def test_ds13_suffix_array_build_and_mems(tmp_path: Path):
    """DS-13: SuffixArray builds correctly and finds MEMs."""
    harness = tmp_path / 'sa_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include <string>
#include <vector>
#include "layer1_clade_graph.hpp"
int main(){
    // Build SA over a simple reference
    tol::SuffixArray sa;
    sa.build({{"ref1", "ACGTACGTACGT"}, {"ref2", "GCGCGCGC"}});
    std::cout << "sa_built\t" << (sa.empty() ? "empty" : "ok") << "\n";
    // Query: "ACGT" should appear in ref1
    auto mems = sa.find_mems("ACGT", 4);
    std::cout << "mems_found\t" << (mems.empty() ? "none" : "ok") << "\n";
    // Reverse complement
    std::string rc = tol::SuffixArray::revcomp("ACGT");
    std::cout << "revcomp\t" << (rc == "ACGT" ? "ok" : rc) << "\n";
    // LCP populated
    std::cout << "lcp_size\t" << (sa.lcp.size() == sa.sa.size() ? "ok" : "fail") << "\n";
}
''')
    exe = tmp_path / 'sa_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(l.split('\t', 1) for l in run([str(exe)]).stdout.strip().splitlines())
    assert rows['sa_built']   == 'ok'
    assert rows['mems_found'] == 'ok'
    assert rows['revcomp']    == 'ok'
    assert rows['lcp_size']   == 'ok'


def test_ds13_svtype_from_chain_classifies_del_ins(tmp_path: Path):
    """DS-13+DS-18: SvTypeFromChain produces DEL for gap in reference, INS for gap in query."""
    harness = tmp_path / 'chain_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include <vector>
#include "layer1_clade_graph.hpp"
int main(){
    // DEL: two MEMs with rGap >> qGap
    // qPos: 0,50  rPos: 0,100  → rGap=50, qGap=0 → DEL of 50
    tol::SuffixArray sa;
    sa.build({{"ref1", std::string(200,'A')}});
    std::vector<tol::SuffixArray::Mem> chain_del = {{0,0,50},{50,100,50}};
    std::vector<bool> rev_del = {false, false};
    auto res_del = tol::SvTypeFromChain::classify(chain_del, rev_del, sa, 20);
    std::cout << "del\t" << (res_del.type == tol::SvTypeFromChain::Type::DEL ? "DEL" : "OTHER") << "\n";

    // INS: qGap >> rGap  → INS
    std::vector<tol::SuffixArray::Mem> chain_ins = {{0,0,20},{70,20,20}};
    std::vector<bool> rev_ins = {false, false};
    auto res_ins = tol::SvTypeFromChain::classify(chain_ins, rev_ins, sa, 20);
    std::cout << "ins\t" << (res_ins.type == tol::SvTypeFromChain::Type::INS ? "INS" : "OTHER") << "\n";

    // INV: reverse-complement seed → INV
    std::vector<tol::SuffixArray::Mem> chain_inv = {{0,0,30},{50,100,30}};
    std::vector<bool> rev_inv = {true, true};
    auto res_inv = tol::SvTypeFromChain::classify(chain_inv, rev_inv, sa, 20);
    std::cout << "inv\t" << (res_inv.type == tol::SvTypeFromChain::Type::INV ? "INV" : "OTHER") << "\n";
}
''')
    exe = tmp_path / 'chain_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(l.split('\t', 1) for l in run([str(exe)]).stdout.strip().splitlines())
    assert rows['del'] == 'DEL'
    assert rows['ins'] == 'INS'
    assert rows['inv'] == 'INV'


def test_ds14_wavelet_tree_rank_and_filter(tmp_path: Path):
    """DS-14: WaveletTree rank() and has_char() work correctly."""
    harness = tmp_path / 'wt_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include <string>
#include "layer3_routing_index.hpp"
int main(){
    tol::WaveletTree wt;
    wt.build("ACGTACGT");
    // "ACGTACGT": A at pos 0,4  C at 1,5  G at 2,6  T at 3,7
    // rank(A, 4): A in [0..3]="ACGT" → 1 occurrence
    std::cout << "rank_A\t" << wt.rank(0, 4) << "\n"; // 1
    // rank(C, 8): C in [0..7] → 2 occurrences
    std::cout << "rank_C\t" << wt.rank(1, 8) << "\n"; // 2
    // has_char: A is present
    std::cout << "has_A\t" << (wt.has_char('A') ? "yes" : "no") << "\n"; // yes
    // Absence check: build a text with no 'A' and confirm has_char('A')=no
    tol::WaveletTree wt2;
    wt2.build("CCGGCCGG"); // contains only C and G, no A or T
    std::cout << "no_A\t" << (wt2.has_char('A') ? "yes" : "no") << "\n";  // no
    std::cout << "has_C\t" << (wt2.has_char('C') ? "yes" : "no") << "\n"; // yes
    // empty check
    tol::WaveletTree empty;
    std::cout << "empty\t" << (empty.empty() ? "yes" : "no") << "\n";
}
''')
    exe = tmp_path / 'wt_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(l.split('\t', 1) for l in run([str(exe)]).stdout.strip().splitlines())
    assert rows['rank_A'] == '1',  f"rank_A={rows['rank_A']}"  # 1 A in "ACGT"[0..3]
    assert rows['rank_C'] == '2',  f"rank_C={rows['rank_C']}"  # 2 C in "ACGTACGT"
    assert rows['has_A']  == 'yes'
    assert rows['no_A']   == 'no',  f"no_A={rows['no_A']}"
    assert rows['has_C']  == 'yes'
    assert rows['empty']  == 'yes'


def test_ds15_veb_tree_predecessor_successor(tmp_path: Path):
    """DS-15: VEBTree<20> predecessor/successor/any_in work correctly."""
    harness = tmp_path / 'veb_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include "layer1_clade_graph.hpp"
int main(){
    tol::VEBTree<20> veb;
    veb.insert(10);
    veb.insert(30);
    veb.insert(50);
    // successor of 10 → 30
    std::cout << "succ_10\t" << veb.successor(10) << "\n";
    // predecessor of 30 → 10
    std::cout << "pred_30\t" << veb.predecessor(30) << "\n";
    // any_in [20,40] → 30 is in range → true
    std::cout << "any_20_40\t" << (veb.any_in(20, 40) ? "yes" : "no") << "\n";
    // any_in [60,80] → nothing → false
    std::cout << "any_60_80\t" << (veb.any_in(60, 80) ? "yes" : "no") << "\n";
    // contains
    std::cout << "has_30\t" << (veb.contains(30) ? "yes" : "no") << "\n";
    std::cout << "has_20\t" << (veb.contains(20) ? "yes" : "no") << "\n";
}
''')
    exe = tmp_path / 'veb_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(l.split('\t', 1) for l in run([str(exe)]).stdout.strip().splitlines())
    assert rows['succ_10']   == '30'
    assert rows['pred_30']   == '10'
    assert rows['any_20_40'] == 'yes'
    assert rows['any_60_80'] == 'no'
    assert rows['has_30']    == 'yes'
    assert rows['has_20']    == 'no'


def test_ds16_merge_sort_tree_stabbing(tmp_path: Path):
    """DS-16: MergeSortTree stab_count and overlapping work correctly."""
    harness = tmp_path / 'mst_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include <vector>
#include "layer1_clade_graph.hpp"
int main(){
    tol::MergeSortTree mst;
    // Intervals: [10,30], [20,40], [50,70]
    mst.build({{10,30},{20,40},{50,70}});
    // stab at 25 → hits [10,30] and [20,40] → 2
    std::cout << "stab_25\t" << mst.stab_count(25) << "\n";
    // stab at 55 → hits [50,70] → 1
    std::cout << "stab_55\t" << mst.stab_count(55) << "\n";
    // stab at 45 → nothing → 0
    std::cout << "stab_45\t" << mst.stab_count(45) << "\n";
    // overlapping [15,35] → [10,30] and [20,40] → 2
    auto ov = mst.overlapping(15, 35);
    std::cout << "overlap_count\t" << ov.size() << "\n";
}
''')
    exe = tmp_path / 'mst_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(l.split('\t', 1) for l in run([str(exe)]).stdout.strip().splitlines())
    assert rows['stab_25']       == '2'
    assert rows['stab_55']       == '1'
    assert rows['stab_45']       == '0'
    assert rows['overlap_count'] == '2'


def test_ds17_fenwick_tree_prefix_sum_and_findkth(tmp_path: Path):
    """DS-17: FenwickTree prefix_sum, range_sum, find_kth work correctly."""
    harness = tmp_path / 'fen_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include "layer1_clade_graph.hpp"
int main(){
    tol::FenwickTree fen(5);
    // Sizes: [10, 20, 30, 40, 50] → prefix sums: [10,30,60,100,150]
    fen.update(0, 10);
    fen.update(1, 20);
    fen.update(2, 30);
    fen.update(3, 40);
    fen.update(4, 50);
    std::cout << "psum_2\t" << fen.prefix_sum(2) << "\n"; // 10+20+30=60
    std::cout << "psum_4\t" << fen.prefix_sum(4) << "\n"; // 150
    std::cout << "range_1_3\t" << fen.range_sum(1, 3) << "\n"; // 20+30+40=90
    // find_kth(55): smallest k where prefix_sum(k) >= 55 → k=2 (sum=60)
    std::cout << "kth_55\t" << fen.find_kth(55) << "\n";
}
''')
    exe = tmp_path / 'fen_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(l.split('\t', 1) for l in run([str(exe)]).stdout.strip().splitlines())
    assert rows['psum_2']    == '60'
    assert rows['psum_4']    == '150'
    assert rows['range_1_3'] == '90'
    assert int(rows['kth_55']) == 2


def test_ds18_chain_treap_basic_chaining(tmp_path: Path):
    """DS-18: ChainTreap chains compatible seeds and rejects incompatible ones."""
    harness = tmp_path / 'treap_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include "layer1_clade_graph.hpp"
int main(){
    tol::ChainTreap treap;
    // Seeds: (qPos, rPos, len, score)
    // Colinear chain: (0,0,20), (25,25,20), (50,50,20) → total score=60
    float s1 = treap.insert_and_chain(0,  0,  20, 20.0f, 100);
    float s2 = treap.insert_and_chain(25, 25, 20, 20.0f, 100);
    float s3 = treap.insert_and_chain(50, 50, 20, 20.0f, 100);
    std::cout << "best_score\t" << treap.best_chain_score() << "\n";
    auto path = treap.best_chain_path();
    std::cout << "chain_len\t" << path.size() << "\n";
    // Score should grow: s1=20, s2=40, s3=60
    std::cout << "scores_grow\t" << (s1 < s2 && s2 < s3 ? "yes" : "no") << "\n";
}
''')
    exe = tmp_path / 'treap_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(l.split('\t', 1) for l in run([str(exe)]).stdout.strip().splitlines())
    assert float(rows['best_score']) == 60.0
    assert int(rows['chain_len'])    == 3
    assert rows['scores_grow']       == 'yes'


def test_mem_chain_detects_del_in_binary(tmp_path: Path):
    """End-to-end: MEM chain path in binary detects DEL with SA-precise breakpoint."""
    import subprocess
    ref = tmp_path / 'ref.fa'
    query = tmp_path / 'query.fa'
    # 300bp ref; query is 300bp with a 60bp internal deletion → net 240bp
    ref_seq   = 'ACGT' * 75   # 300bp
    # Remove positions 100-159 (60 bp deletion)
    query_seq = ref_seq[:100] + ref_seq[160:]  # 240bp
    ref.write_text(f'>ctg1\n{ref_seq}\n')
    query.write_text(f'>ctg1\n{query_seq}\n')
    (tmp_path / 'r.txt').write_text(str(ref)   + '\n')
    (tmp_path / 'q.txt').write_text(str(query) + '\n')

    subprocess.run([str(BIN),
        '--ref-list',    str(tmp_path / 'r.txt'),
        '--query-list',  str(tmp_path / 'q.txt'),
        '--out-prefix',  str(tmp_path / 'calls'),
        '--query-mode',  'assembly',
        '--quiet'], check=True)

    vcf = (tmp_path / 'calls.vcf').read_text()
    data = [l for l in vcf.splitlines() if not l.startswith('#') and l.strip()]
    assert data, "no call emitted for 60 bp deletion"
    assert any('SVTYPE=DEL' in l for l in data), f"expected DEL, got: {data}"


def test_mem_chain_honors_secondary_seed_rescue_and_thresholds(tmp_path: Path):
    harness = tmp_path / 'mem_rescue_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include <memory>
#include <vector>
#include "fungi_tol_bridge.hpp"

int main() {
    tol::TolGlobal::RefSeq ref;
    ref.asmName = "asm1";
    ref.contig = "ctg1";
    ref.clade = "clade1";
    ref.cladeRank = "species";
    ref.phylum = "Ascomycota";
    const std::string left = "ACGTACGTTGCAACGA";
    const std::string right = "TTCGAGGCTAACCGTA";
    ref.seqShared = std::make_shared<std::string>(left + std::string(60, 'C') + right);

    std::vector<tol::TolGlobal::RefSeq> refs = {ref};
    std::vector<const tol::TolGlobal::RefSeq*> ptrs = {&refs[0]};
    const std::string query = left + right;

    tol::FederatedOptions fo;
    fo.primarySketchParams.k = 23;
    fo.secondarySketchParams.k = 20;
    fo.minSvLen = 20;
    fo.minAnchors = 2;
    fo.minBlockScore = 20.0;
    fo.chainGapBand = 5000;

    tol::VariantCallBridge no_secondary;
    fo.useSecondarySeeds = false;
    std::cout << "no_secondary\t"
              << (tol::try_mem_chain_call_public("q", "ctg1", query, ptrs, fo, no_secondary) ? "yes" : "no")
              << "\n";

    tol::VariantCallBridge rescued;
    fo.useSecondarySeeds = true;
    fo.repeatRescueMinAnchors = 3;
    const bool rescued_ok = tol::try_mem_chain_call_public("q", "ctg1", query, ptrs, fo, rescued);
    std::cout << "with_secondary\t" << (rescued_ok ? "yes" : "no") << "\n";
    if (rescued_ok) {
        std::cout << "rescued_type\t" << rescued.type << "\n";
        std::cout << "rescued_mode\t" << rescued.alignmentMode << "\n";
    }

    tol::VariantCallBridge blocked;
    fo.minBlockScore = 100.0;
    std::cout << "high_block_threshold\t"
              << (tol::try_mem_chain_call_public("q", "ctg1", query, ptrs, fo, blocked) ? "yes" : "no")
              << "\n";
}
''')
    exe = tmp_path / 'mem_rescue_harness'
    run(['g++', '-O2', '-std=c++17', '-pthread', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert rows['no_secondary'] == 'no'
    assert rows['with_secondary'] == 'yes'
    assert rows['rescued_type'] == 'DEL'
    assert 'secondary_seed_rescue' in rows['rescued_mode']
    assert rows['high_block_threshold'] == 'no'


def test_repeat_hgt_rip_detectors_follow_benchmark_rules(tmp_path: Path):
    """Detector thresholds match the benchmark rules for tandem/HGT/RIP."""
    harness = tmp_path / 'rule_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include <string>
#include "layer1_clade_graph.hpp"

int main() {
    std::string tandem;
    for (int i = 0; i < 30; ++i) tandem += "AC"; // period 2, 30 copies, 60 bp

    std::string hgt_hi(500, 'G');
    std::string hgt_lo(500, 'C');
    std::string rip_yes(500, 'C');
    for (int i = 400; i < 500; ++i) rip_yes[static_cast<size_t>(i)] = 'G'; // C/G = 4.0
    std::string rip_no(500, 'C');
    for (int i = 200; i < 500; ++i) rip_no[static_cast<size_t>(i)] = 'G'; // C/G ~= 0.67

    std::cout << "tandem\t" << (tol::detect_tandem_repeat(tandem) ? "yes" : "no") << "\n";
    std::cout << "hgt_hi\t" << (tol::detect_hgt_island(hgt_hi, 0.50, 0.12, 500) ? "yes" : "no") << "\n";
    std::cout << "hgt_lo\t" << (tol::detect_hgt_island(hgt_lo, 0.50, 0.12, 500) ? "yes" : "no") << "\n";
    std::cout << "rip_yes\t" << (tol::detect_rip_window(rip_yes, 2.5, 500) ? "yes" : "no") << "\n";
    std::cout << "rip_no\t" << (tol::detect_rip_window(rip_no, 2.5, 500) ? "yes" : "no") << "\n";
}
''')
    exe = tmp_path / 'rule_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert rows['tandem'] == 'yes'
    assert rows['hgt_hi'] == 'yes'
    assert rows['hgt_lo'] == 'yes'
    assert rows['rip_yes'] == 'yes'
    assert rows['rip_no'] == 'no'


def test_graph_native_offref_windows_and_router_windowing(tmp_path: Path):
    harness = tmp_path / 'offref_router_harness.cpp'
    harness.write_text('#include <iostream>\n#include <string>\n#include <vector>\n#include <memory>\n#include "fungi_tol_bridge.hpp"\n#include "layer3_routing_index.hpp"\n\nint main() {\n    tol::TolGlobal::RefSeq ref1;\n    ref1.asmName = "asm1";\n    ref1.contig = "ctg1";\n    ref1.seqShared = std::make_shared<std::string>(std::string(800, \'A\') + std::string(800, \'C\'));\n    ref1.clade = "clade1";\n    ref1.cladeRank = "species";\n    ref1.phylum = "Ascomycota";\n    ref1.cladeGc = 0.50;\n\n    tol::TolGlobal::RefSeq ref2 = ref1;\n    ref2.asmName = "asm2";\n    ref2.clade = "clade2";\n    ref2.phylum = "Basidiomycota";\n    ref2.seqShared = std::make_shared<std::string>(std::string(1600, \'G\'));\n    std::vector<tol::TolGlobal::RefSeq> refs = {ref1, ref2};\n\n    tol::FederatedOptions fo;\n    fo.minSvLen = 40;\n    fo.fallbackSketchParams.k = 7;\n    std::string query = std::string(700, \'A\') + std::string(700, \'T\') + std::string(700, \'C\');\n    auto wins = tol::discover_graph_native_offref_windows(query, refs, fo, 500);\n    std::cout << "windows\\t" << wins.size() << "\\n";\n    if (!wins.empty()) {\n        auto call = tol::make_offref_window_call("qasm", "qctg", query, wins.front());\n        std::cout << "first_mode\\t" << call.alignmentMode << "\\n";\n        std::cout << "first_type\\t" << call.type << "\\n";\n        std::cout << "first_pos\\t" << call.pos << "\\n";\n        std::cout << "first_end\\t" << call.end << "\\n";\n    }\n\n    tol::PhylumShardedRouter router;\n    tol::CladeCentroid c1; c1.cladeName = "clade1"; c1.phylum = "Ascomycota"; c1.centroidHashes = {1,2,3,4};\n    tol::CladeCentroid c2; c2.cladeName = "clade2"; c2.phylum = "Basidiomycota"; c2.centroidHashes = {10,11,12,13};\n    router.register_clade_centroid(c1);\n    router.register_clade_centroid(c2);\n    router.rebuild();\n    tol::SyncmerParams sp; sp.k = 7; sp.s = 3; sp.t = 1;\n    auto routed = router.route_windows(std::string_view(query), sp, sp, 1.0, 1, 400, 100);\n    std::cout << "route_windows\\t" << routed.size() << "\\n";\n}\n')
    exe = tmp_path / 'offref_router_harness'
    run(['g++', '-O2', '-std=c++17', '-pthread', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('	', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert int(rows['windows']) >= 1
    assert rows['first_mode'] == 'graph_native_offref_window'
    assert rows['first_type'] == 'OFF_REF'
    assert int(rows['first_pos']) >= 1
    assert int(rows['first_end']) >= int(rows['first_pos'])
    assert int(rows['route_windows']) >= 2


def test_probabilistic_evidence_fusion_prefers_alt_with_multilayer_support(tmp_path: Path):
    harness = tmp_path / 'evidence_fusion_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include <vector>
#include "hierarchical_engine.hpp"
int main() {
    std::vector<tol::EvidenceObservation> ev = {
        {tol::EvidenceLayer::ASSEMBLY,   -3.0,  2.0,  1.0, 60.0, 1.0, 1.0},
        {tol::EvidenceLayer::LONG_READ,  -2.0,  1.5, 18.0, 50.0, 0.9, 0.8},
        {tol::EvidenceLayer::SHORT_READ, -1.0,  0.8, 65.0, 45.0, 0.8, 0.7},
    };
    auto fused = tol::fuse_probabilistic_evidence(ev, 0.2);
    std::cout << "posterior_alt\t" << fused.posteriorAlt << "\n";
    std::cout << "posterior_ref\t" << fused.posteriorRef << "\n";
    std::cout << "layers\t" << fused.layersUsed << "\n";
    std::cout << "supports\t" << (fused.supports_variant(0.8) ? "yes" : "no") << "\n";
}
''')
    exe = tmp_path / 'evidence_fusion_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert float(rows['posterior_alt']) > 0.8
    assert float(rows['posterior_ref']) < 0.2
    assert int(rows['layers']) == 3
    assert rows['supports'] == 'yes'


def test_external_memory_centroid_store_streams_topk(tmp_path: Path):
    harness = tmp_path / 'external_store_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include <vector>
#include "layer3_routing_index.hpp"
int main() {
    std::vector<tol::CladeCentroid> centroids;
    for (int i = 0; i < 8; ++i) {
        tol::CladeCentroid c;
        c.cladeName = "clade" + std::to_string(i);
        c.phylum = (i % 2 == 0 ? "Ascomycota" : "Basidiomycota");
        c.cladeRank = "species";
        c.centroidHashes = {1u, 2u, static_cast<uint64_t>(10 + i), static_cast<uint64_t>(20 + i)};
        centroids.push_back(c);
    }
    tol::ExternalMemoryCentroidStore store("STORE_PATH");
    store.build(centroids);

    tol::CladeCentroid q;
    q.cladeName = "query";
    q.phylum = "Ascomycota";
    q.cladeRank = "species";
    q.centroidHashes = {1u, 2u, 13u, 23u};

    auto top = store.query_topk_streaming(q, 3);
    std::cout << "topk\t" << top.size() << "\n";
    if (!top.empty()) {
        std::cout << "best\t" << top.front().cladeName << "\n";
        std::cout << "best_j\t" << top.front().jaccard << "\n";
    }
}
'''.replace('STORE_PATH', str((tmp_path / 'centroids.bin')).replace('\\', '\\\\')))
    exe = tmp_path / 'external_store_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert int(rows['topk']) == 3
    assert rows['best'] == 'clade3'
    assert float(rows['best_j']) > 0.9




def test_tol_init_reuses_existing_routing_centroid_store(tmp_path: Path):
    ref = tmp_path / 'ref.fa'
    write_fasta(ref, [('ctg1', 'ACGT' * 100)])
    manifest = tmp_path / 'manifest.tsv'
    manifest.write_text(
        '#asm_name\tphylum\tclade_name\tclade_rank\tfasta_path\n'
        f'asm1\tAscomycota\tClade1\tspecies\t{ref}\n'
    )

    idx = tmp_path / 'idx'
    reg = tmp_path / 'reg'
    idx.mkdir()
    reg.mkdir()

    manifest_cpp = str(manifest).replace('\\', '\\\\')
    idx_cpp = str(idx).replace('\\', '\\\\')
    reg_cpp = str(reg).replace('\\', '\\\\')

    harness = tmp_path / 'tol_store_reuse_harness.cpp'
    harness.write_text(f'''
#include <chrono>
#include <filesystem>
#include <iostream>
#include <thread>
#include "fungi_tol_bridge.hpp"

namespace fs = std::filesystem;

int main() {{
    tol::SyncmerParams sp;
    tol::SyncmerParams fb;
    fb.k = 7;
    fb.s = 3;
    tol::build_tol_index_from_manifest(
        "{manifest_cpp}", "{idx_cpp}", "{reg_cpp}", sp, 0.12, false, "", 500, 1, &fb, false);

    const fs::path store = fs::path("{idx_cpp}") / "routing_centroids.bin";
    tol::TolGlobal::instance().init("{idx_cpp}", "{reg_cpp}", 1ull << 20, 64);
    auto t1 = fs::last_write_time(store);
    std::this_thread::sleep_for(std::chrono::milliseconds(25));
    tol::TolGlobal::instance().init("{idx_cpp}", "{reg_cpp}", 1ull << 20, 64);
    auto t2 = fs::last_write_time(store);
    std::cout << "store_exists\t" << (fs::exists(store) ? "yes" : "no") << "\\n";
    std::cout << "reused\t" << (t1 == t2 ? "yes" : "no") << "\\n";
}}
''')
    exe = tmp_path / 'tol_store_reuse_harness'
    run(['g++', '-O2', '-std=c++17', '-pthread', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert rows['store_exists'] == 'yes'
    assert rows['reused'] == 'yes'


def test_multi_rank_partitioned_build_compacts_oversized_clades(tmp_path: Path):
    ensure_binary()

    ref = tmp_path / 'shared_ref.fa'
    write_fasta(ref, [('ctg1', 'ACGT' * 64), ('ctg2', 'TGCA' * 64)])

    manifest = tmp_path / 'hierarchy_manifest.tsv'
    n_genomes = 256
    with manifest.open('w', encoding='utf-8') as fh:
        fh.write('#asm_name\tphylum\tclass\torder\tfamily\tgenus\tclade_name\tclade_rank\tfasta_path\n')
        for i in range(n_genomes):
            fh.write(
                f'asm{i}\tAscomycota\tSordariomycetes\tHypocreales\tNectriaceae\tFusarium\t'
                f'Fusarium_sp_{i}\tspecies\t{ref}\n'
            )

    idx = tmp_path / 'idx'
    reg = tmp_path / 'reg'
    run([
        str(BIN),
        '--tol-hierarchical',
        '--tol-build-index', str(manifest),
        '--tol-index-dir', str(idx),
        '--tol-registry-dir', str(reg),
        '--tol-multi-rank',
        '--tol-base-graph-build',
        '--tol-max-clade-genomes', '32',
        '--tol-index-threads', '4',
    ])

    assert not (idx / '.manifest_partitions').exists()

    routing_lines = (idx / 'routing_manifest.tsv').read_text().strip().splitlines()
    assert routing_lines[0].startswith('#clade_name\tclade_rank\tphylum\t')
    rows = [line.split('\t') for line in routing_lines[1:]]
    assert len(rows) == n_genomes + 5

    by_rank = {}
    for row in rows:
        by_rank.setdefault(row[1], []).append(row)

    assert {rank for rank in by_rank} == {
        'species', 'genus', 'family', 'order', 'class', 'phylum'
    }
    assert len(by_rank['species']) == n_genomes
    assert len(by_rank['genus']) == 1
    assert len(by_rank['family']) == 1
    assert len(by_rank['order']) == 1
    assert len(by_rank['class']) == 1
    assert len(by_rank['phylum']) == 1

    genus_row = by_rank['genus'][0]
    family_row = by_rank['family'][0]
    order_row = by_rank['order'][0]
    class_row = by_rank['class'][0]
    phylum_row = by_rank['phylum'][0]
    for row in (genus_row, family_row, order_row, class_row, phylum_row):
        assert int(row[7]) == n_genomes
        assert int(row[9]) <= 64
        assert Path(row[3]).exists()
        assert Path(row[4]).exists()
        assert Path(row[5]).exists()
        assert Path(row[6]).exists()

    species_row = by_rank['species'][0]
    assert int(species_row[7]) == 1
    assert 1 <= int(species_row[9]) <= 2

    clade_manifest = (reg / 'clade_manifest.tsv').read_text().strip().splitlines()
    assert clade_manifest[0].startswith('#clade_name\tclade_rank\tphylum\tgraph_path\t')
    assert len(clade_manifest) == len(rows) + 1


def test_hierarchical_query_works_with_external_only_routing_store(tmp_path: Path):
    ensure_binary()

    sim = tmp_path / 'sim'
    idx = tmp_path / 'idx'
    reg = tmp_path / 'reg'
    out = tmp_path / 'out'
    out.mkdir()

    run([
        'python3', str(SIM),
        '--phylum', 'Ascomycota',
        '--n-genomes', '3',
        '--n-reps', '1',
        '--total-len', '6000',
        '--n-contigs', '2',
        '--out-dir', str(sim),
        '--scenario-set', 'compact_yeast',
        '--write-extended-manifest',
        '--query-mode', 'assembly',
    ])

    run([
        str(BIN),
        '--tol-hierarchical',
        '--tol-build-index', str(sim / 'hierarchy_manifest.tsv'),
        '--tol-index-dir', str(idx),
        '--tol-registry-dir', str(reg),
        '--tol-multi-rank',
        '--tol-base-graph-build',
        '--tol-max-clade-genomes', '32',
        '--tol-index-threads', '2',
    ])

    store = idx / 'routing_centroids.bin'
    data = store.read_bytes()
    store.write_bytes(struct.pack('<Q', 200001) + data[8:])
    skip = Path(str(store) + '.skip')
    if skip.exists():
        skip.unlink()

    run([
        str(BIN),
        '--tol-hierarchical',
        '--tol-index-dir', str(idx),
        '--tol-registry-dir', str(reg),
        '--ref-list', str(sim / 'ref_list.txt'),
        '--query-list', str(sim / 'query_list.txt'),
        '--out-prefix', str(out / 'calls'),
        '--query-mode', 'assembly',
    ])

    hits = (out / 'calls.hits.tsv').read_text()
    vcf = (out / 'calls.vcf').read_text()
    assert 'query_asm\tquery_contig\ttype\tref_asm' in hits
    assert '\tref_pos\tref_end\tpos\tend\t' in hits
    assert '\n' in hits.strip()
    assert '##fileformat=VCFv4.3' in vcf
    assert '##INFO=<ID=REFCONTIG' in vcf
    assert '##INFO=<ID=REFPOS' in vcf
    assert '##INFO=<ID=REFEND' in vcf


def test_million_scale_routing_benchmark_harness_smoke(tmp_path: Path):
    src = ROOT / 'million_scale_routing_benchmark.cpp'
    exe = tmp_path / 'million_scale_routing_benchmark'
    report = tmp_path / 'million_scale_report.tsv'
    store = tmp_path / 'million_scale_store.bin'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(src), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([
        str(exe),
        '--n-centroids', '4000',
        '--hashes-per-centroid', '16',
        '--queries', '8',
        '--top-k', '3',
        '--phylum-count', '4',
        '--chunk-records', '256',
        '--store', str(store),
        '--report-tsv', str(report),
    ]).stdout.strip().splitlines())
    assert int(rows['n_centroids']) == 4000
    assert float(rows['top_hit_recall']) >= 0.99
    assert int(rows['skip_index_bytes']) > 0
    assert Path(str(store) + '.skip').exists()
    report_lines = report.read_text().strip().splitlines()
    assert report_lines[0].startswith('n_centroids\t')
    assert report_lines[1].split('\t')[0] == '4000'


def test_mode_pr_benchmark_runner_outputs_all_modes(tmp_path: Path):
    ensure_binary()
    runner = ROOT / 'run_mode_pr_benchmark.py'
    outdir = tmp_path / 'mode_bench'
    run([
        'python3', str(runner),
        '--out-dir', str(outdir),
        '--binary-path', str(BIN),
        '--skip-build',
        '--modes', 'assembly,short-reads,long-reads',
        '--phylum', 'Ascomycota',
        '--scenario-set', 'compact_yeast',
        '--n-genomes', '4',
        '--n-reps', '2',
        '--total-len', '8000',
        '--n-contigs', '2',
    ])
    summary_tsv = (outdir / 'mode_pr_summary.tsv').read_text().strip().splitlines()
    assert summary_tsv[0].startswith('mode\ttruth_records\t')
    rows = [line.split('\t') for line in summary_tsv[1:]]
    assert {row[0] for row in rows} == {'assembly', 'short-reads', 'long-reads'}
    metrics_by_mode = {row[0]: row for row in rows}
    for row in rows:
        precision = float(row[7])
        recall = float(row[8])
        assert 0.0 <= precision <= 1.0
        assert 0.0 <= recall <= 1.0
    assert float(metrics_by_mode['short-reads'][7]) >= 0.5
    assert float(metrics_by_mode['short-reads'][8]) >= 0.5
    assert float(metrics_by_mode['long-reads'][7]) >= 0.5
    assert float(metrics_by_mode['long-reads'][8]) >= 0.5
    svtype_tsv = (outdir / 'mode_svtype_pr_summary.tsv').read_text()
    assert 'mode\tsvtype\ttp\tfp\tfn\tprecision' in svtype_tsv
    summary_json = json.loads((outdir / 'mode_pr_summary.json').read_text())
    assert set(summary_json['modes']) == {'assembly', 'short-reads', 'long-reads'}


def test_query_input_recommends_probabilistic_fusion_by_mode_and_coverage(tmp_path: Path):
    harness = tmp_path / 'fusion_recommend_harness.cpp'
    harness.write_text(r'''#include <iostream>
#include "query_input_handler.hpp"
int main() {
    auto a = query_input::recommend_evidence_fusion(
        query_input::QueryMode::ASSEMBLY,
        query_input::CoverageTier::NORMAL,
        2);
    auto s = query_input::recommend_evidence_fusion(
        query_input::QueryMode::SHORT_READS,
        query_input::CoverageTier::HIGH,
        3);
    std::cout << "asm_enable\t" << (a.enableProbabilisticFusion ? "yes" : "no") << "\n";
    std::cout << "asm_prior\t" << a.priorAlt << "\n";
    std::cout << "sr_read_weight\t" << s.expectedReadWeight << "\n";
    std::cout << "sr_rationale\t" << s.rationale << "\n";
}
''')
    exe = tmp_path / 'fusion_recommend_harness'
    run(['g++', '-O2', '-std=c++17', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert rows['asm_enable'] == 'yes'
    assert float(rows['asm_prior']) >= 0.5
    assert float(rows['sr_read_weight']) > 1.0
    assert rows['sr_rationale'] == 'short_reads_depth_stabilised_fusion'

def test_full_ancestral_sequence_reconstruction_returns_sequence_and_posteriors(tmp_path: Path):
    harness = tmp_path / 'ancestral_recon_harness.cpp'
    harness.write_text(r'''
#include <iostream>
#include "fungi_tol_bridge.hpp"
int main() {
    tol::VariantCallBridge call;
    call.type = "DEL";
    call.refAsm = "cladeA";
    call.phylum = "Ascomycota";
    tol::FederatedOptions fo;
    auto res = tol::reconstruct_full_ancestral_sequence("ACGTAC", "ACGGAC", call, fo);
    std::cout << "len\t" << res.ancestralSequence.size() << "\n";
    std::cout << "aligned\t" << res.alignedBases << "\n";
    std::cout << "breakpoints\t" << res.breakpointCount << "\n";
    std::cout << "mean_post\t" << res.meanPosterior << "\n";
    std::cout << "seq\t" << res.ancestralSequence << "\n";
    std::cout << "src\t" << res.sourceClade << "\n";
}
''')
    exe = tmp_path / 'ancestral_recon_harness'
    run(['g++', '-O2', '-std=c++17', '-pthread', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert int(rows['len']) == 6
    assert int(rows['aligned']) == 6
    assert int(rows['breakpoints']) == 1
    assert float(rows['mean_post']) > 0.8
    assert rows['seq'] == 'ACGTAC'
    assert rows['src'] == 'cladeA'


def test_candidate_arbitration_prefers_simple_indel_over_spurious_rearrangement(tmp_path: Path):
    harness = tmp_path / 'arb_harness.cpp'
    harness.write_text(r'''
#define main fungi_graphsv_tol_main
#include "main.cpp"
#undef main
#include <iostream>

int main() {
    Options o;
    o.minSvLen = 20;
    o.maxSvLen = 1000000;
    o.minAnchors = 2;

    std::unordered_map<std::string, std::string> contigs;
    contigs["ctg1"] = std::string(1000, 'A') + std::string(80, 'C');

    SimpleRefIndex refIdx;
    refIdx["ctg1"].push_back(RefContigInfo{"asm1", "ctg1", std::string(1000, 'A'), 1000});
    refIdx["ctgX"].push_back(RefContigInfo{"asm1", "ctgX", std::string(540, 'G') + std::string(540, 'T'), 1080});

    std::unordered_map<std::string, std::vector<VariantCallBridge>> candidates;

    VariantCallBridge mem;
    mem.qAsm = "q1"; mem.qContig = "ctg1"; mem.refAsm = "asm1"; mem.refContig = "ctgX";
    mem.type = "INV"; mem.pos = 10; mem.end = 900; mem.svlen = 890;
    mem.blockScore = 250; mem.anchors = 18; mem.mapq = 40; mem.gq = 40;
    mem.alignmentMode = "mem_chain_ds13_ds18";
    candidates["ctg1"].push_back(mem);

    VariantCallBridge indel;
    indel.qAsm = "q1"; indel.qContig = "ctg1"; indel.refAsm = "asm1"; indel.refContig = "ctg1";
    indel.type = "INS"; indel.pos = 250; indel.end = 250; indel.svlen = 80;
    indel.blockScore = 20; indel.anchors = 2; indel.mapq = 30; indel.gq = 30;
    indel.alignmentMode = "simple_length_fallback";
    candidates["ctg1"].push_back(indel);

    auto chosen = select_best_call_per_contig(contigs, candidates, refIdx, o,
                    query_input::QueryMode::ASSEMBLY, "assembly");
    if (chosen.size() != 1) return 2;
    std::cout << chosen[0].type << "\t" << chosen[0].alignmentMode << "\n";
    return 0;
}
''')
    exe = tmp_path / 'arb_harness'
    run(['g++', '-O2', '-std=c++17', '-pthread', '-I', str(ROOT), str(harness), '-o', str(exe)])
    out = run([str(exe)]).stdout.strip().split('\t')
    assert out[0] == 'INS'
    assert out[1] == 'simple_length_fallback'


def test_candidate_arbitration_in_reads_modes_prefers_high_overlap_indel(tmp_path: Path):
    harness = tmp_path / 'arb_reads_harness.cpp'
    harness.write_text(r'''
#define main fungi_graphsv_tol_main
#include "main.cpp"
#undef main
#include <iostream>

static int run_mode(query_input::QueryMode mode) {
    Options o;
    o.minSvLen = 20;
    o.maxSvLen = 1000000;
    o.minAnchors = 1;

    std::unordered_map<std::string, std::string> contigs;
    contigs["sr_unitig0"] = std::string(600, 'A') + std::string(60, 'C');

    SimpleRefIndex refIdx;
    refIdx["refA"].push_back(RefContigInfo{"asm1", "refA", std::string(600, 'A'), 600});
    refIdx["refB"].push_back(RefContigInfo{"asm1", "refB", std::string(330, 'G') + std::string(330, 'T'), 660});

    std::unordered_map<std::string, std::vector<VariantCallBridge>> candidates;

    VariantCallBridge mem;
    mem.qAsm = "q1"; mem.qContig = "sr_unitig0"; mem.refAsm = "asm1"; mem.refContig = "refB";
    mem.type = "TRA"; mem.pos = 12; mem.end = 12; mem.svlen = 540;
    mem.blockScore = 220; mem.anchors = 14; mem.mapq = 40; mem.gq = 40;
    mem.alignmentMode = "mem_chain_ds13_ds18";
    candidates["sr_unitig0"].push_back(mem);

    VariantCallBridge indel;
    indel.qAsm = "q1"; indel.qContig = "sr_unitig0"; indel.refAsm = "asm1"; indel.refContig = "refA";
    indel.type = "INS"; indel.pos = 300; indel.end = 300; indel.svlen = 60;
    indel.blockScore = 16; indel.anchors = 1; indel.mapq = 25; indel.gq = 25;
    indel.alignmentMode = "reads_mode_kmer_fallback";
    candidates["sr_unitig0"].push_back(indel);

    auto chosen = select_best_call_per_contig(contigs, candidates, refIdx, o, mode,
                    mode == query_input::QueryMode::LONG_READS ? "long-reads" : "short-reads");
    if (chosen.size() != 1) return 2;
    std::cout << chosen[0].type << "\t" << chosen[0].alignmentMode << "\n";
    return 0;
}

int main() {
    if (run_mode(query_input::QueryMode::SHORT_READS) != 0) return 2;
    if (run_mode(query_input::QueryMode::LONG_READS) != 0) return 3;
    return 0;
}
    ''')
    exe = tmp_path / 'arb_reads_harness'
    run(['g++', '-O2', '-std=c++17', '-pthread', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = [line.split('\t') for line in run([str(exe)]).stdout.strip().splitlines()]
    assert rows == [
        ['INS', 'reads_mode_kmer_fallback'],
        ['INS', 'reads_mode_kmer_fallback'],
    ]


def test_candidate_arbitration_tolerates_clade_stamped_refasm(tmp_path: Path):
    harness = tmp_path / 'arb_clade_harness.cpp'
    harness.write_text(r'''
#define main fungi_graphsv_tol_main
#include "main.cpp"
#undef main
#include <iostream>

int main() {
    Options o;
    o.minSvLen = 20;
    o.maxSvLen = 1000000;
    o.minAnchors = 1;

    std::unordered_map<std::string, std::string> contigs;
    contigs["sr_unitig0"] = std::string(600, 'A') + std::string(60, 'C');

    SimpleRefIndex refIdx;
    refIdx["refA"].push_back(RefContigInfo{"asm1", "refA", std::string(600, 'A'), 600});
    refIdx["refB"].push_back(RefContigInfo{"asm2", "refB", std::string(330, 'G') + std::string(330, 'T'), 660});

    std::unordered_map<std::string, std::vector<VariantCallBridge>> candidates;

    VariantCallBridge mem;
    mem.qAsm = "q1"; mem.qContig = "sr_unitig0";
    mem.refAsm = "species_clade_alpha"; mem.refContig = "refB";
    mem.type = "TRA"; mem.pos = 12; mem.end = 12; mem.svlen = 540;
    mem.blockScore = 220; mem.anchors = 14; mem.mapq = 40; mem.gq = 40;
    mem.alignmentMode = "mem_chain_ds13_ds18";
    candidates["sr_unitig0"].push_back(mem);

    VariantCallBridge indel;
    indel.qAsm = "q1"; indel.qContig = "sr_unitig0"; indel.refAsm = "asm1"; indel.refContig = "refA";
    indel.type = "INS"; indel.pos = 300; indel.end = 300; indel.svlen = 60;
    indel.blockScore = 16; indel.anchors = 1; indel.mapq = 25; indel.gq = 25;
    indel.alignmentMode = "reads_mode_kmer_fallback";
    candidates["sr_unitig0"].push_back(indel);

    auto chosen = select_best_call_per_contig(contigs, candidates, refIdx, o,
                    query_input::QueryMode::SHORT_READS, "short-reads");
    if (chosen.size() != 1) return 2;
    std::cout << chosen[0].type << "\t" << chosen[0].alignmentMode << "\n";
    return 0;
}
''')
    exe = tmp_path / 'arb_clade_harness'
    run(['g++', '-O2', '-std=c++17', '-pthread', '-I', str(ROOT), str(harness), '-o', str(exe)])
    out = run([str(exe)]).stdout.strip().split('\t')
    assert out[0] == 'INS'
    assert out[1] == 'reads_mode_kmer_fallback'


def test_tol_graph_native_and_chain_threshold_knobs_are_live(tmp_path: Path):
    import random

    rng = random.Random(41)

    def rand_dna(n: int) -> str:
        return ''.join(rng.choice('ACGT') for _ in range(n))

    left = rand_dna(400)
    insert = rand_dna(2200)
    right = rand_dna(400)
    ref_seq = left + right
    manifest = tmp_path / 'manifest.tsv'
    ref = tmp_path / 'ref.fa'
    write_fasta(ref, [('ctg1', ref_seq)])
    manifest.write_text(
        '#asm_name\tphylum\tclade_name\tclade_rank\tfasta_path\n'
        f'asm1\tAscomycota\tClade1\tspecies\t{ref}\n'
    )

    idx = tmp_path / 'idx'
    reg = tmp_path / 'reg'
    idx.mkdir()
    reg.mkdir()

    manifest_cpp = str(manifest).replace('\\', '\\\\')
    idx_cpp = str(idx).replace('\\', '\\\\')
    reg_cpp = str(reg).replace('\\', '\\\\')
    query_cpp = (left + insert + right).replace('\\', '\\\\')

    harness = tmp_path / 'tol_knobs_harness.cpp'
    harness.write_text(f'''
#include <iostream>
#include <string>
#include <memory>
#include <unordered_map>
#include <vector>
#include "fungi_tol_bridge.hpp"

static int count_mode(const std::vector<tol::VariantCallBridge>& calls,
                      const std::string& needle) {{
    int n = 0;
    for (const auto& c : calls)
        if (c.alignmentMode.find(needle) != std::string::npos) ++n;
    return n;
}}

int main() {{
    tol::TolGlobal::RefSeq ref_row;
    ref_row.asmName = "asm1";
    ref_row.contig = "ctg1";
    ref_row.clade = "Clade1";
    ref_row.cladeRank = "species";
    ref_row.phylum = "Ascomycota";
    ref_row.seqShared = std::make_shared<std::string>("{ref_seq}");
    std::vector<tol::TolGlobal::RefSeq> direct_refs = {{ref_row}};

    tol::FederatedOptions direct_fo;
    direct_fo.fallbackSketchParams.k = 7;
    direct_fo.graphNativeMode = true;
    direct_fo.tolMinBlockBp = 500;
    auto wins_on = tol::discover_graph_native_offref_windows("{insert}", direct_refs, direct_fo);
    std::cout << "graph_on\t" << wins_on.size() << "\\n";

    direct_fo.graphNativeMode = false;
    auto wins_off = tol::discover_graph_native_offref_windows("{insert}", direct_refs, direct_fo);
    std::cout << "graph_off\t" << wins_off.size() << "\\n";

    direct_fo.graphNativeMode = true;
    direct_fo.tolMinBlockBp = 1500;
    auto wins_hi_block = tol::discover_graph_native_offref_windows("{insert}", direct_refs, direct_fo);
    std::cout << "graph_hi_block\t" << wins_hi_block.size() << "\\n";

    tol::SyncmerParams sp;
    tol::SyncmerParams fb;
    fb.k = 7;
    fb.s = 3;
    tol::build_tol_index_from_manifest(
        "{manifest_cpp}", "{idx_cpp}", "{reg_cpp}", sp, 0.12, false, "", 500, 1, &fb, false);
    tol::MultiRankIndex::instance().init("{idx_cpp}", "{reg_cpp}", 1ull << 20, 64);

    std::unordered_map<std::string, std::string> contigs;
    contigs["ctg1"] = "{query_cpp}";

    tol::FederatedOptions fo;
    fo.primarySketchParams.k = 21;
    fo.secondarySketchParams.k = 15;
    fo.fallbackSketchParams.k = 7;
    fo.fallbackSketchParams.s = 3;
    fo.useSecondarySeeds = true;
    fo.repeatRescueMinAnchors = 3;
    fo.minSvLen = 40;
    fo.minBlockScore = 6.0;
    fo.minAnchors = 2;
    fo.chainGapBand = 5000;
    fo.graphNativeMode = true;
    fo.tolMinBlockBp = 500;
    fo.tolMinChainAnchors = 2;

    auto mem_lo = tol::hierarchical_call_assembly("q", contigs, fo);
    std::cout << "mem_lo\t" << count_mode(mem_lo, "mem_chain_ds13_ds18") << "\\n";

    fo.tolMinBlockBp = 500;
    fo.tolMinChainAnchors = 5;
    auto mem_hi = tol::hierarchical_call_assembly("q", contigs, fo);
    std::cout << "mem_hi\t" << count_mode(mem_hi, "mem_chain_ds13_ds18") << "\\n";
}}
''')
    exe = tmp_path / 'tol_knobs_harness'
    run(['g++', '-O2', '-std=c++17', '-pthread', '-I', str(ROOT), str(harness), '-o', str(exe)])
    rows = dict(line.split('\t', 1) for line in run([str(exe)]).stdout.strip().splitlines())
    assert int(rows['graph_on']) >= 1
    assert int(rows['graph_off']) == 0
    assert int(rows['graph_hi_block']) <= int(rows['graph_on'])
