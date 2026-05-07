#!/usr/bin/env python3
"""
Detect candidate intron-containing 16S rRNA sequences from BLASTN HSP results.

Reference modes:
  1. --blast: use an existing BLASTN outfmt 6 table.
  2. --db: use an existing BLAST database and run BLASTN internally.
  3. --ref-fasta: preprocess a SILVA-style FASTA/FASTA.gz reference, build a
     BLAST database, and run BLASTN internally.

SILVA-style header example:
  >AB000393.1.1510 Bacteria;Pseudomonadota;...;Vibrio;Vibrio halioticoli

Core algorithm: hsp-gap-support
  query:    HSP1 ---- candidate intron ---- HSP2
  subject:  HSP1 ---- nearly continuous ---- HSP2

Default important parameters:
  --min-pident              75.0
  --min-hsp-len             100
  --min-intron-len          25
  --max-intron-len          2000
  --max-ref-gap             30
  --breakpoint-window       30
  --top-subjects            100
  --ref-domains             Archaea,Bacteria
  --ref-per-species         5
  --threads                 4
"""

from __future__ import annotations

import argparse
import concurrent.futures as futures
import csv
import gzip
import logging
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, TextIO, Tuple

from rich.console import Console
from rich.logging import RichHandler
from rich.progress import BarColumn, Progress, SpinnerColumn, TaskProgressColumn, TextColumn, TimeElapsedColumn, TimeRemainingColumn
from rich.table import Table
from rich_argparse import RichHelpFormatter


CONSOLE = Console(stderr=True)
LOG = logging.getLogger("detect_16s_introns_from_blast")


@dataclass(frozen=True)
class HSP:
    """One BLASTN high-scoring pair."""

    qseqid: str
    sseqid: str
    pident: float
    length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: str
    bitscore: float

    @property
    def qlo(self) -> int:
        return min(self.qstart, self.qend)

    @property
    def qhi(self) -> int:
        return max(self.qstart, self.qend)

    @property
    def slo(self) -> int:
        return min(self.sstart, self.send)

    @property
    def shi(self) -> int:
        return max(self.sstart, self.send)

    @property
    def qdir(self) -> int:
        return 1 if self.qend >= self.qstart else -1

    @property
    def sdir(self) -> int:
        return 1 if self.send >= self.sstart else -1

    @property
    def orientation(self) -> int:
        return self.qdir * self.sdir


@dataclass(frozen=True)
class SupportPair:
    """Best intron-like HSP pair for one query-subject comparison."""

    query_id: str
    subject_id: str
    hsp1: HSP
    hsp2: HSP
    query_gap: int
    subject_gap: int
    intron_start: int
    intron_end: int
    intron_len: int
    pair_score: float
    taxonomy: str
    taxon_at_rank: str


@dataclass
class QueryResult:
    """Final intron detection result for one query."""

    query_id: str
    query_len: int
    classification: str
    confidence: str
    intron_start: int = 0
    intron_end: int = 0
    intron_len: int = 0
    exon1_start: int = 0
    exon1_end: int = 0
    exon2_start: int = 0
    exon2_end: int = 0
    intron_free_len: int = 0
    support_subjects: int = 0
    support_taxa: int = 0
    support_species: int = 0
    support_genera: int = 0
    median_subject_gap: float = 0.0
    median_pident: float = 0.0
    mean_bitscore: float = 0.0
    best_subject: str = ""
    best_subject_taxonomy: str = ""
    reasons: Optional[List[str]] = None
    support_pairs: Optional[List[SupportPair]] = None

    def __post_init__(self) -> None:
        if self.reasons is None:
            self.reasons = []
        if self.support_pairs is None:
            self.support_pairs = []


def setup_logging(verbose: bool) -> None:
    """Configure rich colored logging."""

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(console=CONSOLE, rich_tracebacks=True, show_path=False)],
    )


def output_path(prefix: Path, suffix: str) -> Path:
    """Return PREFIX + suffix without replacing dotted suffixes."""

    return Path(str(prefix) + suffix)


def fasta_label(path: Path) -> str:
    """Create a compact label from FASTA/FASTA.gz filename."""

    name = path.name
    for suffix in [".fasta.gz", ".fa.gz", ".fna.gz", ".faa.gz", ".fasta", ".fa", ".fna", ".faa", ".gz"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    chars = [c if (c.isalnum() or c in {".", "_", "-"}) else "_" for c in name]
    return "".join(chars).strip("._-") or "query"


def setup_output_paths(args: argparse.Namespace) -> None:
    """Create output directories and deterministic output file paths."""

    args.reference_dir = args.outdir / "reference"
    args.blast_dir = args.outdir / "blast"
    args.results_dir = args.outdir / "results"
    for directory in [args.outdir, args.reference_dir, args.blast_dir, args.results_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    args.query_label = fasta_label(args.query)
    args.raw_ref_fa = args.reference_dir / "raw_reference.fa"
    args.raw_ref_tax = args.reference_dir / "raw_reference.tax.tsv"
    args.raw_ref_db = args.reference_dir / "raw_reference_db"
    args.ref_self_blast = args.reference_dir / "reference_self.blastn.tsv"
    args.ref_self_clean_prefix = args.reference_dir / "reference_self_clean"
    args.cleaned_ref_fa = args.reference_dir / "cleaned_reference.fa"
    args.cleaned_ref_tax = args.reference_dir / "cleaned_reference.tax.tsv"
    args.cleaned_ref_db = args.reference_dir / "cleaned_reference_db"

    args.query_blast = args.blast_dir / f"{args.query_label}.vs_reference.blastn.tsv"
    args.result_prefix = args.results_dir / args.query_label
    args.summary_tsv = output_path(args.result_prefix, ".summary.tsv")
    args.supporting_tsv = output_path(args.result_prefix, ".supporting_hsps.tsv")
    args.intron_free_fa = output_path(args.result_prefix, ".intron_free.fa")
    args.introns_fa = output_path(args.result_prefix, ".introns.fa")
    args.exons_bed = output_path(args.result_prefix, ".exons.bed")
    args.introns_bed = output_path(args.result_prefix, ".introns.bed")
    args.report_md = output_path(args.result_prefix, ".report.md")


def open_text_auto(path: Path) -> TextIO:
    """Open plain text or gzip-compressed text."""

    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def wrap_fasta(seq: str, width: int = 80) -> Iterable[str]:
    """Yield wrapped FASTA sequence lines."""

    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def write_fasta_record(handle: TextIO, seq_id: str, seq: str) -> None:
    """Write one FASTA record safely."""

    print(f">{seq_id}", file=handle)
    for part in wrap_fasta(seq):
        print(part, file=handle)


def read_fasta(path: Path) -> Dict[str, str]:
    """Read FASTA or FASTA.gz into {sequence_id: sequence}."""

    seqs: Dict[str, List[str]] = {}
    current_id: Optional[str] = None
    with open_text_auto(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                seqs[current_id] = []
            elif current_id is not None:
                seqs[current_id].append(line)
    return {seq_id: "".join(parts).upper() for seq_id, parts in seqs.items()}


def parse_silva_header(header: str) -> Tuple[str, str]:
    """Parse SILVA FASTA header into subject ID and taxonomy."""

    text = header[1:] if header.startswith(">") else header
    parts = text.strip().split(maxsplit=1)
    seq_id = parts[0]
    taxonomy = parts[1].strip() if len(parts) > 1 else ""
    return seq_id, taxonomy


def taxonomy_domain(taxonomy: str) -> str:
    """Return first taxonomy field."""

    parts = [p.strip() for p in taxonomy.split(";") if p.strip()]
    return parts[0] if parts else "NA"


def species_key_from_taxonomy(taxonomy: str) -> str:
    """Return species-level key from taxonomy."""

    parts = [p.strip() for p in taxonomy.split(";") if p.strip()]
    if len(parts) >= 7:
        return ";".join(parts[:7])
    return ";".join(parts) if parts else "NA"


def clean_dna_sequence(seq_parts: List[str]) -> str:
    """Join sequence lines, normalize U to T, and remove whitespace."""

    seq = "".join(seq_parts).upper().replace("U", "T")
    return "".join(seq.split())


def is_strict_atgc(seq: str) -> bool:
    """Return True only if sequence is non-empty and contains only A/T/G/C."""

    return bool(seq) and set(seq) <= {"A", "T", "G", "C"}


def preprocess_reference(args: argparse.Namespace) -> Tuple[Path, Path, Path]:
    """Filter SILVA reference, write clean FASTA/taxonomy, and build BLAST DB."""

    allowed_domains = {x.strip() for x in args.ref_domains.split(",") if x.strip()}
    per_species_counts: Dict[str, int] = defaultdict(int)
    kept = 0
    skipped_domain = 0
    skipped_species_cap = 0
    skipped_non_atgc = 0
    skipped_empty = 0
    current_header = ""
    current_seq: List[str] = []

    with open_text_auto(args.ref_fasta) as ref_in, args.raw_ref_fa.open("wt", encoding="utf-8") as fa_out, args.raw_ref_tax.open("wt", encoding="utf-8", newline="") as tax_file:
        tax_writer = csv.writer(tax_file, delimiter=chr(9), lineterminator=chr(10))

        def flush_record() -> None:
            nonlocal kept, skipped_domain, skipped_species_cap, skipped_non_atgc, skipped_empty, current_header, current_seq
            if not current_header:
                return
            seq_id, taxonomy = parse_silva_header(current_header)
            domain = taxonomy_domain(taxonomy)
            if domain not in allowed_domains:
                skipped_domain += 1
                return
            seq = clean_dna_sequence(current_seq)
            if not seq:
                skipped_empty += 1
                return
            if not is_strict_atgc(seq):
                skipped_non_atgc += 1
                return
            species_key = species_key_from_taxonomy(taxonomy)
            if args.ref_per_species > 0 and per_species_counts[species_key] >= args.ref_per_species:
                skipped_species_cap += 1
                return
            per_species_counts[species_key] += 1
            write_fasta_record(fa_out, seq_id, seq)
            tax_writer.writerow([seq_id, taxonomy])
            kept += 1

        for raw in ref_in:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_record()
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        flush_record()

    LOG.info("Reference preprocessing kept %s sequences", kept)
    LOG.info("Skipped by domain filter: %s", skipped_domain)
    LOG.info("Skipped by per-species cap: %s", skipped_species_cap)
    LOG.info("Skipped empty sequences: %s", skipped_empty)
    LOG.info("Skipped sequences containing non-ATGC characters: %s", skipped_non_atgc)
    make_blast_db(args.raw_ref_fa, args.raw_ref_db, args.makeblastdb_bin)
    return args.raw_ref_fa, args.raw_ref_tax, args.raw_ref_db


def make_blast_db(fasta: Path, db_prefix: Path, makeblastdb_bin: str) -> None:
    """Build nucleotide BLAST database."""

    cmd = [makeblastdb_bin, "-in", str(fasta), "-dbtype", "nucl", "-out", str(db_prefix)]
    LOG.info("Running makeblastdb: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_taxonomy(path: Optional[Path]) -> Dict[str, str]:
    """Read subject taxonomy mapping."""

    if path is None:
        return {}
    taxonomy: Dict[str, str] = {}
    with path.open("rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=chr(9))
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 2:
                continue
            taxonomy[row[0]] = row[1]
    LOG.info("Parsed taxonomy mappings: %s subjects", len(taxonomy))
    return taxonomy


def run_blastn_to_file(query: Path, db: Path, out_file: Path, args: argparse.Namespace, max_targets: int, max_hsps: int, label: str) -> Path:
    """Run BLASTN and return output table path."""

    cmd = [
        args.blastn_bin,
        "-query",
        str(query),
        "-db",
        str(db),
        "-outfmt",
        "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
        "-max_target_seqs",
        str(max_targets),
        "-max_hsps",
        str(max_hsps),
        "-num_threads",
        str(args.threads),
        "-out",
        str(out_file),
    ]
    if args.blast_task:
        cmd.extend(["-task", args.blast_task])
    if args.blast_evalue:
        cmd.extend(["-evalue", str(args.blast_evalue)])
    LOG.info("Running %s: %s", label, " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out_file


def run_query_blastn(args: argparse.Namespace) -> Path:
    """Run query BLASTN and return table path."""

    return run_blastn_to_file(
        query=args.query,
        db=args.db,
        out_file=args.query_blast,
        args=args,
        max_targets=args.blast_max_target_seqs,
        max_hsps=args.blast_max_hsps,
        label="BLASTN",
    )


def parse_blast(path: Path, min_pident: float, min_hsp_len: int) -> Dict[str, Dict[str, List[HSP]]]:
    """Parse BLAST outfmt 6 and group HSPs by query and subject."""

    grouped: Dict[str, Dict[str, List[HSP]]] = defaultdict(lambda: defaultdict(list))
    n_lines = 0
    n_kept = 0
    with path.open("rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=chr(9))
        for parts in reader:
            if not parts or parts[0].startswith("#"):
                continue
            n_lines += 1
            if len(parts) < 10:
                LOG.debug("Skipping BLAST row with <10 columns: %s", parts)
                continue
            pident = float(parts[2])
            length = int(parts[3])
            if pident < min_pident or length < min_hsp_len:
                continue
            hsp = HSP(
                qseqid=parts[0],
                sseqid=parts[1],
                pident=pident,
                length=length,
                qstart=int(parts[4]),
                qend=int(parts[5]),
                sstart=int(parts[6]),
                send=int(parts[7]),
                evalue=parts[8],
                bitscore=float(parts[9]),
            )
            grouped[hsp.qseqid][hsp.sseqid].append(hsp)
            n_kept += 1
    LOG.info("Parsed BLAST HSPs: %s rows, %s kept after filters", n_lines, n_kept)
    return grouped


def get_taxon_at_rank(taxonomy: str, rank: str) -> str:
    """Extract taxon at a rank from semicolon-delimited taxonomy."""

    if not taxonomy:
        return "NA"
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    parts = [p.strip() for p in taxonomy.split(";") if p.strip()]
    index = {name: i for i, name in enumerate(ranks)}.get(rank, 6)
    if index < len(parts):
        return parts[index]
    return parts[-1] if parts else "NA"


def subject_gap_by_query_order(left: HSP, right: HSP) -> Optional[int]:
    """Compute subject gap for two HSPs ordered by query coordinate."""

    if left.orientation != right.orientation:
        return None
    if left.orientation == 1:
        if left.slo > right.slo:
            return None
        return right.slo - left.shi - 1
    if left.orientation == -1:
        if left.slo < right.slo:
            return None
        return left.slo - right.shi - 1
    return None


def hsp_pair_supports_intron(h1: HSP, h2: HSP, args: argparse.Namespace) -> Optional[Tuple[int, int, int, int, int]]:
    """Return intron geometry if two HSPs support an intron-like insertion."""

    left, right = sorted([h1, h2], key=lambda h: (h.qlo, h.qhi))
    query_gap = right.qlo - left.qhi - 1
    if query_gap < 0 and abs(query_gap) > args.max_query_overlap:
        return None
    if query_gap < args.min_intron_len or query_gap > args.max_intron_len:
        return None
    subject_gap = subject_gap_by_query_order(left, right)
    if subject_gap is None:
        return None
    if abs(subject_gap) > args.max_ref_gap:
        return None
    intron_start = left.qhi + 1
    intron_end = right.qlo - 1
    intron_len = intron_end - intron_start + 1
    if intron_len < args.min_intron_len or intron_len > args.max_intron_len:
        return None
    return query_gap, subject_gap, intron_start, intron_end, intron_len


def best_support_pair_for_subject(query_id: str, subject_id: str, hsps: List[HSP], taxonomy: Dict[str, str], args: argparse.Namespace) -> Optional[SupportPair]:
    """Find best intron-like HSP pair for one query-subject pair."""

    if len(hsps) < 2:
        return None
    best_pair: Optional[SupportPair] = None
    tax = taxonomy.get(subject_id, "")
    taxon = get_taxon_at_rank(tax, args.tax_rank)
    hsps_sorted = sorted(hsps, key=lambda h: (-h.bitscore, h.qlo, h.qhi))
    for i in range(len(hsps_sorted)):
        for j in range(i + 1, len(hsps_sorted)):
            support = hsp_pair_supports_intron(hsps_sorted[i], hsps_sorted[j], args)
            if support is None:
                continue
            query_gap, subject_gap, intron_start, intron_end, intron_len = support
            pair_score = hsps_sorted[i].bitscore + hsps_sorted[j].bitscore - 2.0 * abs(subject_gap)
            candidate = SupportPair(
                query_id=query_id,
                subject_id=subject_id,
                hsp1=hsps_sorted[i],
                hsp2=hsps_sorted[j],
                query_gap=query_gap,
                subject_gap=subject_gap,
                intron_start=intron_start,
                intron_end=intron_end,
                intron_len=intron_len,
                pair_score=pair_score,
                taxonomy=tax,
                taxon_at_rank=taxon,
            )
            if best_pair is None or candidate.pair_score > best_pair.pair_score:
                best_pair = candidate
    return best_pair


def top_subjects_by_bitscore(subject_hsps: Dict[str, List[HSP]], top_subjects: int) -> Dict[str, List[HSP]]:
    """Keep top subjects by summed HSP bitscore."""

    ranked = sorted(subject_hsps.items(), key=lambda item: sum(h.bitscore for h in item[1]), reverse=True)
    return dict(ranked[:top_subjects])


def cluster_support_pairs(pairs: List[SupportPair], breakpoint_window: int) -> List[List[SupportPair]]:
    """Cluster support pairs by similar query intron coordinates."""

    clusters: List[List[SupportPair]] = []
    for pair in sorted(pairs, key=lambda p: (p.intron_start, p.intron_end, -p.pair_score)):
        assigned = False
        for cluster in clusters:
            med_start = int(median([p.intron_start for p in cluster]))
            med_end = int(median([p.intron_end for p in cluster]))
            if abs(pair.intron_start - med_start) <= breakpoint_window and abs(pair.intron_end - med_end) <= breakpoint_window:
                cluster.append(pair)
                assigned = True
                break
        if not assigned:
            clusters.append([pair])
    return clusters


def unique_taxa_at_rank(pairs: List[SupportPair], rank: str) -> int:
    """Count unique taxa at a rank."""

    taxa = set()
    for pair in pairs:
        taxon = get_taxon_at_rank(pair.taxonomy, rank)
        if taxon != "NA":
            taxa.add(taxon)
    return len(taxa)


def classify_confidence(support_subjects: int, support_taxa: int, args: argparse.Namespace) -> Tuple[str, str, List[str]]:
    """Assign confidence label from support counts."""

    if support_subjects >= args.high_support_subjects and support_taxa >= args.high_support_taxa:
        return "HIGH_CONFIDENCE_16S_INTRON", "HIGH", [f"high_support_subjects>={args.high_support_subjects}", f"high_support_taxa>={args.high_support_taxa}"]
    if support_subjects >= args.medium_support_subjects and support_taxa >= args.medium_support_taxa:
        return "MEDIUM_CONFIDENCE_16S_INTRON", "MEDIUM", [f"medium_support_subjects>={args.medium_support_subjects}", f"medium_support_taxa>={args.medium_support_taxa}"]
    if support_subjects >= args.min_support_subjects:
        return "LOW_CONFIDENCE_16S_INTRON", "LOW", [f"min_support_subjects>={args.min_support_subjects}"]
    return "NO_INTRON_SIGNAL", "NONE", ["no_supported_intron_cluster"]


def analyze_query(query_id: str, subject_hsps: Dict[str, List[HSP]], query_len: int, taxonomy: Dict[str, str], args: argparse.Namespace) -> QueryResult:
    """Analyze one query sequence."""

    top_subjects = top_subjects_by_bitscore(subject_hsps, args.top_subjects)
    support_pairs: List[SupportPair] = []
    for subject_id, hsps in top_subjects.items():
        pair = best_support_pair_for_subject(query_id, subject_id, hsps, taxonomy, args)
        if pair is not None:
            support_pairs.append(pair)
    if not support_pairs:
        return QueryResult(query_id=query_id, query_len=query_len, classification="NO_INTRON_SIGNAL", confidence="NONE", reasons=["no_subject_supported_hsp_gap_pattern"])

    clusters = cluster_support_pairs(support_pairs, args.breakpoint_window)
    clusters.sort(key=lambda c: (len(c), len({p.taxon_at_rank for p in c if p.taxon_at_rank != "NA"}), sum(p.pair_score for p in c)), reverse=True)
    best_cluster = clusters[0]
    best_pair = sorted(best_cluster, key=lambda p: p.pair_score, reverse=True)[0]
    intron_start = best_pair.intron_start
    intron_end = best_pair.intron_end
    intron_len = intron_end - intron_start + 1
    support_subjects = len({p.subject_id for p in best_cluster})
    support_taxa = len({p.taxon_at_rank for p in best_cluster if p.taxon_at_rank != "NA"})
    support_species = unique_taxa_at_rank(best_cluster, "species")
    support_genera = unique_taxa_at_rank(best_cluster, "genus")
    classification, confidence, reasons = classify_confidence(support_subjects, support_taxa, args)
    reasons.extend([
        f"algorithm={args.algorithm}",
        "coordinate_source=best_subject_pair",
        f"support_subjects={support_subjects}",
        f"support_taxa_at_{args.tax_rank}={support_taxa}",
        f"breakpoint_window={args.breakpoint_window}bp",
    ])
    median_subject_gap = float(median([p.subject_gap for p in best_cluster]))
    median_pident = float(median([(p.hsp1.pident + p.hsp2.pident) / 2.0 for p in best_cluster]))
    mean_bitscore = sum(p.hsp1.bitscore + p.hsp2.bitscore for p in best_cluster) / len(best_cluster)
    exon1_start = 1
    exon1_end = intron_start - 1
    exon2_start = intron_end + 1
    exon2_end = query_len
    intron_free_len = max(0, exon1_end - exon1_start + 1) + max(0, exon2_end - exon2_start + 1)
    return QueryResult(
        query_id=query_id,
        query_len=query_len,
        classification=classification,
        confidence=confidence,
        intron_start=intron_start,
        intron_end=intron_end,
        intron_len=intron_len,
        exon1_start=exon1_start,
        exon1_end=exon1_end,
        exon2_start=exon2_start,
        exon2_end=exon2_end,
        intron_free_len=intron_free_len,
        support_subjects=support_subjects,
        support_taxa=support_taxa,
        support_species=support_species,
        support_genera=support_genera,
        median_subject_gap=median_subject_gap,
        median_pident=median_pident,
        mean_bitscore=mean_bitscore,
        best_subject=best_pair.subject_id,
        best_subject_taxonomy=best_pair.taxonomy,
        reasons=reasons,
        support_pairs=best_cluster,
    )


def confidence_rank(label: str) -> int:
    """Rank confidence labels."""

    return {"NONE": 0, "LOW": 1, "MEDIUM": 2, "HIGH": 3}.get(label, 0)


def is_intron_result(result: QueryResult, min_confidence: str) -> bool:
    """Return True if a result passes output confidence threshold."""

    ranks = {"LOW": 1, "MEDIUM": 2, "HIGH": 3}
    return ranks.get(result.confidence, 0) >= ranks.get(min_confidence, 1)


def summary_row(result: QueryResult, tax_rank: str) -> Dict[str, str]:
    """Convert one query result to a TSV row."""

    return {
        "query_id": result.query_id,
        "query_len": str(result.query_len),
        "classification": result.classification,
        "confidence": result.confidence,
        "intron_start": str(result.intron_start) if result.intron_start else "",
        "intron_end": str(result.intron_end) if result.intron_end else "",
        "intron_len": str(result.intron_len) if result.intron_len else "",
        "exon1": f"{result.exon1_start}-{result.exon1_end}" if result.exon1_start else "",
        "exon2": f"{result.exon2_start}-{result.exon2_end}" if result.exon2_start else "",
        "intron_free_len": str(result.intron_free_len) if result.intron_free_len else "",
        "support_subjects": str(result.support_subjects),
        f"support_taxa_at_{tax_rank}": str(result.support_taxa),
        "support_species": str(result.support_species),
        "support_genera": str(result.support_genera),
        "median_subject_gap": f"{result.median_subject_gap:.2f}",
        "median_pident": f"{result.median_pident:.2f}",
        "mean_bitscore": f"{result.mean_bitscore:.2f}",
        "best_subject": result.best_subject,
        "best_subject_taxonomy": result.best_subject_taxonomy,
        "reasons": "|".join(result.reasons or []),
    }


def support_row(pair: SupportPair) -> Dict[str, str]:
    """Convert one support pair to a TSV row."""

    h1, h2 = sorted([pair.hsp1, pair.hsp2], key=lambda h: (h.qlo, h.qhi))
    return {
        "query_id": pair.query_id,
        "subject_id": pair.subject_id,
        "taxonomy": pair.taxonomy,
        "taxon_at_rank": pair.taxon_at_rank,
        "intron_start": str(pair.intron_start),
        "intron_end": str(pair.intron_end),
        "intron_len": str(pair.intron_len),
        "query_gap": str(pair.query_gap),
        "subject_gap": str(pair.subject_gap),
        "hsp1_q": f"{h1.qlo}-{h1.qhi}",
        "hsp2_q": f"{h2.qlo}-{h2.qhi}",
        "hsp1_s": f"{h1.slo}-{h1.shi}",
        "hsp2_s": f"{h2.slo}-{h2.shi}",
        "hsp1_pident": f"{h1.pident:.2f}",
        "hsp2_pident": f"{h2.pident:.2f}",
        "hsp1_len": str(h1.length),
        "hsp2_len": str(h2.length),
        "bitscore_sum": f"{h1.bitscore + h2.bitscore:.2f}",
        "pair_score": f"{pair.pair_score:.2f}",
    }


def write_summary(path: Path, results: List[QueryResult], tax_rank: str) -> None:
    """Write query-level summary TSV."""

    fields = [
        "query_id",
        "query_len",
        "classification",
        "confidence",
        "intron_start",
        "intron_end",
        "intron_len",
        "exon1",
        "exon2",
        "intron_free_len",
        "support_subjects",
        f"support_taxa_at_{tax_rank}",
        "support_species",
        "support_genera",
        "median_subject_gap",
        "median_pident",
        "mean_bitscore",
        "best_subject",
        "best_subject_taxonomy",
        "reasons",
    ]
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=chr(9), lineterminator=chr(10))
        writer.writeheader()
        for result in results:
            writer.writerow(summary_row(result, tax_rank))


def write_supporting_hsps(path: Path, results: List[QueryResult]) -> None:
    """Write subject-level supporting HSP pairs."""

    fields = [
        "query_id",
        "subject_id",
        "taxonomy",
        "taxon_at_rank",
        "intron_start",
        "intron_end",
        "intron_len",
        "query_gap",
        "subject_gap",
        "hsp1_q",
        "hsp2_q",
        "hsp1_s",
        "hsp2_s",
        "hsp1_pident",
        "hsp2_pident",
        "hsp1_len",
        "hsp2_len",
        "bitscore_sum",
        "pair_score",
    ]
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=chr(9), lineterminator=chr(10))
        writer.writeheader()
        for result in results:
            for pair in result.support_pairs or []:
                writer.writerow(support_row(pair))


def write_fasta_outputs(args: argparse.Namespace, results: List[QueryResult], seqs: Dict[str, str]) -> None:
    """Write intron-free 16S and intron FASTA outputs."""

    with args.intron_free_fa.open("wt", encoding="utf-8") as free_out, args.introns_fa.open("wt", encoding="utf-8") as intron_out:
        for result in results:
            if not is_intron_result(result, args.min_output_confidence):
                continue
            seq = seqs.get(result.query_id, "")
            if not seq:
                continue
            exon1 = seq[result.exon1_start - 1 : result.exon1_end]
            exon2 = seq[result.exon2_start - 1 : result.exon2_end]
            intron = seq[result.intron_start - 1 : result.intron_end]
            write_fasta_record(free_out, f"{result.query_id}|intron_free|intron={result.intron_start}-{result.intron_end}|confidence={result.confidence}", exon1 + exon2)
            write_fasta_record(intron_out, f"{result.query_id}|intron|{result.intron_start}-{result.intron_end}|len={result.intron_len}|confidence={result.confidence}", intron)


def write_bed_outputs(args: argparse.Namespace, results: List[QueryResult]) -> None:
    """Write BED files."""

    with args.introns_bed.open("wt", encoding="utf-8", newline="") as intron_file, args.exons_bed.open("wt", encoding="utf-8", newline="") as exon_file:
        intron_writer = csv.writer(intron_file, delimiter=chr(9), lineterminator=chr(10))
        exon_writer = csv.writer(exon_file, delimiter=chr(9), lineterminator=chr(10))
        for result in results:
            if not is_intron_result(result, args.min_output_confidence):
                continue
            intron_writer.writerow([result.query_id, result.intron_start - 1, result.intron_end, f"{result.query_id}|intron|{result.intron_start}-{result.intron_end}|confidence={result.confidence}", ".", "+"])
            exon_writer.writerow([result.query_id, result.exon1_start - 1, result.exon1_end, f"{result.query_id}|exon1|{result.exon1_start}-{result.exon1_end}|confidence={result.confidence}", ".", "+"])
            exon_writer.writerow([result.query_id, result.exon2_start - 1, result.exon2_end, f"{result.query_id}|exon2|{result.exon2_start}-{result.exon2_end}|confidence={result.confidence}", ".", "+"])


def write_report(path: Path, results: List[QueryResult], args: argparse.Namespace) -> None:
    """Write Markdown report."""

    counts = Counter(r.confidence for r in results)
    class_counts = Counter(r.classification for r in results)
    intron_lens = [r.intron_len for r in results if r.intron_len > 0]

    with path.open("wt", encoding="utf-8") as handle:
        print("# BLASTN-based 16S intron detection report", file=handle)
        print("", file=handle)
        print("## Inputs", file=handle)
        print("", file=handle)
        print(f"- Query FASTA: `{args.query}`", file=handle)
        print(f"- BLAST table: `{args.blast}`", file=handle)
        print(f"- BLAST database: `{args.db}`", file=handle)
        print(f"- Reference FASTA: `{args.ref_fasta}`", file=handle)
        print(f"- Taxonomy table: `{args.taxonomy}`", file=handle)
        print(f"- Output directory: `{args.outdir}`", file=handle)
        print("", file=handle)
        print("## Algorithm", file=handle)
        print("", file=handle)
        print(f"- Algorithm: `{args.algorithm}`", file=handle)
        print("- Detects HSP1--query insertion--HSP2 patterns where reference subject coordinates are nearly continuous.", file=handle)
        print("", file=handle)
        print("## Parameters", file=handle)
        print("", file=handle)
        for name in [
            "min_pident",
            "min_hsp_len",
            "min_intron_len",
            "max_intron_len",
            "max_ref_gap",
            "max_query_overlap",
            "breakpoint_window",
            "top_subjects",
            "ref_domains",
            "ref_per_species",
            "clean_ref_introns",
            "ref_clean_min_confidence",
            "ref_self_blast_max_target_seqs",
            "ref_self_blast_max_hsps",
            "blast_max_target_seqs",
            "blast_max_hsps",
            "blast_task",
            "blast_evalue",
            "tax_rank",
            "min_support_subjects",
            "medium_support_subjects",
            "medium_support_taxa",
            "high_support_subjects",
            "high_support_taxa",
            "min_output_confidence",
            "threads",
        ]:
            print(f"- `{name}`: `{getattr(args, name)}`", file=handle)
        print("", file=handle)
        print("## Summary", file=handle)
        print("", file=handle)
        print(f"- Total query results: `{len(results)}`", file=handle)
        for label in ["HIGH", "MEDIUM", "LOW", "NONE"]:
            print(f"- {label}: `{counts.get(label, 0)}`", file=handle)
        print("", file=handle)
        print("## Classification counts", file=handle)
        print("", file=handle)
        for cls, count in class_counts.most_common():
            print(f"- `{cls}`: `{count}`", file=handle)
        if intron_lens:
            print("", file=handle)
            print("## Candidate intron length distribution", file=handle)
            print("", file=handle)
            print(f"- Min: `{min(intron_lens)}` bp", file=handle)
            print(f"- Median: `{int(round(median(intron_lens)))}` bp", file=handle)
            print(f"- Max: `{max(intron_lens)}` bp", file=handle)
        print("", file=handle)
        print("## Top candidates", file=handle)
        print("", file=handle)
        for result in sorted(results, key=lambda r: (confidence_rank(r.confidence), r.support_subjects, r.support_taxa), reverse=True)[:30]:
            if result.confidence == "NONE":
                continue
            print(f"- `{result.query_id}`: {result.confidence}, intron={result.intron_start}-{result.intron_end} ({result.intron_len} bp), support_subjects={result.support_subjects}, support_taxa={result.support_taxa}, best_subject=`{result.best_subject}`", file=handle)


def clean_reference_introns(args: argparse.Namespace, ref_fa: Path, tax_tsv: Path, db_prefix: Path) -> Tuple[Path, Path, Path]:
    """Self-BLAST reference, remove candidate introns, and build cleaned DB."""

    self_blast = run_blastn_to_file(
        query=ref_fa,
        db=db_prefix,
        out_file=args.ref_self_blast,
        args=args,
        max_targets=args.ref_self_blast_max_target_seqs,
        max_hsps=args.ref_self_blast_max_hsps,
        label="reference self-BLASTN",
    )
    ref_seqs = read_fasta(ref_fa)
    ref_taxonomy = parse_taxonomy(tax_tsv)
    ref_blast_by_query = parse_blast(self_blast, args.min_pident, args.min_hsp_len)
    ref_results: List[QueryResult] = []
    for query_id in sorted(ref_seqs):
        ref_results.append(analyze_query(query_id, ref_blast_by_query.get(query_id, {}), len(ref_seqs.get(query_id, "")), ref_taxonomy, args))
    ref_results.sort(key=lambda r: (confidence_rank(r.confidence), r.support_subjects, r.support_taxa, r.query_id), reverse=True)
    write_summary(output_path(args.ref_self_clean_prefix, ".summary.tsv"), ref_results, args.tax_rank)
    write_supporting_hsps(output_path(args.ref_self_clean_prefix, ".supporting_hsps.tsv"), ref_results)
    write_report(output_path(args.ref_self_clean_prefix, ".report.md"), ref_results, args)

    remove_by_id = {r.query_id: r for r in ref_results if is_intron_result(r, args.ref_clean_min_confidence)}
    LOG.info("Reference sequences with candidate introns to remove: %s", len(remove_by_id))
    with args.cleaned_ref_fa.open("wt", encoding="utf-8") as fa_out, args.cleaned_ref_tax.open("wt", encoding="utf-8", newline="") as tax_file:
        tax_writer = csv.writer(tax_file, delimiter=chr(9), lineterminator=chr(10))
        for seq_id in sorted(ref_seqs):
            seq = ref_seqs[seq_id]
            result = remove_by_id.get(seq_id)
            if result is not None:
                seq = seq[result.exon1_start - 1 : result.exon1_end] + seq[result.exon2_start - 1 : result.exon2_end]
            write_fasta_record(fa_out, seq_id, seq)
            tax_writer.writerow([seq_id, ref_taxonomy.get(seq_id, "")])
    make_blast_db(args.cleaned_ref_fa, args.cleaned_ref_db, args.makeblastdb_bin)
    return args.cleaned_ref_fa, args.cleaned_ref_tax, args.cleaned_ref_db


def print_summary(results: List[QueryResult], args: argparse.Namespace) -> None:
    """Print terminal summary."""

    counts = Counter(r.confidence for r in results)
    table = Table(title="BLASTN-based 16S intron detection summary")
    table.add_column("Metric", style="bold")
    table.add_column("Value", justify="right")
    table.add_row("Queries analyzed", f"{len(results):,}")
    table.add_row("High confidence", f"{counts.get('HIGH', 0):,}")
    table.add_row("Medium confidence", f"{counts.get('MEDIUM', 0):,}")
    table.add_row("Low confidence", f"{counts.get('LOW', 0):,}")
    table.add_row("No signal", f"{counts.get('NONE', 0):,}")
    table.add_row("Summary TSV", str(args.summary_tsv))
    table.add_row("Intron-free FASTA", str(args.intron_free_fa))
    table.add_row("Introns FASTA", str(args.introns_fa))
    CONSOLE.print(table)


def positive_int(value: str) -> int:
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("Value must be positive")
    return ivalue


def nonnegative_int(value: str) -> int:
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("Value must be non-negative")
    return ivalue


def probability_percent(value: str) -> float:
    fvalue = float(value)
    if fvalue < 0 or fvalue > 100:
        raise argparse.ArgumentTypeError("Value must be between 0 and 100")
    return fvalue


def build_parser() -> argparse.ArgumentParser:
    """Create argument parser."""

    parser = argparse.ArgumentParser(
        description="Detect candidate intron-containing 16S rRNA sequences from BLASTN HSP results.",
        formatter_class=RichHelpFormatter,
    )

    required = parser.add_argument_group("Required inputs")
    required.add_argument("--query", required=True, type=Path, help="Query 16S FASTA/FASTA.gz file. No default.")
    required.add_argument("--blast", default=None, type=Path, help="Existing BLASTN outfmt 6 table. Default: not used; provide exactly one of --blast, --db, or --ref-fasta.")
    required.add_argument("--db", default=None, type=Path, help="Existing BLAST database prefix. Default: not used; provide exactly one of --blast, --db, or --ref-fasta.")
    required.add_argument("--ref-fasta", default=None, type=Path, help="SILVA-style reference FASTA/FASTA.gz with taxonomy in headers. Default: not used; provide exactly one of --blast, --db, or --ref-fasta.")
    required.add_argument("--outdir", required=True, type=Path, help="Output directory. No default.")

    algorithm = parser.add_argument_group("Algorithm")
    algorithm.add_argument("--algorithm", default="hsp-gap-support", choices=["hsp-gap-support"], help="Detection algorithm. Default: hsp-gap-support.")

    filters = parser.add_argument_group("Reference preprocessing and HSP filters")
    filters.add_argument("--ref-domains", default="Archaea,Bacteria", help="Comma-separated SILVA domains retained with --ref-fasta. Default: Archaea,Bacteria.")
    filters.add_argument("--ref-per-species", default=5, type=nonnegative_int, help="Maximum sequences per species in reference preprocessing. Use 0 to disable. Default: 5.")
    filters.add_argument("--clean-ref-introns", action="store_true", help="Self-BLAST reference and remove candidate introns before query BLAST. Default: disabled.")
    filters.add_argument("--ref-clean-min-confidence", default="LOW", choices=["LOW", "MEDIUM", "HIGH"], help="Minimum confidence required to remove a reference intron. Default: LOW.")
    filters.add_argument("--ref-self-blast-max-target-seqs", default=100, type=positive_int, help="Reference self-BLAST -max_target_seqs. Default: 100.")
    filters.add_argument("--ref-self-blast-max-hsps", default=20, type=positive_int, help="Reference self-BLAST -max_hsps. Default: 20.")
    filters.add_argument("--min-pident", default=75.0, type=probability_percent, help="Minimum HSP percent identity. Default: 75.0.")
    filters.add_argument("--min-hsp-len", default=100, type=positive_int, help="Minimum HSP length in bp. Default: 100.")
    filters.add_argument("--top-subjects", default=100, type=positive_int, help="Top subjects retained per query after parsing BLAST. Default: 100.")
    filters.add_argument("--makeblastdb-bin", default="makeblastdb", help="makeblastdb executable. Default: makeblastdb.")
    filters.add_argument("--blastn-bin", default="blastn", help="blastn executable. Default: blastn.")
    filters.add_argument("--blast-max-target-seqs", default=100, type=positive_int, help="Query BLASTN -max_target_seqs. Default: 100.")
    filters.add_argument("--blast-max-hsps", default=20, type=positive_int, help="Query BLASTN -max_hsps. Default: 20.")
    filters.add_argument("--blast-task", default="blastn", choices=["blastn", "megablast", "dc-megablast", "blastn-short"], help="BLASTN task. Default: blastn.")
    filters.add_argument("--blast-evalue", default="1e-20", help="BLASTN e-value. Default: 1e-20.")

    intron = parser.add_argument_group("Candidate intron geometry")
    intron.add_argument("--min-intron-len", default=25, type=nonnegative_int, help="Minimum query gap size. Default: 25 bp.")
    intron.add_argument("--max-intron-len", default=2000, type=positive_int, help="Maximum query gap size. Default: 2000 bp.")
    intron.add_argument("--max-ref-gap", default=30, type=nonnegative_int, help="Maximum absolute reference gap/overlap. Default: 30 bp.")
    intron.add_argument("--max-query-overlap", default=20, type=nonnegative_int, help="Maximum allowed query HSP overlap. Default: 20 bp.")
    intron.add_argument("--breakpoint-window", default=30, type=nonnegative_int, help="Breakpoint clustering window. Default: 30 bp.")

    taxonomy = parser.add_argument_group("Taxonomy support")
    taxonomy.add_argument("--taxonomy", default=None, type=Path, help="Optional subject taxonomy TSV. If --ref-fasta is used, generated automatically.")
    taxonomy.add_argument("--tax-rank", default="species", choices=["domain", "phylum", "class", "order", "family", "genus", "species"], help="Taxonomic rank used for confidence support. Default: species.")

    confidence = parser.add_argument_group("Confidence thresholds")
    confidence.add_argument("--min-support-subjects", default=1, type=positive_int, help="Minimum subjects for LOW confidence. Default: 1.")
    confidence.add_argument("--medium-support-subjects", default=3, type=positive_int, help="Minimum subjects for MEDIUM confidence. Default: 3.")
    confidence.add_argument("--medium-support-taxa", default=2, type=positive_int, help="Minimum taxa for MEDIUM confidence. Default: 2.")
    confidence.add_argument("--high-support-subjects", default=10, type=positive_int, help="Minimum subjects for HIGH confidence. Default: 10.")
    confidence.add_argument("--high-support-taxa", default=3, type=positive_int, help="Minimum taxa for HIGH confidence. Default: 3.")
    confidence.add_argument("--min-output-confidence", default="LOW", choices=["LOW", "MEDIUM", "HIGH"], help="Minimum confidence written to FASTA/BED outputs. Default: LOW.")

    runtime = parser.add_argument_group("Runtime and logging")
    runtime.add_argument("--threads", default=4, type=positive_int, help="Worker threads. BLASTN also uses this value. Default: 4.")
    runtime.add_argument("--verbose", action="store_true", help="Print debug logs. Default: disabled.")

    return parser


def validate_args(args: argparse.Namespace) -> None:
    """Validate input arguments."""

    if not args.query.exists():
        raise FileNotFoundError(f"Query FASTA not found: {args.query}")
    supplied_modes = sum(x is not None for x in [args.blast, args.db, args.ref_fasta])
    if supplied_modes != 1:
        raise ValueError("Provide exactly one of --blast, --db, or --ref-fasta")
    if args.blast is not None and not args.blast.exists():
        raise FileNotFoundError(f"BLAST table not found: {args.blast}")
    if args.ref_fasta is not None and not args.ref_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {args.ref_fasta}")
    if args.taxonomy is not None and not args.taxonomy.exists():
        raise FileNotFoundError(f"Taxonomy table not found: {args.taxonomy}")
    if args.min_intron_len > args.max_intron_len:
        raise ValueError("--min-intron-len must be <= --max-intron-len")


def main() -> int:
    """CLI entry point."""

    parser = build_parser()
    args = parser.parse_args()
    setup_logging(args.verbose)
    validate_args(args)
    setup_output_paths(args)

    progress = Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=CONSOLE,
    )

    with progress:
        task = progress.add_task("Reading query FASTA", total=None)
        seqs = read_fasta(args.query)
        progress.update(task, completed=1, total=1)
        LOG.info("Loaded query sequences: %s", len(seqs))

        if args.ref_fasta is not None:
            task = progress.add_task("Preprocessing reference and building BLAST DB", total=None)
            generated_ref, generated_taxonomy, generated_db = preprocess_reference(args)
            args.taxonomy = generated_taxonomy
            args.db = generated_db
            progress.update(task, completed=1, total=1)

            if args.clean_ref_introns:
                task = progress.add_task("Cleaning reference introns by self-BLAST", total=None)
                _, cleaned_taxonomy, cleaned_db = clean_reference_introns(args, generated_ref, generated_taxonomy, generated_db)
                args.taxonomy = cleaned_taxonomy
                args.db = cleaned_db
                progress.update(task, completed=1, total=1)

        task = progress.add_task("Reading taxonomy", total=None)
        taxonomy = parse_taxonomy(args.taxonomy)
        progress.update(task, completed=1, total=1)

        task = progress.add_task("Preparing BLAST HSPs", total=None)
        blast_path = args.blast if args.blast is not None else run_query_blastn(args)
        args.blast = blast_path
        blast_by_query = parse_blast(blast_path, args.min_pident, args.min_hsp_len)
        progress.update(task, completed=1, total=1)

        query_ids = sorted(seqs.keys())
        task = progress.add_task("Detecting introns", total=len(query_ids))

        def worker(query_id: str) -> QueryResult:
            return analyze_query(query_id, blast_by_query.get(query_id, {}), len(seqs.get(query_id, "")), taxonomy, args)

        results: List[QueryResult] = []
        if args.threads == 1:
            for query_id in query_ids:
                results.append(worker(query_id))
                progress.advance(task)
        else:
            with futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
                for result in executor.map(worker, query_ids, chunksize=128):
                    results.append(result)
                    progress.advance(task)

        results.sort(key=lambda r: (confidence_rank(r.confidence), r.support_subjects, r.support_taxa, r.query_id), reverse=True)

        task = progress.add_task("Writing outputs", total=6)
        write_summary(args.summary_tsv, results, args.tax_rank)
        progress.advance(task)
        write_supporting_hsps(args.supporting_tsv, results)
        progress.advance(task)
        write_fasta_outputs(args, results, seqs)
        progress.advance(task)
        write_bed_outputs(args, results)
        progress.advance(task)
        write_report(args.report_md, results, args)
        progress.advance(task)
        print_summary(results, args)
        progress.advance(task)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
