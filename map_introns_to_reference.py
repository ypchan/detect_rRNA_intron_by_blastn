#!/usr/bin/env python3
"""
Map detected rRNA intron positions onto a standard full-length reference.

Inputs are designed to work directly with detect_rRNA_intron_by_blastn.py:
  1. *.summary.tsv
  2. *.intron_free.fa
  3. one standard full-length reference FASTA

The script aligns each intron-free sequence to the standard reference, maps the
intron insertion boundary through that alignment, and writes reference-relative
coordinates for each detected intron.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, TextIO, Tuple


CONFIDENCE_RANK = {"NONE": 0, "LOW": 1, "MEDIUM": 2, "HIGH": 3}


@dataclass(frozen=True)
class FastaRecord:
    seq_id: str
    header: str
    seq: str


@dataclass(frozen=True)
class IntronCall:
    query_id: str
    classification: str
    confidence: str
    intron_start: int
    intron_end: int
    intron_len: int
    raw_row: Dict[str, str]


@dataclass(frozen=True)
class AlignmentResult:
    strand: str
    score: int
    identity_pct: float
    matches: int
    aligned_pairs: int
    ref_start: int
    ref_end: int
    q_to_ref: List[Optional[int]]


@dataclass(frozen=True)
class BoundaryMapping:
    ref_left: Optional[int]
    ref_right: Optional[int]
    site_type: str
    coordinate: str
    bed_start: Optional[int]
    bed_end: Optional[int]


def open_text_auto(path: Path) -> TextIO:
    """Open plain text or gzip-compressed text."""

    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def clean_sequence(seq_parts: Iterable[str]) -> str:
    """Normalize FASTA sequence text for DNA-style alignment."""

    seq = "".join(seq_parts).upper().replace("U", "T")
    return "".join(seq.split())


def read_fasta(path: Path) -> List[FastaRecord]:
    """Read FASTA or FASTA.gz records."""

    records: List[FastaRecord] = []
    current_header: Optional[str] = None
    current_seq: List[str] = []

    def flush_record() -> None:
        nonlocal current_header, current_seq
        if current_header is None:
            return
        seq_id = current_header.split()[0]
        records.append(FastaRecord(seq_id=seq_id, header=current_header, seq=clean_sequence(current_seq)))
        current_header = None
        current_seq = []

    with open_text_auto(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush_record()
                current_header = line[1:]
            else:
                current_seq.append(line)
        flush_record()

    return records


def intron_free_query_id(record: FastaRecord) -> str:
    """Recover the original query ID from an intron-free FASTA header."""

    marker = "|intron_free|"
    if marker in record.seq_id:
        return record.seq_id.split(marker, 1)[0]
    if marker in record.header:
        return record.header.split(marker, 1)[0].split()[0]
    return record.seq_id


def parse_int_field(row: Dict[str, str], field: str) -> Optional[int]:
    """Parse a positive integer field, returning None for empty values."""

    value = row.get(field, "").strip()
    if not value:
        return None
    return int(value)


def read_intron_calls(path: Path, min_confidence: str) -> List[IntronCall]:
    """Read detector summary TSV and keep intron calls at or above threshold."""

    calls: List[IntronCall] = []
    min_rank = CONFIDENCE_RANK[min_confidence]
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            confidence = row.get("confidence", "NONE").strip() or "NONE"
            if CONFIDENCE_RANK.get(confidence, 0) < min_rank:
                continue
            intron_start = parse_int_field(row, "intron_start")
            intron_end = parse_int_field(row, "intron_end")
            if intron_start is None or intron_end is None:
                continue
            intron_len = parse_int_field(row, "intron_len") or (intron_end - intron_start + 1)
            calls.append(
                IntronCall(
                    query_id=row.get("query_id", "").strip(),
                    classification=row.get("classification", "").strip(),
                    confidence=confidence,
                    intron_start=intron_start,
                    intron_end=intron_end,
                    intron_len=intron_len,
                    raw_row=row,
                )
            )
    return calls


def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence with common ambiguity codes."""

    table = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
    return seq.translate(table)[::-1].upper()


def fitting_align(
    query: str,
    reference: str,
    strand: str,
    match_score: int,
    mismatch_score: int,
    gap_score: int,
) -> AlignmentResult:
    """
    Align the full query to the best-fitting reference interval.

    The full query is penalized end-to-end, while unaligned leading/trailing
    reference bases are free. This maps partial or full-length intron-free rRNA
    sequences into the coordinate system of the standard reference.
    """

    n = len(query)
    m = len(reference)
    if n == 0 or m == 0:
        raise ValueError("Query and reference sequences must be non-empty")

    prev = [0] * (m + 1)
    trace: List[bytearray] = [bytearray(m + 1)]

    for i in range(1, n + 1):
        curr = [i * gap_score] + [0] * m
        trace_row = bytearray(m + 1)
        trace_row[0] = 1
        qbase = query[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (match_score if qbase == reference[j - 1] else mismatch_score)
            up = prev[j] + gap_score
            left = curr[j - 1] + gap_score
            if diag >= up and diag >= left:
                curr[j] = diag
                trace_row[j] = 0
            elif up >= left:
                curr[j] = up
                trace_row[j] = 1
            else:
                curr[j] = left
                trace_row[j] = 2
        prev = curr
        trace.append(trace_row)

    end_j = max(range(m + 1), key=lambda j: prev[j])
    best_score = prev[end_j]
    i = n
    j = end_j
    aligned_query: List[str] = []
    aligned_ref: List[str] = []

    while i > 0:
        if j == 0:
            move = 1
        else:
            move = trace[i][j]
        if move == 0:
            aligned_query.append(query[i - 1])
            aligned_ref.append(reference[j - 1])
            i -= 1
            j -= 1
        elif move == 1:
            aligned_query.append(query[i - 1])
            aligned_ref.append("-")
            i -= 1
        else:
            aligned_query.append("-")
            aligned_ref.append(reference[j - 1])
            j -= 1

    ref_start = j + 1
    ref_end = end_j
    aligned_query.reverse()
    aligned_ref.reverse()

    q_to_ref: List[Optional[int]] = [None] * (n + 1)
    q_pos = 0
    r_pos = ref_start - 1
    matches = 0
    aligned_pairs = 0
    for qchar, rchar in zip(aligned_query, aligned_ref):
        if rchar != "-":
            r_pos += 1
        if qchar != "-":
            q_pos += 1
            if rchar != "-":
                q_to_ref[q_pos] = r_pos
                aligned_pairs += 1
                if qchar == rchar:
                    matches += 1

    identity_pct = 100.0 * matches / aligned_pairs if aligned_pairs else 0.0
    return AlignmentResult(
        strand=strand,
        score=best_score,
        identity_pct=identity_pct,
        matches=matches,
        aligned_pairs=aligned_pairs,
        ref_start=ref_start,
        ref_end=ref_end,
        q_to_ref=q_to_ref,
    )


def best_alignment(
    query: str,
    reference: str,
    strand_mode: str,
    match_score: int,
    mismatch_score: int,
    gap_score: int,
) -> Tuple[AlignmentResult, str]:
    """Align plus strand or both strands and return the best result."""

    candidates: List[Tuple[AlignmentResult, str]] = [
        (fitting_align(query, reference, "+", match_score, mismatch_score, gap_score), query)
    ]
    if strand_mode == "both":
        rc_query = reverse_complement(query)
        candidates.append((fitting_align(rc_query, reference, "-", match_score, mismatch_score, gap_score), rc_query))
    return max(candidates, key=lambda item: (item[0].score, item[0].identity_pct))


def nearest_mapped_left(q_to_ref: List[Optional[int]], after_query_pos: int) -> Optional[int]:
    """Find nearest mapped reference position at or before a query boundary."""

    for pos in range(after_query_pos, 0, -1):
        ref_pos = q_to_ref[pos]
        if ref_pos is not None:
            return ref_pos
    return None


def nearest_mapped_right(q_to_ref: List[Optional[int]], after_query_pos: int) -> Optional[int]:
    """Find nearest mapped reference position after a query boundary."""

    for pos in range(after_query_pos + 1, len(q_to_ref)):
        ref_pos = q_to_ref[pos]
        if ref_pos is not None:
            return ref_pos
    return None


def map_query_boundary_to_reference(q_to_ref: List[Optional[int]], after_query_pos: int, reference_id: str) -> BoundaryMapping:
    """Project a query boundary onto reference flanking coordinates."""

    ref_left = nearest_mapped_left(q_to_ref, after_query_pos)
    ref_right = nearest_mapped_right(q_to_ref, after_query_pos)

    if ref_left is None and ref_right is None:
        return BoundaryMapping(None, None, "unmapped", "unmapped", None, None)
    if ref_left is None:
        return BoundaryMapping(ref_left, ref_right, "before_reference_position", f"{reference_id}:before:{ref_right}", ref_right - 1, ref_right - 1)
    if ref_right is None:
        return BoundaryMapping(ref_left, ref_right, "after_reference_position", f"{reference_id}:after:{ref_left}", ref_left, ref_left)
    if ref_right == ref_left + 1:
        return BoundaryMapping(ref_left, ref_right, "between_adjacent_reference_bases", f"{reference_id}:between:{ref_left}-{ref_right}", ref_left, ref_left)
    if ref_right > ref_left + 1:
        return BoundaryMapping(ref_left, ref_right, "between_nonadjacent_reference_bases", f"{reference_id}:span:{ref_left + 1}-{ref_right - 1}", ref_left, ref_right - 1)
    return BoundaryMapping(ref_left, ref_right, "complex_mapping", f"{reference_id}:complex:{ref_left}-{ref_right}", min(ref_left, ref_right) - 1, max(ref_left, ref_right))


def choose_reference(records: List[FastaRecord], reference_id: Optional[str]) -> FastaRecord:
    """Pick the requested reference record, or the first record if unambiguous."""

    if not records:
        raise ValueError("Reference FASTA contains no records")
    if reference_id is None:
        if len(records) > 1:
            raise ValueError("Reference FASTA has multiple records; provide --reference-id")
        return records[0]
    for record in records:
        if record.seq_id == reference_id:
            return record
    available = ", ".join(record.seq_id for record in records[:10])
    raise ValueError(f"Reference ID not found: {reference_id}. Available examples: {available}")


def default_out_tsv(summary_path: Path) -> Path:
    """Return default TSV path next to the input summary."""

    name = summary_path.name
    suffix = ".summary.tsv"
    if name.endswith(suffix):
        return summary_path.with_name(name[: -len(suffix)] + ".reference_introns.tsv")
    return Path(str(summary_path) + ".reference_introns.tsv")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Map detected rRNA introns from query coordinates onto a standard full-length reference sequence.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--summary", required=True, type=Path, help="Summary TSV from detect_rRNA_intron_by_blastn.py.")
    parser.add_argument("--intron-free-fasta", required=True, type=Path, help="Intron-free FASTA from detect_rRNA_intron_by_blastn.py.")
    parser.add_argument("--reference", required=True, type=Path, help="Standard full-length reference FASTA/FASTA.gz.")
    parser.add_argument("--reference-id", default=None, help="Reference sequence ID to use when --reference contains multiple records.")
    parser.add_argument("--out-tsv", default=None, type=Path, help="Output TSV path.")
    parser.add_argument("--out-bed", default=None, type=Path, help="Optional BED output path with reference-relative insertion intervals.")
    parser.add_argument("--min-confidence", default="LOW", choices=["LOW", "MEDIUM", "HIGH"], help="Minimum intron confidence to map.")
    parser.add_argument("--strand", default="both", choices=["plus", "both"], help="Align only plus strand or choose the better of plus/reverse-complement.")
    parser.add_argument("--min-alignment-identity", default=0.0, type=float, help="Flag mappings below this percent identity as low identity.")
    parser.add_argument("--match-score", default=2, type=int, help="Alignment match score.")
    parser.add_argument("--mismatch-score", default=-3, type=int, help="Alignment mismatch score.")
    parser.add_argument("--gap-score", default=-5, type=int, help="Linear gap score.")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    for path_name in ["summary", "intron_free_fasta", "reference"]:
        path = getattr(args, path_name)
        if not path.exists():
            raise FileNotFoundError(f"Input file not found: {path}")
    if args.min_alignment_identity < 0 or args.min_alignment_identity > 100:
        raise ValueError("--min-alignment-identity must be between 0 and 100")


def output_fields() -> List[str]:
    return [
        "query_id",
        "reference_id",
        "status",
        "classification",
        "confidence",
        "query_intron_start",
        "query_intron_end",
        "query_intron_len",
        "intron_free_len",
        "intron_free_insertion_after_pos_original_orientation",
        "aligned_intron_free_insertion_after_pos",
        "strand",
        "ref_left_pos_1based",
        "ref_right_pos_1based",
        "reference_coordinate",
        "site_type",
        "bed_start_0based",
        "bed_end_0based",
        "alignment_score",
        "alignment_identity_pct",
        "alignment_matches",
        "alignment_aligned_pairs",
        "aligned_reference_start_1based",
        "aligned_reference_end_1based",
    ]


def blank_mapping_row(call: IntronCall, reference_id: str, status: str) -> Dict[str, str]:
    row = {field: "" for field in output_fields()}
    row.update(
        {
            "query_id": call.query_id,
            "reference_id": reference_id,
            "status": status,
            "classification": call.classification,
            "confidence": call.confidence,
            "query_intron_start": str(call.intron_start),
            "query_intron_end": str(call.intron_end),
            "query_intron_len": str(call.intron_len),
        }
    )
    return row


def map_call(
    call: IntronCall,
    intron_free_records: Dict[str, FastaRecord],
    reference: FastaRecord,
    args: argparse.Namespace,
) -> Dict[str, str]:
    """Map one intron call to the standard reference."""

    record = intron_free_records.get(call.query_id)
    if record is None:
        return blank_mapping_row(call, reference.seq_id, "missing_intron_free_sequence")

    intron_free_len = len(record.seq)
    insertion_after_original = call.intron_start - 1
    if insertion_after_original < 0 or insertion_after_original > intron_free_len:
        return blank_mapping_row(call, reference.seq_id, "invalid_intron_boundary_for_intron_free_sequence")

    alignment, aligned_query = best_alignment(
        query=record.seq,
        reference=reference.seq,
        strand_mode=args.strand,
        match_score=args.match_score,
        mismatch_score=args.mismatch_score,
        gap_score=args.gap_score,
    )

    if alignment.strand == "+":
        insertion_after_aligned = insertion_after_original
    else:
        insertion_after_aligned = intron_free_len - insertion_after_original

    boundary = map_query_boundary_to_reference(alignment.q_to_ref, insertion_after_aligned, reference.seq_id)
    status = "mapped"
    if boundary.site_type == "unmapped":
        status = "unmapped_boundary"
    elif alignment.identity_pct < args.min_alignment_identity:
        status = "mapped_low_identity"

    return {
        "query_id": call.query_id,
        "reference_id": reference.seq_id,
        "status": status,
        "classification": call.classification,
        "confidence": call.confidence,
        "query_intron_start": str(call.intron_start),
        "query_intron_end": str(call.intron_end),
        "query_intron_len": str(call.intron_len),
        "intron_free_len": str(intron_free_len),
        "intron_free_insertion_after_pos_original_orientation": str(insertion_after_original),
        "aligned_intron_free_insertion_after_pos": str(insertion_after_aligned),
        "strand": alignment.strand,
        "ref_left_pos_1based": "" if boundary.ref_left is None else str(boundary.ref_left),
        "ref_right_pos_1based": "" if boundary.ref_right is None else str(boundary.ref_right),
        "reference_coordinate": boundary.coordinate,
        "site_type": boundary.site_type,
        "bed_start_0based": "" if boundary.bed_start is None else str(boundary.bed_start),
        "bed_end_0based": "" if boundary.bed_end is None else str(boundary.bed_end),
        "alignment_score": str(alignment.score),
        "alignment_identity_pct": f"{alignment.identity_pct:.2f}",
        "alignment_matches": str(alignment.matches),
        "alignment_aligned_pairs": str(alignment.aligned_pairs),
        "aligned_reference_start_1based": str(alignment.ref_start),
        "aligned_reference_end_1based": str(alignment.ref_end),
    }


def write_bed(path: Path, rows: List[Dict[str, str]]) -> None:
    """Write reference-relative BED intervals for mapped intron positions."""

    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in rows:
            if row.get("status") not in {"mapped", "mapped_low_identity"}:
                continue
            if not row.get("bed_start_0based") or not row.get("bed_end_0based"):
                continue
            name = (
                f"{row['query_id']}|intron={row['query_intron_start']}-{row['query_intron_end']}"
                f"|ref={row['reference_coordinate']}|confidence={row['confidence']}"
            )
            writer.writerow([row["reference_id"], row["bed_start_0based"], row["bed_end_0based"], name, ".", row["strand"]])


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args)

    out_tsv = args.out_tsv or default_out_tsv(args.summary)
    reference = choose_reference(read_fasta(args.reference), args.reference_id)
    intron_free_records = {intron_free_query_id(record): record for record in read_fasta(args.intron_free_fasta)}
    calls = read_intron_calls(args.summary, args.min_confidence)

    rows = [map_call(call, intron_free_records, reference, args) for call in calls]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    if args.out_bed is not None:
        args.out_bed.parent.mkdir(parents=True, exist_ok=True)
        write_bed(args.out_bed, rows)

    mapped = sum(1 for row in rows if row["status"] in {"mapped", "mapped_low_identity"})
    print(f"Mapped introns: {mapped}/{len(rows)}")
    print(f"Output TSV: {out_tsv}")
    if args.out_bed is not None:
        print(f"Output BED: {args.out_bed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

