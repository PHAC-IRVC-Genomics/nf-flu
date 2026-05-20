#!/usr/bin/env python3
"""
dvg_filter.py — Streaming SAM filter for influenza Defective Viral Genome (DVG) exclusion.

Reads SAM format from stdin, applies up to four independent filter criteria to each
alignment record, and writes passing records to stdout. Header lines are always passed
through unchanged. Filter statistics are written to a JSON file on completion.

If no reads pass filtering, a warning is emitted and only the SAM header is written
to stdout (producing a valid but empty BAM when piped through samtools view -b).

Filter criteria (applied sequentially; first failure exits early):
    1. Minimum aligned length   — mapping confidence floor (M/=/X CIGAR ops)
    2. Maximum single N/D op    — primary DVG deletion signature filter
    3. Single alignment block   — split-alignment / DVG structural guard
    4. Maximum soft-clip frac   — internally-mapping DVG read filter

Notable reads — passing reads flagged for inspection:
    Reads that pass all filter criteria but contain a single D/N operation in the
    range [NOTABLE_DELETION_MIN, max_skip) are flagged as 'notable_deletion'.
    These represent biologically ambiguous events — larger than typical nanopore
    indel noise but below the DVG threshold — and are recorded in the JSON stats
    for downstream inspection without being removed from the filtered BAM.

Usage:
    samtools view -h input.bam \\
        | dvg_filter.py \\
            --min-len 890 \\
            --max-skip 350 \\
            --max-clip-frac 0.15 \\
            --segment NS \\
            --sample sample_01 \\
            --stats-out stats.json \\
        | samtools view -b | samtools sort > filtered.bam

References:
    Frensing (2015): Defective interfering viruses and their impact on disease and treatment.
    Felt & Vignuzzi (mBio, 2022): DelVG definition and packaging signal retention.
    Romero-Brey et al. (PMC, 2022): 50-150 nt terminal retention across segments.
    Genome Medicine (2025): Nanopore-specific DVG sizes 400-500 nt for polymerase segments.
"""

import argparse
import json
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional


# Minimum deletion size to consider flagging as notable.
# Below this is consistent with normal nanopore indel noise.
NOTABLE_DELETION_MIN = 50


# ── Data structures ────────────────────────────────────────────────────────────

@dataclass
class CigarMetrics:
    """Metrics derived from parsing a single CIGAR string."""
    aligned_len: int   # Total M/=/X bases — the aligned query length
    max_deletion: int  # Largest single N or D operation
    soft_clip: int     # Total S and H bases
    num_blocks: int    # Number of contiguous alignment blocks


@dataclass
class FilterResult:
    """Outcome of applying all filter criteria to one read."""
    passed: bool
    fail_reason: Optional[str] = None  # None if passed; criterion name if failed


@dataclass
class FilterStats:
    """Accumulated counts across all records processed."""
    segment: str = ""
    sample: str = ""
    thresholds: dict = field(default_factory=dict)
    total: int = 0
    passed: int = 0
    fail_minlen: int = 0
    fail_maxskip: int = 0
    fail_blocks: int = 0
    fail_clipfrac: int = 0
    notable_reads: List[Dict] = field(default_factory=list)
    warning: Optional[str] = None

    @property
    def failed(self) -> int:
        return self.total - self.passed

    def to_dict(self) -> dict:
        d = {
            "sample": self.sample,
            "segment": self.segment,
            "thresholds": self.thresholds,
            "counts": {
                "total": self.total,
                "pass": self.passed,
                "fail_total": self.failed,
                "fail_minlen": self.fail_minlen,
                "fail_maxskip": self.fail_maxskip,
                "fail_blocks": self.fail_blocks,
                "fail_clipfrac": self.fail_clipfrac,
            },
            "notable_reads": self.notable_reads,
        }
        if self.warning:
            d["warning"] = self.warning
        return d


# ── CIGAR parsing ──────────────────────────────────────────────────────────────

def parse_cigar(cigar: str, max_skip: int) -> Optional[CigarMetrics]:
    """
    Parse a CIGAR string and return alignment metrics.

    Args:
        cigar:    CIGAR string from SAM field 6. '*' returns None.
        max_skip: Gap size at or above which a D/N op is considered block-breaking.
                  Pass 0 to disable block tracking entirely.

    Returns:
        CigarMetrics, or None if cigar is '*'.

    Notes:
        - M, =, X ops accumulate into aligned_len and track block continuity.
        - S and H ops accumulate into soft_clip (query length but not alignment).
        - N and D ops: largest value tracked as max_deletion; ops >= max_skip
          break the current alignment block (in_block reset to 0).
        - I and P ops are intentionally ignored — they do not affect reference
          coordinates or block continuity.
    """
    if cigar == '*':
        return None

    aligned_len = 0
    max_deletion = 0
    soft_clip = 0
    num_blocks = 0
    in_block = False
    num_str = []

    for ch in cigar:
        if ch.isdigit():
            num_str.append(ch)
        else:
            if num_str:
                num_val = int(''.join(num_str))
                num_str = []

                if ch in ('M', '=', 'X'):
                    aligned_len += num_val
                    if not in_block:
                        num_blocks += 1
                        in_block = True

                elif ch in ('S', 'H'):
                    soft_clip += num_val

                elif ch in ('N', 'D'):
                    if num_val > max_deletion:
                        max_deletion = num_val
                    # A gap at or above max_skip breaks alignment block continuity.
                    # Threshold is unified with max_skip so block-break and deletion-size
                    # filters are internally consistent.
                    if max_skip > 0 and num_val >= max_skip:
                        in_block = False

                # I (insertion) and P (padding): intentionally ignored.

    return CigarMetrics(
        aligned_len=aligned_len,
        max_deletion=max_deletion,
        soft_clip=soft_clip,
        num_blocks=num_blocks,
    )


# ── Filter logic ───────────────────────────────────────────────────────────────

def filter_read(
    metrics: CigarMetrics,
    min_len: int,
    max_skip: int,
    max_clip_frac: float,
) -> FilterResult:
    """
    Apply all filter criteria to a read's CIGAR metrics.

    Criteria are applied sequentially; the first failure exits early.
    This means counts reflect reads that *reached* each filter, not
    the total read pool — documented in stats output accordingly.

    Args:
        metrics:       Parsed CIGAR metrics for this read.
        min_len:       Minimum aligned bases required (0 = disabled).
        max_skip:      Maximum allowed single N/D op in bp (0 = disabled).
        max_clip_frac: Maximum allowed soft/hard clip fraction (0.0 = disabled).

    Returns:
        FilterResult with passed=True or passed=False and the failing criterion name.
    """
    # 1. Minimum aligned length — mapping confidence floor
    if min_len > 0 and metrics.aligned_len < min_len:
        return FilterResult(passed=False, fail_reason='fail_minlen')

    # 2. Maximum single N/D operation — primary DVG deletion signature
    if max_skip > 0 and metrics.max_deletion >= max_skip:
        return FilterResult(passed=False, fail_reason='fail_maxskip')

    # 3. Single contiguous alignment block — split-alignment / DVG structural guard.
    #    Only evaluated when max_skip > 0; when max_skip == 0 the block-break
    #    threshold is undefined so this check is meaningless.
    if max_skip > 0 and metrics.num_blocks != 1:
        return FilterResult(passed=False, fail_reason='fail_blocks')

    # 4. Soft-clip fraction — internally-mapping DVG read filter.
    #    Legitimate clipping should only represent primer overhangs (~20-30 bp/end).
    if max_clip_frac > 0.0:
        total_query = metrics.aligned_len + metrics.soft_clip
        clip_frac = metrics.soft_clip / total_query if total_query > 0 else 0.0
        if clip_frac > max_clip_frac:
            return FilterResult(passed=False, fail_reason='fail_clipfrac')

    return FilterResult(passed=True)


def is_notable(metrics: CigarMetrics, max_skip: int) -> Optional[str]:
    """
    Determine whether a passing read should be flagged as notable.

    A read is notable if it contains a single D/N operation in the range
    [NOTABLE_DELETION_MIN, max_skip) — larger than typical nanopore indel
    noise but below the DVG threshold. These are biologically ambiguous
    events worth recording for downstream inspection.

    Args:
        metrics:  Parsed CIGAR metrics for this read.
        max_skip: The active DVG deletion threshold. If 0, notability
                  cannot be evaluated and None is returned.

    Returns:
        A flag string if notable, None otherwise.
    """
    if max_skip > 0 and NOTABLE_DELETION_MIN <= metrics.max_deletion < max_skip:
        return 'notable_deletion'
    return None


# ── SAM streaming ──────────────────────────────────────────────────────────────

def process_sam(
    in_stream,
    out_stream,
    min_len: int,
    max_skip: int,
    max_clip_frac: float,
    segment: str,
    sample: str,
    stats_out: str,
) -> FilterStats:
    """
    Stream SAM records from in_stream, filter, and write passing records to out_stream.

    Header lines (@) are always passed through unchanged and not counted.
    Records with CIGAR '*' are silently skipped (guards against unmapped reads
    that slipped through upstream -F filtering).

    Passing reads with a deletion in [NOTABLE_DELETION_MIN, max_skip) are flagged
    as notable and recorded in the JSON stats for downstream inspection.

    Args:
        in_stream:     Readable stream of SAM-format text (typically sys.stdin).
        out_stream:    Writable stream for filtered SAM output (typically sys.stdout).
        min_len:       Minimum aligned bases threshold (0 = disabled).
        max_skip:      Maximum N/D op threshold in bp (0 = disabled).
        max_clip_frac: Maximum clip fraction threshold (0.0 = disabled).
        segment:       Segment name for logging/stats context only.
        sample:        Sample name for logging/stats context only.
        stats_out:     Path to write JSON stats file.

    Returns:
        Populated FilterStats instance.
    """
    stats = FilterStats(
        segment=segment,
        sample=sample,
        thresholds={
            "min_len": min_len,
            "max_skip": max_skip,
            "max_clip_frac": max_clip_frac,
            "notable_deletion_min": NOTABLE_DELETION_MIN,
        }
    )

    for line in in_stream:
        line = line.rstrip('\n')

        # Pass SAM header lines through unchanged
        if line.startswith('@'):
            out_stream.write(line + '\n')
            continue

        fields = line.split('\t')
        if len(fields) < 6:
            continue

        qname = fields[0]
        cigar = fields[5]

        # Skip reads with no CIGAR
        metrics = parse_cigar(cigar, max_skip)
        if metrics is None:
            continue

        stats.total += 1
        result = filter_read(metrics, min_len, max_skip, max_clip_frac)

        if result.passed:
            stats.passed += 1
            out_stream.write(line + '\n')

            # Check whether this passing read should be flagged as notable
            flag = is_notable(metrics, max_skip)
            if flag:
                total_query = metrics.aligned_len + metrics.soft_clip
                clip_frac = (
                    metrics.soft_clip / total_query if total_query > 0 else 0.0
                )
                stats.notable_reads.append({
                    "qname": qname,
                    "cigar": cigar,
                    "aligned_len": metrics.aligned_len,
                    "max_deletion": metrics.max_deletion,
                    "clip_frac": round(clip_frac, 4),
                    "flag": flag,
                })
        else:
            # Increment the appropriate per-criterion counter
            setattr(stats, result.fail_reason, getattr(stats, result.fail_reason) + 1)

    # ── Post-processing ────────────────────────────────────────────────────────
    if stats.total > 0 and stats.passed == 0:
        msg = (
            f"WARNING: No reads passed DVG filtering for "
            f"sample={sample}, segment={segment}. "
            f"Header-only BAM written. Check filter thresholds."
        )
        stats.warning = msg
        print(msg, file=sys.stderr)

    # Write JSON stats
    with open(stats_out, 'w') as fh:
        json.dump(stats.to_dict(), fh, indent=2)

    return stats


# ── Entry point ────────────────────────────────────────────────────────────────

def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--min-len',
        type=int,
        default=0,
        help='Minimum aligned bases (M/=/X) required to retain a read. 0 = disabled.',
    )
    parser.add_argument(
        '--max-skip',
        type=int,
        default=0,
        help=(
            'Maximum allowed single N or D CIGAR operation in bp. '
            'Reads at or above this threshold are discarded as likely DVGs. '
            '0 = disabled.'
        ),
    )
    parser.add_argument(
        '--max-clip-frac',
        type=float,
        default=0.0,
        help=(
            'Maximum allowed soft/hard clip fraction of total query length. '
            'Targets internally-mapping DVG reads with heavily clipped flanks. '
            '0.0 = disabled.'
        ),
    )
    parser.add_argument(
        '--segment',
        type=str,
        default='UNKNOWN',
        help='Influenza segment name (e.g. NS, PB2). Used for logging only.',
    )
    parser.add_argument(
        '--sample',
        type=str,
        default='UNKNOWN',
        help='Sample name. Used for logging only.',
    )
    parser.add_argument(
        '--stats-out',
        type=str,
        default='dvg_filter_stats.json',
        help='Path to write JSON filter statistics file.',
    )
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_args(argv)

    print(
        f"DVG filter: sample={args.sample}, segment={args.segment}, "
        f"min_len={args.min_len}, max_skip={args.max_skip}, "
        f"max_clip_frac={args.max_clip_frac}",
        file=sys.stderr,
    )

    stats = process_sam(
        in_stream=sys.stdin,
        out_stream=sys.stdout,
        min_len=args.min_len,
        max_skip=args.max_skip,
        max_clip_frac=args.max_clip_frac,
        segment=args.segment,
        sample=args.sample,
        stats_out=args.stats_out,
    )

    print(
        f"DVG filter complete: total={stats.total}, pass={stats.passed}, "
        f"fail_minlen={stats.fail_minlen}, fail_maxskip={stats.fail_maxskip}, "
        f"fail_blocks={stats.fail_blocks}, fail_clipfrac={stats.fail_clipfrac}, "
        f"notable={len(stats.notable_reads)}",
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
