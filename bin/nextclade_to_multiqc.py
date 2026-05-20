#!/usr/bin/env python3
"""
Convert Nextclade TSV output to MultiQC custom content format.

This script reads a combined Nextclade TSV file, deduplicates entries,
selects the best QC result per sample/segment, and outputs a MultiQC-compatible
TSV file with custom headers.

When multiple datasets produce 'good' QC results for the same sample/segment
(e.g., H5 samples matching h5n1, h5n1-cattle-outbreak, and h5nx), all 'good'
results are retained. If no 'good' results exist, the best available is used.

Use --split-segments to generate separate MultiQC sections for specific segments
(e.g., --split-segments HA,NA). Each segment gets its own table in the report.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

__version__ = "1.5.1"

QC_RANK = {'good': 0, 'mediocre': 1, 'bad': 2, '': 3}

# All columns (used for "all segments" output)
COLUMNS_ALL = [
    'sample', 'subtype', 'segment', 'legacy-clade', 'short-clade', 'clade',
    'cdsCoverage', 'totalNonACGTNs', 'nonACGTNs', 'qc.frameShifts.frameShifts',
    'qc.stopCodons.stopCodons', 'dataset_tag'
]

# Segment-specific columns (exclude 'segment' since it's redundant)
COLUMNS_SEGMENT = [
    'sample', 'subtype', 'legacy-clade', 'short-clade', 'clade',
    'cdsCoverage', 'totalNonACGTNs', 'nonACGTNs', 'qc.frameShifts.frameShifts',
    'qc.stopCodons.stopCodons', 'dataset_tag'
]

COLUMN_LABELS = {
    'sample':                      'Sample',
    'subtype':                     'Subtype',
    'segment':                     'Segment',
    'legacy-clade':                'Legacy Clade',
    'short-clade':                 'Short Clade',
    'clade':                       'Clade',
    'cdsCoverage':                 'cds Coverage',
    'totalNonACGTNs':              '# Mixed Sites',
    'nonACGTNs':                   'Mixed Sites',
    'qc.frameShifts.frameShifts':  'Frameshifts',
    'qc.stopCodons.stopCodons':    'Stop Codons',
    'dataset_tag':                 'Dataset Tag',
}

# Human-readable segment names for MultiQC section titles
SEGMENT_NAMES = {
    'HA':  'Hemagglutinin (HA)',
    'NA':  'Neuraminidase (NA)',
    'MP':  'Matrix Protein (MP)',
    'NP':  'Nucleoprotein (NP)',
    'NS':  'Non-Structural (NS)',
    'PA':  'Polymerase Acidic (PA)',
    'PB1': 'Polymerase Basic 1 (PB1)',
    'PB2': 'Polymerase Basic 2 (PB2)',
}


def parse_dataset(dataset_name):
    """Extract subtype and segment from dataset name."""
    parts = dataset_name.split('/')
    return parts[2].upper(), parts[3].upper()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Convert Nextclade TSV to MultiQC custom content format'
    )
    parser.add_argument(
        'input_tsv',
        help='Input Nextclade TSV file'
    )
    parser.add_argument(
        'output_tsv',
        help='Output MultiQC TSV file (used for "all segments" when --split-segments is specified)'
    )
    parser.add_argument(
        '--split-segments',
        type=str,
        default=None,
        metavar='SEGMENTS',
        help='Comma-separated list of segments to split into separate MultiQC sections '
             '(e.g., --split-segments HA,NA). Each segment gets its own file: '
             'nextclade_<segment>_mqc.tsv'
    )
    parser.add_argument(
        '--version',
        action='version',
        version=f'%(prog)s {__version__}'
    )
    return parser.parse_args(argv)


def process_rows(rows):
    """Process input rows and return deduplicated, filtered output rows."""
    # Deduplicate by sample/subtype/segment
    seen = set()
    deduplicated = []
    for r in rows:
        subtype, segment = parse_dataset(r['dataset_name'])
        key = (r['sample'], subtype, segment)
        if key not in seen:
            seen.add(key)
            deduplicated.append((subtype, segment, r))

    # Group by sample/segment
    segment_candidates = defaultdict(list)
    for subtype, segment, r in deduplicated:
        segment_candidates[(r['sample'], segment)].append((subtype, r))

    output_rows = []
    for (sample, segment), candidates in segment_candidates.items():
        # Filter for 'good' QC results
        good_candidates = [(subtype, r) for subtype, r in candidates if r['qc.overallStatus'] == 'good']

        if good_candidates:
            # Keep all 'good' results
            results_to_output = good_candidates
        else:
            # No 'good' results, fall back to best available
            best = min(candidates, key=lambda x: QC_RANK.get(x[1]['qc.overallStatus'], 3))
            results_to_output = [best]

        for subtype, r in results_to_output:
            # Create unique ID for MultiQC by combining sample, segment, and subtype
            sample_segment_id = f"{sample}_{segment}_{subtype}"
            output_rows.append({
                'id':                          sample_segment_id,
                'sample':                      sample,
                'subtype':                     subtype,
                'segment':                     segment,
                'legacy-clade':                r['legacy-clade'],
                'short-clade':                 r['short-clade'],
                'clade':                       r['clade'],
                'cdsCoverage':                 r['cdsCoverage'],
                'totalNonACGTNs':              r['totalNonACGTNs'],
                'nonACGTNs':                   r['nonACGTNs'],
                'qc.frameShifts.frameShifts':  r['qc.frameShifts.frameShifts'],
                'qc.stopCodons.stopCodons':    r['qc.stopCodons.stopCodons'],
                'dataset_tag':                 r['dataset_tag'],
            })

    output_rows.sort(key=lambda x: (x['sample'], x['segment'], x['subtype']))
    return output_rows


def write_multiqc_tsv(output_rows, output_path, section_id, section_name, description, columns):
    """Write output rows to a MultiQC-compatible TSV file."""
    with open(output_path, 'w', newline='') as f:
        f.write(f"# id: '{section_id}'\n")
        f.write(f"# section_name: '{section_name}'\n")
        f.write(f"# description: '{description}'\n")
        f.write("# format: 'tsv'\n")
        f.write("# plot_type: 'table'\n")
        f.write(f"# anchor: '{section_id}'\n")
        f.write("# headers:\n")
        for col in columns:
            f.write(f"#   {col}:\n")
            f.write(f"#     title: '{COLUMN_LABELS[col]}'\n")
        # Use 'id' as first column (MultiQC row identifier) but don't show it
        writer = csv.DictWriter(f, fieldnames=['id'] + columns, delimiter='\t',
                                extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(output_rows)


def main(argv=None):
    args = parse_args(argv)

    # Read input TSV
    with open(args.input_tsv) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    # Process all rows
    output_rows = process_rows(rows)

    if args.split_segments:
        # Parse comma-separated segment list and normalize to uppercase
        segments_to_split = [s.strip().upper() for s in args.split_segments.split(',')]
        
        # Determine output directory
        output_dir = os.path.dirname(args.output_tsv) or '.'

        # Write a separate file for each specified segment
        for segment in segments_to_split:
            segment_rows = [r for r in output_rows if r['segment'].upper() == segment]
            segment_lower = segment.lower()
            segment_output = os.path.join(output_dir, f'nextclade_{segment_lower}_mqc.tsv')
            
            # Get human-readable segment name, fallback to just the segment code
            segment_display = SEGMENT_NAMES.get(segment, segment)
            
            write_multiqc_tsv(
                segment_rows,
                segment_output,
                section_id=f'nextclade_{segment_lower}',
                section_name=f'Nextclade {segment_display}',
                description=f'Nextclade clade assignment and QC results for {segment_display} segment',
                columns=COLUMNS_SEGMENT
            )

        # Write all segments file
        write_multiqc_tsv(
            output_rows,
            args.output_tsv,
            section_id='nextclade_all',
            section_name='Nextclade All Segments',
            description='Nextclade clade assignment and QC results for all segments',
            columns=COLUMNS_ALL
        )
    else:
        # Single output file with all segments
        write_multiqc_tsv(
            output_rows,
            args.output_tsv,
            section_id='nextclade',
            section_name='Nextclade',
            description='Nextclade clade assignment and QC results per sample and segment',
            columns=COLUMNS_ALL
        )


if __name__ == '__main__':
    sys.exit(main())
