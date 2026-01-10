#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

def filter_and_dedup(df: pd.DataFrame, max_evalue: float) -> pd.DataFrame:
    # compute the midpoint of aligned region on the target
    df['midpoint'] = (df['tstart'] + df['tend']) / 2

    # filter: midpoint must be in the last 50% of the target, and evalue <= threshold
    df_filtered = df[
        (df['midpoint'] > df['tlen'] / 2) &
        (df['evalue'] <= max_evalue)
    ].copy()

    # sort by bits in descending order
    df_sorted = df_filtered.sort_values('bits', ascending=False)

    # deduplicate by target, keeping the first (highest bits)
    df_dedup = df_sorted.drop_duplicates(subset=['target'], keep='first')

    # remove temporary column and restore original structure
    return df_dedup.drop(columns=['midpoint'])

def main():
    parser = argparse.ArgumentParser(
        description='Filter and deduplicate alignment results (TSV format)'
    )
    parser.add_argument(
        '-i', '--input', required=True,
        help='Input TSV file path (must include header)'
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Output TSV file path'
    )
    parser.add_argument(
        '-e', '--evalue', type=float, default=0.05,
        help='Maximum allowed evalue (default 0.05)'
    )
    args = parser.parse_args()

    # 1. read TSV
    df = pd.read_csv(
        args.input, sep='\t', comment='#', dtype={
            'query': str, 'target': str,
            'fident': float, 'alnlen': int, 'nident': int,
            'mismatch': int, 'gapopen': int,
            'qstart': int, 'qend': int,
            'tstart': int, 'tend': int,
            'evalue': float, 'qtmscore': float,
            'bits': float, 'tlen': int,
            'qaln': str, 'taln': str, 'tseq': str
        }
    )

    # 2–5. filter (with evalue), sort, deduplicate
    df_result = filter_and_dedup(df, args.evalue)

    # 6. write TSV without index
    df_result.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()
