#!/usr/bin/env python

import pandas as pd
import argparse
import numpy as np
import os

def main(args):

    # read sbarrier scores:
    scores = pd.read_csv(args.barrier_scores, sep=args.delimiter, encoding='utf-8')

    #define metrics for aggregation:
    metrics = ['all_paths_barrier_fraction', 'all_paths_adjacent_fraction']

    # group dataframe by image and source cell used in barrier calculation:
    grouped = scores.groupby([args.image_id_col, 'source_cell'])
    
    # define aggregation functions for summary table:
    agg_funcs = [np.mean, np.median, np.std]
    summary_data = grouped[metrics].agg(agg_funcs)

    # join multiindex columns to single:
    summary_data.columns = [' '.join(col).strip().replace(' ', '_') for col in summary_data.columns.values]
    summary_data = summary_data.reset_index()

    # save output:
    fname = os.path.split(args.barrier_scores)[1]
    fname = fname.replace('_barrier_scores.csv', '_barrier_summary.csv')
    fname = fname.replace(' ', '_')
    spath = os.path.join(args.outdir, fname)
    summary_data.to_csv(spath, index=False, sep='\t')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Aggregate barrier scores')
    parser.add_argument('--barrier_scores', help='Path to aggregate barrier scores for all images', required=True)
    parser.add_argument('--delimiter', help='Delimiter for barrier scores', default='\t')
    parser.add_argument('--outdir', help='Output directory', required=True, default='.')
    parser.add_argument('--image_id_col', help='Column name for image id', default='imagename')

    args = parser.parse_args()
    main(args)
