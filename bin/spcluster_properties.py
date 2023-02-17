#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns

def get_clustered(df, clustered_cellType):
    return df[df[f'{clustered_cellType}_spatial_cluster_id']>=0]

def get_clustered_cellType(df):
    clustering_col = [x for x in list(df) if '_spatial_cluster_id' in x][0]
    clustered_cellType = clustering_col.split('_spatial_cluster_id')[0]
    return clustered_cellType


def plot(df, phenotyping_column, imagename, clustered_cellType, outdir):
    
    fig, ax = plt.subplots(nrows=1,ncols=2, figsize=(10.5,10.5))
    sns.set_style('white')

    order = get_order(df, phenotyping_column, 'intracluster_density')
    g = sns.boxplot(data=df, x=phenotyping_column, y='intracluster_density', showcaps=False, order=order, palette='coolwarm', ax=ax[0])
    g = sns.stripplot(data=df, x=phenotyping_column, y='intracluster_density', order=order, color='k', s=5, ax=ax[0])

    ax[0].set_yscale('log')
    ax[0].tick_params(labelrotation=90)

    g = sns.boxplot(data=df, x=phenotyping_column, y='intracluster_fraction', showcaps=False, order=order, palette='coolwarm', ax=ax[1])
    g = sns.stripplot(data=df, x=phenotyping_column, y='intracluster_fraction', order=order, color='k', s=5, ax=ax[1])

    ax[1].set_yscale('log')
    ax[1].tick_params(labelrotation=90)
    plt.suptitle(imagename)
    plt.tight_layout()
    
    plt.savefig(os.path.join(outdir, f'{imagename}_{clustered_cellType}_{phenotyping_column}_intracluster_properties.png'))


def get_order(df, cat, val, ordering='median'):
    """
    Get ordering for categorical plots based on e.g. increasing mean or median.
    """
    grouped = df.groupby(cat)
    if ordering == 'median':
        print(grouped.median()[val].sort_values())
        order = grouped.median()[val].sort_values().index.tolist()   
    if ordering == 'mean':
        print(grouped.mean()[val].sort_values())
        order = grouped.mean()[val].sort_values().index.tolist() 
    return order


def main(args):

    # get args:
    
    data = pd.read_csv(args.clustered_data, sep=args.delimiter)
    phenotyping_column = args.phenotyping_column
    imagename = args.imagename
    outdir = args.outdir
    clustered_cellType = get_clustered_cellType(data)

    # calculate intracluster cdensity of cells:
    clustered = get_clustered(data, clustered_cellType)

    if len(clustered) > 0:
        grouped = clustered.groupby([f'{clustered_cellType}_spatial_cluster_id', f'{clustered_cellType}_cluster_area', phenotyping_column])
        grouped = grouped.agg({'imagename': 'count'}).rename(columns={'imagename':'cells_per_cluster'}).reset_index()
        grouped['intracluster_density'] = grouped['cells_per_cluster'] / grouped[f'{clustered_cellType}_cluster_area']

        # calculate total cells per cluster andcell type fractions per cluster:
        total_cells_per_cluster = grouped.groupby([f'{clustered_cellType}_spatial_cluster_id']).agg({'cells_per_cluster':'sum'}).rename(columns={'cells_per_cluster':'total_cells_per_cluster'}).reset_index()
        grouped = pd.merge(grouped, total_cells_per_cluster, on=f'{clustered_cellType}_spatial_cluster_id')
        grouped['intracluster_fraction'] = grouped['cells_per_cluster'] / grouped['total_cells_per_cluster']
        grouped['imagename'] = imagename
        grouped['clustered_cellType'] = clustered_cellType

        # save output:
        grouped.to_csv(os.path.join(outdir, f'{imagename}_{clustered_cellType}_{phenotyping_column}_intracluster_densities.csv'))
        plot(grouped, phenotyping_column, imagename, clustered_cellType, outdir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate cell specific properties (e.g. density) of spatial clusters.')
    parser.add_argument('--clustered_data', help='Path to clustered cell data.', required=True)
    parser.add_argument('--phenotyping_column', help='Phenotyping level', required=True)
    parser.add_argument('--imagename', help='ID of image', required=True)
    parser.add_argument('--outdir', help='Output directory', required=True, default='.')
    parser.add_argument('--delimiter', help='Delimiter of input file', required=False, default='\t')
    args = parser.parse_args()
    main(args)