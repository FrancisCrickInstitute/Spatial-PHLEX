#!/usr/bin/env python

import pandas as pd
import numpy as np
from sklearn.cluster import MiniBatchKMeans
import argparse

import matplotlib
import matplotlib.pyplot as plt
import os

import seaborn as sns

def count_Type_occurrences(neighbourhood_window, type_dict):
    """Count the frequency of each cell type in a neighbourhood window

    Args:
        neighbourhood_window (np.ndarray): The N nearest neighbours of a cell
        type_dict (dict): Initialised dictionary of unique cell types in dataset.   

    Returns:
        np.ndarray: The frequency of each cell type in the neighbourhood window as numpy array
    """    
    
    for cell in neighbourhood_window:
        type_dict[cell] += 1
    return np.fromiter(type_dict.values(), dtype=int)

def cell_frequencies(data, cellTypes, type_cols):
    """Function to turn a dataframe of nearest neighbourhood windows into a dataframe of cell type frequencies

    Args:
        data (pd.DataFrame): Dataframe of nearest neighbourhood windows. Rows should correspond to every cell in the dataset. Columns need to contain the columns specified in type_cols.
        cellTypes (list): list of unique cell types in the dataset.
        type_cols (list): List of column headers in data that contain the cell type in formation for each cell\s nearest neighbour. e.g. ['cellType', 'cellType_Neighbour1', 'cellType_Neighbour2', 'cellType_Neighbour3', 'cellType_Neighbour4', 'cellType_Neighbour5']

    Returns:
        pd.DataFrame: DataFrame of cell frequencies. Rows correspond to every neighbourhood window in the dataset. Columns correspond to the unique cell types in the dataset.
    """


    data_values = data[type_cols].astype(str).values
    frequencies = np.zeros((data_values.shape[0], len(cellTypes)))
    for i in range(data_values.shape[0]):
        type_dict = dict(zip(cellTypes, np.zeros(len(cellTypes))))
        frequencies[i,:] =  count_Type_occurrences(data_values[i,:], type_dict)
        # frequencies[i,:] = occs/np.sum(occs)

    return pd.DataFrame(frequencies, columns=cellTypes)


def find_communities(frequencies, n_communities=10):
    """Function to cluster cell type frequencies into communities through the method of Schürch et al.

    Args:
        frequencies (): _description_

    Returns:
        _type_: _description_
    """    
    clusterer = MiniBatchKMeans(n_clusters=n_communities, random_state=10)
    community_clusters = clusterer.fit_predict(frequencies)
    centers = clusterer.cluster_centers_
    return community_clusters, centers



def community_plot(community_centers, image_data, cellTypes, communities_type):
    """Generates a spatial scatter plot of the community centers and a heatmap of the community centers.

    Args:
        community_centers (np.ndarray): Community centers identified with sklearn kmeans clustering.
        image_data (pd.DataFrame): Dataframe containing the single cell image data, with columns 'centerX', 'centerY' and 'community'.
    """    

    sns.set_style('white')
    g = sns.clustermap(community_centers, method='ward', metric='euclidean', cmap='coolwarm', 
                    figsize=(15, 8), col_cluster=False, row_cluster=False)

    # set the gridspec to only cover half of the figure
    g.gs.update(left=0.05, right=0.45)
    g.ax_cbar.set_position((0.5, 0.25, .01, .4))
    g.ax_heatmap.set_title('Community centers')
    g.ax_heatmap.set_ylabel('Community')
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.ax_heatmap.set_xlabel(f'{communities_type}')
    g.ax_heatmap.set_xticklabels(cellTypes, rotation=90)

    #create new gridspec for the right hand side
    gs2 = matplotlib.gridspec.GridSpec(1,1, left=0.6)
    # create axes within this new gridspec
    ax2 = g.fig.add_subplot(gs2[0], aspect='equal')

    # move ax2 down a bit:
    pos1 = ax2.get_position() # get the original position
    pos2 = [pos1.x0, pos1.y0 - 0.1,  pos1.width, pos1.height]
    ax2.set_position(pos2) # set a new position

    # plot the scatter plot
    sns.scatterplot(data=image_data, x=args.x_id_col, y = args.y_id_col, hue= 'community', palette='tab20', s=5, ax=ax2)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0,title='Community', frameon=False)
    plt.tight_layout
    plt.show()
    return g

def main(args):

    print('Performing community detection...')
    data = pd.read_csv(args.input)
    communities_type = args.communities_type #'majorType'
    type_cols = [col for col in data.columns if f'{communities_type}' in col]
    data.dropna(subset=type_cols, inplace=True)
    cellTypes = data[communities_type].unique().tolist()

    # count frequency of each cell tye in every neighbourhood window:
    frequencies = cell_frequencies(data, cellTypes, type_cols)

    # assign communities to each cell throughkmeans clustering of neighbourhood window frequencies:
    community_assignments, community_centers = find_communities(frequencies, n_communities=args.n_communities)
    data['community'] = community_assignments

    if not os.path.exists(args.results_outdir):
        os.mkdir(args.results_outdir)
    if not os.path.exists(args.plot_outdir):
        os.mkdir(args.plot_outdir)
    print('Done!')

    print('Saving results...')
    # save community centers:
    center_df = pd.DataFrame(community_centers, columns=cellTypes)
    center_df.to_csv(f'{args.results_outdir}/community_centers.csv', index=False)


    print('Making spatial plots of communities...')
    # make spatial plots of communities:
    for i, imagename in enumerate(data[args.image_id_col].unique()):
        image_data = data[data[args.image_id_col] == imagename]

        fig = community_plot(community_centers, image_data, cellTypes, communities_type)
        fig.savefig(f'{args.plot_outdir}/{imagename}_community_plot.pdf')
        fig.savefig(f'{args.plot_outdir}/{imagename}_community_plot.png')

        # save community assignments:
        image_data.to_csv(f'{args.results_outdir}/{imagename}_community_assignments.csv', index=False)
    print('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, required=True, help='Path to input csv file')
    parser.add_argument('--communities_type', '-c', type=str, required=True, help='Column header prefix for cell type to use for community detection e.g. "cellType"')
    parser.add_argument('--n_communities', '-n', type=int, required=False, help='Number of communities to detect', default=10)
    parser.add_argument('--plot_outdir', '-o', type=str, help='Path to output directory for community plots', default='./community_plots')
    parser.add_argument('--results_outdir', '-r', type=str, help='Path to output directory for community results', default='./community_results')
    parser.add_argument('--image_id_col', '-id', type=str, help='Column in input dataframe specifying image id', default='imagename')
    parser.add_argument('--x_id_col', '-x', type=str, help='Column in input dataframe specifying x coordinate of cell center', default='centerX')
    parser.add_argument('--y_id_col', '-y', type=str, help='Column in input dataframe specifying y coordinate of cell center', default='centerY')
    args = parser.parse_args()
    main(args)