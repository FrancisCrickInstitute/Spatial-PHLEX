#!/usr/bin/env python

import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
import argparse
import os

def calc_neighbours(measurements, n_neighbours=5, object_identifier = 'object', x_id_col='centroid-x', y_id_col='centroid-y', phenotyping_ids=None):
    # calculate object neighbours based on cell centroid position
    """Calculate the nearest neighbours of each object in a dataframe based on the centroid position of the objects. 
    Keep track of the distance and phenotyping information of the neighbours.

    Args:
        measurements (pd.DataFrame): A dataframe of cell spatial measurements from a single image.
        n_neighbours (int, optional): Number of nearest neighbours of the central cell to calculate. Defaults to 5.
        object_identifier (str, optional): Column header corresponding to a unique numerical label for a given cell. Defaults to 'object'.
        x_id_col (str, optional): Column header for x coordinate values. Defaults to 'centroid-x'.
        y_id_col (str, optional): Column header for y coordinate values. Defaults to 'centroid-y'.
        phenotyping_ids (list, optional): A list of column headers that contain qualitative phenotypic information for 
        the cells. e.g. ['majorType', 'cellType']. Creates columns of the phenotypes of the nearest neighbours. Defaults to None.

    Returns:
        pd.DataFrame: The nearest neighbour dataframe with columns for the object identifier, the object identifiers of the nearest neighbours, 
        plus the distance to the nearest neighbours, and phenotyping information of the nearest neighbours if requested.
    """    

    measurements = measurements.reset_index()
    coords = measurements[[x_id_col, y_id_col]].values
    nbrs = NearestNeighbors(n_neighbors=n_neighbours, algorithm='kd_tree').fit(coords)
    distances, indices = nbrs.kneighbors(coords)

    # 
    neighbour_labels = np.zeros_like(indices)
    for i in range(n_neighbours):
        neighbour_labels[:,i] = measurements[object_identifier][indices[:,i]]

    neighbour_columns=[object_identifier] + [f'{object_identifier}_neighbour_{i}' for i in range(1, n_neighbours)]

    neighbour_df = pd.DataFrame(data=neighbour_labels, columns=neighbour_columns)
    distance_columns = [f'distance_neighbour_{i}' for i in range(1, n_neighbours)]
    distance_df = pd.DataFrame(data=distances[:,1:], columns=distance_columns)

    # merge result:
    nearest_neighbours = pd.merge(neighbour_df, distance_df,  left_index=True, right_index=True)

    # add phenotyping columns to neighbour_df
    # create char array for each neighbour
    if phenotyping_ids is not None:  
    
        for phenotyping_column in phenotyping_ids:

            neighbour_columns=[phenotyping_column] + [f'{phenotyping_column}_neighbour_{i}' for i in range(1, n_neighbours)]

            neighbour_ids = np.zeros_like(indices).astype(str)
            for i in range(n_neighbours):
                neighbour_ids[:,i] = measurements[phenotyping_column][indices[:,i]]

            phenotyping_neighbours = pd.DataFrame(data = neighbour_ids, columns=neighbour_columns)
            nearest_neighbours = pd.merge(nearest_neighbours, phenotyping_neighbours,  left_index=True, right_index=True)

    # merge coords back in:
    nearest_neighbours[[x_id_col, y_id_col]] = measurements[[x_id_col, y_id_col]]
    
    return nearest_neighbours


def main(args):

    # load data
    data = pd.read_csv(args.measurements, sep=args.measurements_sep)

    print('Calculating nearest neighbours...')
    print(data.head())

    # calculate neighbours
    all_neighbours = []

    for imagename in data[args.image_id_col].unique():
        print(f'Processing {imagename}...')
        image_data = data[data[args.image_id_col]==imagename]
        nearest_neighbours = calc_neighbours(image_data, 
                                            n_neighbours=args.n_neighbours, 
                                            object_identifier=args.object_identifier, 
                                            x_id_col=args.x_id_col, 
                                            y_id_col=args.y_id_col, 
                                            phenotyping_ids=args.phenotyping_ids)
        
        nearest_neighbours[args.image_id_col] = imagename
        all_neighbours.append(nearest_neighbours)

    nn_data = pd.concat(all_neighbours)
    cols = list(nn_data.columns)
    cols = [cols[-1]] + cols[:-1] # move imagename column to front for convenience
    nn_data = nn_data[cols]

    # save result
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    input_fname = os.path.splitext(os.path.basename(args.measurements))[0]
    save_path = os.path.join(args.outdir, f'{input_fname}_{args.n_neighbours}_nearest_neighbours.csv')
    nn_data.to_csv(save_path, index=False)

    print('Done!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Calculate nearest neighbours of cells based on cell centroid position.')
    parser.add_argument('--measurements', type=str, help='Path to cell measurements file.')
    parser.add_argument('--measurements_sep', type=str, default=',', help='Delimiter for cell measurements file.')
    parser.add_argument('--n_neighbours', type=int, default=10, help='Number of nearest neighbours to calculate.')
    parser.add_argument('--object_identifier', type=str, default='object', help='Column header corresponding to a unique numerical label for a given cell.')
    parser.add_argument('--image_id_col', type=str, default='imagename', help='Column header for image name.')
    parser.add_argument('--x_id_col', type=str, default='centroid-x', help='Column header for x coordinate values.')
    parser.add_argument('--y_id_col', type=str, default='centroid-y', help='Column header for y coordinate values.')
    parser.add_argument('--phenotyping_ids', nargs='+', default=None, help='A list of column headers that contain qualitative phenotypic information for the cells. e.g. --phenotyping_ids majorType cellType. Creates columns of the phenotypes of the nearest neighbours.')
    parser.add_argument('--outdir', type=str, help='Path to output directory.', default='./nearest_neighbours')
    args = parser.parse_args()
    main(args)



'''Specifically, each window was converted to a vector of length 29 containing the frequency of each of the 29 cell types among the 10 neighbors, and the windows were subsequently clustered using Python’s scikit-learn implementation of MiniBatchKMeans with k = 10. Each cell was then allocated to the CN that its surrounding window was. To validate the CN assignment, these allocations were overlaid on the original tissue H&E-stained and fluorescent images. During this process, the CN cluster that contained the imaging artifacts (cellular cluster 29) was removed.'''