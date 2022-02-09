#!/usr/bin/env python

import os
from re import I
import sys

import numpy as np
import pandas as pd
from sklearn.feature_extraction import img_to_graph
from tqdm import *

from spclust_util import (assemble_cluster_polygons, assign_to_spclust,
                          do_clustering, get_image_shape, get_image_shape_from_metadata, plot_clusters,
                          save_cluster_alphashapes)


class config(object):

    def __init__(self, X, c, p, t, o):

        # study cohort (from command line):
        self.IMAGENAME = X

        self.COHORT = c

        # IMC panel (from command line):
        self.PANEL = p
        print(self.COHORT, self.PANEL)

        # what level of phenotyping, majorType or cellType:
        self.PHENOTYPING_LEVEL = t #'cellType' 

        # EPS density parameter:
        self.EPS = 25

        # minimum samples for clustering
        self.MIN_S = 0

        # path to dataframe file containing phenotype locations:
        self.OBJECTS_FILEPATH = o #'/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/{}/imc/outputs/cell_typing/{}_cell_objects_{}_final_{}.txt'.format(self.COHORT, self.COHORT, self.COHORT, self.PANEL)

        # object table separator:
        self.OBJECT_SEP = '\t'

        # path to metadata file:
        self.METADATA = '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'

        # separator of metadata file:
        self.MSEP = '\t'

        # alphashape curvature parameter:
        self.ALPHA = 0.05  

        self.RELEASE_VERSION = '2022-02-02_release'

        # base output directory:
        self.ROOT_OUT_DIR = './single_cell_assignment/{}/{}/{}/dbscan_{}/min_size_{}/alpha_{}'.format(self.COHORT, self.PANEL, self.PHENOTYPING_LEVEL, self.EPS, self.MIN_S, self.ALPHA)

        

    def display(self):
        """Display Configuration values."""
        print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("_____CONFIG:_____")
        for a in dir(self):
            if not a.startswith("__") and not callable(getattr(self, a)):
                print("{:30} {}".format(a, getattr(self, a)))
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("\n")


########
# MAIN #
########

def main(CONFIG):

    # display configuration:
    CONFIG.display()

    # read metadata
    metadata = pd.read_csv(CONFIG.METADATA, sep = CONFIG.MSEP)
    
    # define imagename to be processed:
    imagename = CONFIG.IMAGENAME

    # read df:
    phenotype_df = pd.read_csv(CONFIG.OBJECTS_FILEPATH, sep=CONFIG.OBJECT_SEP)

    print(phenotype_df)

    # define probe cell types and clustering cell types
    probe_cell_types = phenotype_df[CONFIG.PHENOTYPING_LEVEL].unique() 
    clustering_cell_types = phenotype_df[CONFIG.PHENOTYPING_LEVEL].unique() # ['Epithelial cells'] 

    # create base output directory:
    if os.path.exists(CONFIG.ROOT_OUT_DIR) != True:
        os.makedirs(CONFIG.ROOT_OUT_DIR)  

    print("probe cell types:", probe_cell_types)
    print("clustering cell types: ", clustering_cell_types)

    # loop over all cellTypes to be clustered:
    for cType in clustering_cell_types:

        # PROCESS IMAGE #
               
        # create out directory:
        out_dir = os.path.join(CONFIG.ROOT_OUT_DIR, '{}_clustering/{}'.format(cType, imagename))
        if os.path.exists(out_dir) != True:
            os.makedirs(out_dir)
        
        # extract df for imagename:
        image_df = phenotype_df.loc[phenotype_df['imagename'] == imagename]

        print('\nimage_df:\n', image_df)

        # get dataframe of locations of all cells to be compared to clustered cell type:
        all_cells_list = list(zip(image_df['centerX'].values, image_df['centerY'].values)) # list of cells of Type

        # perform DBSCAN clustering for cellType:
        print('\nPerforming spatial clustering of {}.'.format(cType))
        clustering_cells_df = image_df.loc[image_df[CONFIG.PHENOTYPING_LEVEL] == cType]

        if len(clustering_cells_df.index) > 0:

            #obtain cluster labels with dbscan clusterinbg:
            cluster_labels = do_clustering(clustering_cells_df, CONFIG.EPS, CONFIG.MIN_S)

            # add clusters to dataframe:
            clustering_cells_df['dbscan_cluster'] = cluster_labels

            # get list of shapely polygons describing boundaries of DBSCAN clusters:
            spatial_cluster_polygons = assemble_cluster_polygons(clustering_cells_df, cluster_labels, imagename, cType, alpha = CONFIG.ALPHA)
            print('\n {} spatial clusters with > {} cells.'.format(len(spatial_cluster_polygons), CONFIG.MIN_S))

            # save alphashapes to csv dataframe:
            save_cluster_alphashapes(spatial_cluster_polygons, cType, imagename, out_dir)

            # assign + measure distances of single cells to dbscan clusters:
            print('\nPerforming assignment of cells to spatial clusters.')
            spclust_assignments, cluster_areas, nearest_cluster, distances = assign_to_spclust(all_cells_list, spatial_cluster_polygons)

            # append info to cellTyping dataframe:
            image_df['{}_spatial_cluster_id'.format(cType)] = spclust_assignments
            image_df['{}_cluster_area'.format(cType)] = cluster_areas
            image_df['nearest {}_cluster_id'.format(cType)] = nearest_cluster
            image_df['distance to nearest {}_cluster_boundary'.format(cType)] = distances
           
        
        else:
            # if no epithelial cells present, isntantiate empty dbscan cluster column and assign appropriate values to other cell types:
            clustering_cells_df['dbscan_cluster'] = []
            image_df['{}_spatial_cluster'.format(cType)] = -1
            image_df['{}_cluster_area'.format(cType)] = np.nan
            image_df['nearest {}_cluster_id'.format(cType)] = -1
            image_df['distance to nearest {}_cluster_boundary'.format(cType)] = np.nan


        # save dataframes to out_dir:
        spath = os.path.join(out_dir, '{}_object_cluster_assignment.csv'.format(imagename))
        image_df.to_csv(spath)
        print('\n{} done.'.format(imagename))
        

        # PLOT CLUSTERS WITH ALPHASHAPE:
        imshape = get_image_shape_from_metadata(metadata, imagename)  
        plot_clusters(clustering_cells_df, clustering_cells_df['dbscan_cluster'].values, imshape, imagename, cType, out_dir, alphashape_param=CONFIG.ALPHA)

        print('\nDone.')

if __name__ == "__main__":
    # create configuration based on input args:
    CONFIG = config(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

    # pass CONFIG to main():
    main(CONFIG)

