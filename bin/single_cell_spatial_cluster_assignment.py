#!/usr/bin/env python

import os
from re import I
import sys
import argparse

import numpy as np
import pandas as pd
from sklearn.feature_extraction import img_to_graph
from tqdm import *
import skimage.io as io

from spclust_util import (assemble_cluster_polygons, assign_to_spclust,
                          do_clustering, get_image_shape, get_image_shape_from_metadata, plot_clusters,
                          save_cluster_alphashapes, pipeline_make_label)


class config(object):

    def __init__(self, args):

        # study cohort (from command line):
        self.IMAGENAME = args.imagename

        # TODO: move cohort and panel to nextflow script and pass root_out constructed from them through command line
        # self.COHORT = c

        # # IMC panel (from command line):
        # self.PANEL = p
        # print(self.COHORT, self.PANEL)

        # what level of phenotyping, majorType or cellType:
        self.PHENOTYPING_LEVEL = args.phenotyping_level

        # EPS density parameter:
        self.EPS = args.eps

        # minimum samples for clustering
        self.MIN_S = args.min_s

        # path to dataframe file containing phenotype locations:
        self.OBJECTS_FILEPATH = args.objects_filepath

        # object table separator:
        self.OBJECT_SEP = args.objects_sep

        # path to metadata file:
        self.metadata = args.metadata_filepath

        # separator of metadata file:
        self.MSEP = args.metadata_sep

        # alphashape curvature parameter:
        self.ALPHA = 0.05  

        # base output directory:
        self.ROOT_OUTDIR = f'{args.root_outdir}/{args.phenotyping_level}/dbscan_{args.eps}/min_size_{args.min_s}/alpha_{args.alpha}'

        

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
    metadata = pd.read_csv(CONFIG.metadata, sep = CONFIG.MSEP, encoding='latin1')
    
    # define imagename to be processed:
    imagename = CONFIG.IMAGENAME

    imshape = get_image_shape_from_metadata(metadata, imagename) # (1747,1756) #

    # read df:
    phenotype_df = pd.read_csv(CONFIG.OBJECTS_FILEPATH, sep=CONFIG.OBJECT_SEP, encoding='latin1')

    print(phenotype_df)

    # define probe cell types and clustering cell types
    probe_cell_types = phenotype_df[CONFIG.PHENOTYPING_LEVEL].unique() 
    clustering_cell_types = phenotype_df[CONFIG.PHENOTYPING_LEVEL].unique() # ['Epithelial cells']

    # create base output directory:
    if os.path.exists(CONFIG.ROOT_OUTDIR) != True:
        os.makedirs(CONFIG.ROOT_OUTDIR)  

    print("probe cell types:", probe_cell_types)
    print("clustering cell types: ", clustering_cell_types)

    # loop over all cellTypes to be clustered:
    for cType in clustering_cell_types:

        # PROCESS IMAGE #
               
        # create out directory:
        out_dir = os.path.join(CONFIG.ROOT_OUTDIR, '{}_clustering/{}'.format(cType, imagename))
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
            alphashape_df = save_cluster_alphashapes(spatial_cluster_polygons, cType, imagename, out_dir)

            # assign + measure distances of single cells to dbscan clusters:
            print('\nPerforming assignment of cells to spatial clusters.')
            spclust_assignments, cluster_areas, nearest_cluster, distances = assign_to_spclust(all_cells_list, spatial_cluster_polygons)

            # append info to cellTyping dataframe:
            image_df['{}_spatial_cluster_id'.format(cType)] = spclust_assignments
            image_df['{}_cluster_area'.format(cType)] = cluster_areas
            image_df['nearest {}_cluster_id'.format(cType)] = nearest_cluster
            image_df['distance to nearest {}_cluster_boundary'.format(cType)] = distances

            # Save alphashape mask:
            alphashape_spath = os.path.join(out_dir, f'{imagename}_{cType}_alphashape_polygons_label.tiff')
            pipeline_make_label(alphashape_df, imshape, alphashape_spath)

           
        
        else:
            # if no epithelial cells present, isntantiate empty dbscan cluster column and assign appropriate values to other cell types:
            clustering_cells_df['dbscan_cluster'] = []
            image_df['{}_spatial_cluster_id'.format(cType)] = -1
            image_df['{}_cluster_area'.format(cType)] = np.nan
            image_df['nearest {}_cluster_id'.format(cType)] = -1
            image_df['distance to nearest {}_cluster_boundary'.format(cType)] = np.nan

            # empty alphashape mask:
            alphashape_mask = np.zeros(imshape)
            alphashape_spath = os.path.join(out_dir, f'{imagename}_{cType}_alphashape_polygons_label.tiff')
            io.imsave(alphashape_spath, alphashape_mask)


        # save dataframes to out_dir:
        spath = os.path.join(out_dir, '{}_object_cluster_assignment.csv'.format(imagename))
        image_df.to_csv(spath, sep=CONFIG.OBJECT_SEP)
        print('\n{} done.'.format(imagename))
        

        # PLOT CLUSTERS WITH ALPHASHAPE:         
        plot_clusters(clustering_cells_df, clustering_cells_df['dbscan_cluster'].values, imshape, imagename, cType, out_dir, alphashape_param=CONFIG.ALPHA)
        
        print('\nDone.')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--imagename', type=str, default='test_image', help='name of image to be processed')
    parser.add_argument('--objects_filepath', type=str, default='test_objects.csv', help='path to objects file')
    parser.add_argument('--objects_sep', type=str, default=';', help='separator in objects file')
    parser.add_argument('--metadata_filepath', type=str, default='metadata.txt', help='path to metadata file')
    parser.add_argument('--metadata_sep', type=str, default='\t', help='separator in metadata file')
    parser.add_argument('--root_outdir', type=str, default='.', help='path to output directory')
    parser.add_argument('--eps', type=float, default=25, help='eps parameter for dbscan')
    parser.add_argument('--min_s', type=int, default=1, help='min number of cells in a cluster')
    parser.add_argument('--alpha', type=float, default=0.05, help='alpha parameter for alphashape')
    parser.add_argument('--phenotyping_level', type=str, default='cellType', help='what level of phenotyping, majorType or cellType')
    args = parser.parse_args()


    # create configuration based on input args:
    CONFIG = config(args)

    # pass CONFIG to main():
    main(CONFIG)

