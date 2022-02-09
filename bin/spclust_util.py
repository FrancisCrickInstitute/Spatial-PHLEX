#!/usr/bin/env python

import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skimage.io as io
from sklearn import metrics
from sklearn.cluster import DBSCAN, KMeans
from tqdm import *

import alphashape
from descartes import PolygonPatch
from shapely.geometry import MultiPolygon, Point, Polygon

########################
# FUNCTION DEFINITIONS #
########################

def do_clustering(cell_type_df, eps_param, min_samples=1):
    """
    Performs DBSCAN clustering given a cell typing dataframe and returns list of labels of clusters that each cell belongs to.
    """

    n_cells = len(cell_type_df.index)
    print('n_cells = ', n_cells)

    # reset index so iteration on cell type df works:
    cell_type_df = cell_type_df.reset_index()

    # create list of X,Y positions in [[X0,Y0], [X1,Y1]...[Xi,Yi]] format for clustering:
    cell_positions = []
    for i in range(n_cells):
        cell_positions.append([cell_type_df['centerX'][i], cell_type_df['centerY'][i]])

    # perform DBSCAN clustering:
    clustering = DBSCAN(eps=eps_param, min_samples=min_samples).fit(cell_positions)
    labels = clustering.labels_
    return labels


def plot_clusters(clustered_df, cluster_labels, image_shape, sample_name, clustering_cell_type, outdir, alphashape_param = 0.05):
    """
    Plots spatial clusters of cells with an alphashape pertaining to each spatial cluster. 
    Counts all other cell types which lie within alphashape and outputs dataframe.
    """
    if os.path.exists(outdir) != True:
        os.makedirs(outdir)

    # create figure:    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (15,15))

    # get unique cluster labels:
    unique_cluster_labels = np.unique(cluster_labels)

    # loop through labels:
    for label in unique_cluster_labels:

        # get df for each unique label:
        cluster_df = clustered_df[clustered_df['dbscan_cluster'] == label]

        # get points of cluster:
        X = cluster_df['centerX'].values
        Y = cluster_df['centerY'].values       
        points = list(zip(X, Y))

        # only proceed if cells exist:
        if len(X) > 0:

            # create alphashape:
            alpha_shape = alphashape.alphashape(points, alphashape_param)
            
            if alpha_shape.geom_type in ['Polygon', 'MultiPolygon']: # only process Polygon and Multipolygons i.e. ignore lines of cells which cannot contain other cells
              
                # plot points and add patch of alphashape:
                ax.scatter(cluster_df['centerX'].values, cluster_df['centerY'].values, alpha=0.5)
                ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))

            else:
                ax.scatter(cluster_df['centerX'].values, cluster_df['centerY'].values, alpha=0.5)

    # update plot with title etc and save:
    title = '{}_{}_alpha_{}.png'.format(sample_name, clustering_cell_type, alphashape_param)
    ax.set_title(title, fontsize=18)
    ax.set_xlabel('Centroid X µm')
    ax.set_ylabel('Centroid Y µm')
    ax.set_ylim(0,image_shape[0])
    ax.set_xlim(0,image_shape[1])
    plt.savefig('{}/{}_{}_alpha_{}.png'.format(outdir,sample_name, clustering_cell_type, alphashape_param))
    plt.close()

def assemble_cluster_polygons(clustered_df, cluster_labels, sample_name, clustering_cell_type, alpha = 0.05):
    """
    Generate alphashape polygons for all cells with cluster labels for a give image and clustering cell type.
    Returns a list of polygons.
    """

    # get unique cluster labels:
    unique_cluster_labels = np.unique(cluster_labels)

    # ignore negative [i.e. noise] points:
    unique_cluster_labels = unique_cluster_labels[unique_cluster_labels >= 0]

    # create list of all alphashapes:
    all_alphashapes = []

    # loop through labels:
    for label in unique_cluster_labels:

        # get df for each unique label:
        # print('processing label {}'.format(label))
        cluster_df = clustered_df[clustered_df['dbscan_cluster'] == label]

        # define cluster_id:
        cluster_id = '{}_{}_{}'.format(sample_name, clustering_cell_type, label)
        cluster_id = cluster_id.replace(' ', '_')
        # print('\ncluster_id: ', cluster_id)

        # get (x,y) values for cells with label
        X = cluster_df['centerX'].values
        Y = cluster_df['centerY'].values

        # get points of cluster:
        points = list(zip(X, Y))

         
        # only proceed if cells exist:
        if len(X) > 0:

            # create alphashape:
            alpha_shape = alphashape.alphashape(points, alpha)
            
            if alpha_shape.geom_type in ['Polygon', 'MultiPolygon']: # only process Polygon and Multipolygons i.e. ignore lines of cells which cannot contain other cells
                
                # add to list of all alphashapes if already a poly or multipolgyon:
                all_alphashapes.append(alpha_shape)

    return all_alphashapes

def cells_in_alphashape(alphashape, cell_list, cellType):
    """
    Counts the number of cells of each cellType within a given alphashape.
    """
    print('{} {} cells in image'.format(len(cell_list), cellType))
    count = 0
    for i in range(len(cell_list)):
        test = alphashape.contains(Point(cell_list[i]))
        if test == True:
            count += 1
        # print(test)
    return count

def distances_to_alphashape(alpha_shape, cell_list):
    distances = []
    for i in range(len(cell_list)):
        cellPoint = Point(cell_list[i])
        test = alpha_shape.contains(cellPoint)
        if test == True:
            dist = cellPoint.distance(alpha_shape)
            distances.append(dist)
            # dist = distance_to_boundary(alpha_shape, cellPoint)
            print("Distance to boundary of interior cell {} is {}.".format(cellPoint, dist)) # define as negative for interior cell?
        else:
            dist = cellPoint.distance(alpha_shape)
            print("Distance to boundary of external cell {} is {}.".format(cellPoint, dist))
            distances.append(dist)
    
    return distances

def dist_to_alphashape(alpha_shape, cell):

    cellPoint = Point(cell)
    test = alpha_shape.contains(cellPoint)
    if test == True:
        dist = -alpha_shape.boundary.distance(cellPoint) # poly.boundary.distance(pt)
        # dist = distance_to_boundary(alpha_shape, cellPoint)
        # print("Distance to boundary of interior cell {} is {}.".format(cellPoint, dist)) # define as negative for interior cell?
    else:
        dist = cellPoint.distance(alpha_shape)
        # print("Distance to boundary of external cell {} is {}.".format(cellPoint, dist))
    
    return dist

def distance_to_boundary(poly, point):
    distance = poly.exterior.distance(point)
    return distance


def cluster_stats(labels):
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    print('Estimated number of clusters: %d' % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)
    return n_clusters_


def get_image_shape(imagepath):
    im = io.imread(imagepath)
    shape = im.shape
    return shape

def get_image_shape_from_metadata(metadata, imagename):
    width = int(metadata[metadata['imagename'] == imagename]['image_width'].item())
    height = int(metadata[metadata['imagename'] == imagename]['image_height'].item())
    shape = (height, width)
    return shape

def assign_to_spclust(cell_list, spatial_cluster_polygons):
    """
    Given a list of cells and a list of polygons corresponding to the boundary of a spatial cluster, 
    test which (if any) polygon boundary the cell is located within.
    If the cell is within a spatial cluster boundary, measure the (min) distance to the boundary.
    If the cell is not within any spatial cluster, measure the distance to the nearest boundary.
    """
    # loop over all cells in image:
    alphashape_assignments = []
    cluster_areas = []
    distances_to_cluster_boundary = []
    nearest_cluster = []
    print('\nPerforming assignment of cells to spatial clusters.')
    for i in tqdm(range(len(cell_list)), ascii=True):
        cell = Point(cell_list[i])
        
        # loop over all alphashapes and test if cell within boundary:
        j = 0
        while j < len(spatial_cluster_polygons):
            spclust_polygon = spatial_cluster_polygons[j]
            spclust_id = j
            # check if alphashape either contains or touches cell: (this can also be done by using .contains() with .buffer(0.5) but is MUCH slower)
            test = spclust_polygon.contains(cell) or spclust_polygon.touches(cell)
            area = spclust_polygon.area
            if test == True:
                alphashape_assignments.append(spclust_id)
                cluster_areas.append(area)
                # measure distance:
                distance = dist_to_alphashape(spclust_polygon, cell)
                distances_to_cluster_boundary.append(distance)
                nearest_cluster.append(spclust_id)
                break
            j+=1
        else:
            # give cells that do not fall inside any clusters value -1:
            alphashape_assignments.append(-1) 
            cluster_areas.append(np.nan)
            
            # loop over all alphashapes and measure distance from cell to boundary:
            ext_dists = []
            for k in range(len(spatial_cluster_polygons)):
                dist = dist_to_alphashape(spatial_cluster_polygons[k], cell)
                ext_dists.append(dist)

            # add minimum distance to cluster (alphashape) boundary to list, and id of that cluster:
            if len(ext_dists) > 0:
                min_dist = min(ext_dists)
                nearest_clust_id = ext_dists.index(min_dist)
                distances_to_cluster_boundary.append(min_dist)
                nearest_cluster.append(nearest_clust_id)
            else:
                distances_to_cluster_boundary.append(np.nan)
                nearest_cluster.append(np.nan)

    return alphashape_assignments, cluster_areas, nearest_cluster, distances_to_cluster_boundary

def save_cluster_alphashapes(alphashapes, cellType, imagename, out_dir):
    """
    Save the alphashapes of the spatial clusters to csv in well-known text format.
    """
    cluster_ids = np.arange(len(alphashapes)).tolist()
    wkts = [alphashapes[i].wkt for i in range(len(alphashapes))]
    d = {'cluster_id': cluster_ids,
            'well-known_txt': wkts}
    df = pd.DataFrame(d)
    df['cellType'] = cellType
    df['image'] = imagename
    print(df)
    spath = os.path.join(out_dir, '{}_{}_alphashape_polygons_wkt.csv'.format(imagename, cellType))
    df.to_csv(spath)