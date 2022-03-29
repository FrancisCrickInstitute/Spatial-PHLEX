#!/usr/bin/env python

import os
import sys
from collections import Counter

import cugraph as cg
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData
from numpy.random import default_rng
from squidpy.gr._utils import (_assert_categorical_obs,
                               _assert_non_empty_sequence, _get_valid_values)
from tqdm import *


def assign_cell_categories(df, typing = 'new'):
    if typing == 'new':
        ### ASSIGN BROAD CELL CATEGORIES ACCORDING TO CELL CATEGORY DICT:
        cell_category_dict = dict({'Progenitor':
        ['CD4 Tcm' 'Naive CD8 T cells', 'Naive CD4 T cells'],
        'Resident effector':
        ['CD8 Trm', 'CD57+ CD8 Trm'],
        'Effector':
        ['Cytotoxic CD8 T cells', 'Cytotoxic CD4 T cells', 'GZMB+ CD8 T cells', 'GZMB+ CD4 T cells', 'CD57+ CD4 T cells', 'CD57+ CD8 T cells', 'PD1+CD27+ CD8 T cells', 'PD1+CD27+ CD4 T cells', 'CD4 Tem'],
        'Tregs':['Tregs'],
        'Dysfunctional':
        ['CD8 Exhausted TDT', 'CD4 Exhausted TDT'],
        'Other':
        ['CD4 (Other)', 'CD8 (Other)']})

        for key, val in cell_category_dict.items():
            df.loc[df.loc[:, 'cellType'].isin(val), 'cellClass'] = key
        print(cell_category_dict)
        return df

    elif typing == 'old':
    # CELL CATEGORIES [OLD TYPING -- from Katey Enfield]
        cell_category_dict = dict({
        'CD8_progenitor': ['CD8 Tcm', 'Naive CD8 T cells', 'PD1+TCF1+ CD8 T cells'],
        'CD4_progenitor': ['CD4 Tcm', 'PD1+TCF1+ CD4 T cells'],
        'CD8_resident_effector': ['CD8 Trm', 'PD1+CD27+ CD8 Trm', 'CD39+PD1+ CD8 Trm'],
        'CD4_resident_effector': ['CD4 Trm'],
        'CD8_effector': ['Cytotoxic CD8 T cells', 'CD39+CD57+ CD8 T cells'],
        'CD4_effector': ['Cytotoxic CD4 T cells', 'CD39+CD57+ CD4 T cells'],
        'Tregs': ['Tregs'],
        'CD8_dysfunctional': ['CD8 Tdys', 'CD8 Exhausted TDT'],
        'CD4_dysfunctional': ['CD4 Tdys', 'CD4 Exhausted TDT']})
        for key, val in cell_category_dict.items():
            df.loc[df.loc[:, 'cellType'].isin(val), 'cellClass'] = key
        print(cell_category_dict)
        return df
    else:
        raise ValueError('ValueError: typing not recognised.')
 

def make_legible(string):
    string = string.replace('&', '&\n')
    return string


def df_select(df, cat, val):
    return df[df[cat] == val]


def compute_cohort_centralities(objects, markers, phenotyping_level='cellType', radius=35):

    imagenames = objects['imagename'].unique()

    cohort_centralities = []
    for i, imagename in enumerate(imagenames):
            
        # annotate adata with phenotyping info
        image_objects = objects[objects['imagename'] == imagename].astype('category')
        image_cTypes = image_objects[phenotyping_level].dropna().unique()
    

        if len(image_cTypes) > 0: # only proceed if there are cells of the specified phenotype level in the image

            features = image_objects[markers].values
            majorTypes = image_objects['majorType'].values
            cellTypes = image_objects['cellType'].values
            cellClasses = image_objects['cellClass'].values
            positivity = image_objects['positive'].values
            coords = np.asarray(list(zip(image_objects['centerX'].values, image_objects['centerY'].values)))

            # create anndata object:
            adata = AnnData(features, obsm={"spatial": coords}, obs={'majorType':majorTypes, 'cellType':cellTypes, 'positive':positivity, 'cellClass':cellClasses})
            print(adata)
            
            ## derive spatial neighbours graph:
            sq.gr.spatial_neighbors(adata, radius=radius)
            
            ## compute centrality scores for image
            sq.gr.centrality_scores(
                adata,
                cluster_key=phenotyping_level,
            )
            
            ## obtain centrality scores
            image_centralities = adata.uns[f'{phenotyping_level}_centrality_scores']
            image_centralities['imagename'] = imagename
            
            cohort_centralities.append(image_centralities)
    
    cohort_centralities = pd.concat(cohort_centralities)

    return cohort_centralities


def compute_spatial_graph(objects, imagename, markers, phenotyping_level='cellType', radius=35):

     # annotate adata with phenotyping info
    image_objects = objects[objects['imagename'] == imagename].astype('category')
    image_cTypes = image_objects[phenotyping_level].dropna().unique()
    

    if len(image_cTypes) > 0: # only proceed if there are cells of the specified phenotype level in the image

        features = image_objects[markers].values
        majorTypes = image_objects['majorType'].values
        cellTypes = image_objects['cellType'].values
        cellClasses = image_objects['cellClass'].values
        positivity = image_objects['positive'].values
        coords = np.asarray(list(zip(image_objects['centerX'].values, image_objects['centerY'].values)))

        # create anndata object:
        adata = AnnData(features, obsm={"spatial": coords}, obs={'majorType':majorTypes, 'cellType':cellTypes, 'positive':positivity, 'cellClass':cellClasses})
        print(adata)
        
        ## derive spatial neighbours graph:
        sq.gr.spatial_neighbors(adata, radius=radius, coord_type = 'generic')

    return adata


def compute_nn_graph(objects, imagename, markers, phenotyping_level='cellType', n_neighs=6):

     # annotate adata with phenotyping info
    image_objects = objects[objects['imagename'] == imagename].astype('category')
    image_cTypes = image_objects[phenotyping_level].dropna().unique()
    

    if len(image_cTypes) > 0: # only proceed if there are cells of the specified phenotype level in the image

        features = image_objects[markers].values
        majorTypes = image_objects['majorType'].values
        cellTypes = image_objects['cellType'].values
        cellClasses = image_objects['cellClass'].values
        positivity = image_objects['positive'].values
        coords = np.asarray(list(zip(image_objects['centerX'].values, image_objects['centerY'].values)))

        # create anndata object:
        adata = AnnData(features, obsm={"spatial": coords}, obs={'majorType':majorTypes, 'cellType':cellTypes, 'positive':positivity, 'cellClass':cellClasses})
        print(adata)
        
        ## derive spatial neighbours graph:
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighs)

    return adata


def compute_shortest_paths(G, node_ids, source = 0):
    sps = cg.traversal.sssp(G, source=source)
    sps = merge_node_ids_to_shortest_paths(sps, node_ids)
    return sps
    

def merge_node_ids_to_shortest_paths(shortest_paths, node_ids):
    # node_ids = pandas dataframe with nodeid in integer format as column called 'vertex'
    node_ids['vertex'] = np.arange(0,len(node_ids.index))
    shortest = pd.merge(shortest_paths, node_ids, how='left', on='vertex')
    return shortest


def min_path_to_epithelial(shortest_paths):
    min_path = shortest_paths[shortest_paths['cellType'] == 'Epithelial cells']['distance'].min()
    return int(min_path)


def get_predecessor(df, minpath_to_target):
    print(df)
    predecessor = df[(df['distance'] == minpath_to_target) & (df['cellType'] == 'Epithelial cells')]['predecessor'].values.tolist()
    # print('initial_predecessor: ', predecessor)
    return predecessor


def get_chain_predecessor(df, predecessor):
    chain_predecessor = df[(df['vertex'] == predecessor)].iloc[0]['predecessor'].item()
    return chain_predecessor


def followchain(shortest_paths, minpath_to_target, source_cell, source_vertex, target_cell = 'Epithelial cells', phenotype_level='cellType'):
    '''
    The min_path_to_target is e.g. 3 jumps from the source to the target cell; 
    there may be more than one path of the same length.
    
    Todo: return vertex ids for subsequent referencing of e.g. positive marker expression'''
    sorted_shortest_paths = shortest_paths.sort_values(by='distance')
    sorted_shortest_paths = sorted_shortest_paths[(sorted_shortest_paths['distance'] == minpath_to_target) & (sorted_shortest_paths['cellType'] == 'Epithelial cells')]
    print('SORTED SHORTEST PATHS\n', sorted_shortest_paths)
    target_vertex = sorted_shortest_paths['vertex'].iloc[0]
        
    # print(minpath_to_target)
    if minpath_to_target == 1:
        vertexes = [target_vertex, source_vertex]
        chaincells = [target_cell, source_cell]
        # print(chaincells)
        return chaincells, vertexes
    
    else:
        ## calculate internal chain:
        
        if sorted_shortest_paths['predecessor'].max() >= 0: # possible to have a df with target cells that are disconnected
            chaincells = [target_cell] # store cells in chain between starting cell and destination 
            vertexes = [target_vertex]

            predecessors = get_predecessor(sorted_shortest_paths, minpath_to_target) # get list of initial predecessor cells            
            p = predecessors[0]

            for j in range(minpath_to_target-1):
            # take first predecessor if more than one:
            
                if p != -1: ## proceed only if the predecessor is not the source cell
                    predecessor_nodeType = get_predecessor_nodeType(shortest_paths, p, phenotype_level=phenotype_level)
                    chaincells.append(predecessor_nodeType)
                    vertexes.append(p)
                    p = get_chain_predecessor(shortest_paths, p)
                
            chaincells.append(source_cell)
            vertexes.append(source_vertex)
            return chaincells, vertexes

                
def relabel_vertex_chain(vertex_chain, relabel_dict):
    relabelled_chain = [relabel_dict[item] for item in vertex_chain]
    return relabelled_chain


def recursive_predecessor(df, predecessor):
    """This is a recursive function
    to find the factorial of an integer"""
    predecessor = df[(df['vertex'] == predecessor)].iloc[0]['predecessor'].item()
    # print(predecessor)
    if predecessor == -1:
        return 1
    else:
        return recursive_predecessor(df, df[df['predecessor'] == predecessor])


def get_predecessor_nodeType(df, predecessor, phenotype_level='cellType'):
    nodetype = df[df['vertex'] == predecessor].iloc[0][phenotype_level]
    return nodetype


def get_graph_node_ids(spg):
    """
    Produce dataframe of graph node ids from squidpy spatial graph.
    """
    node_ids = pd.DataFrame(spg.obs[['majorType', 'cellType', 'positive', 'cellClass']])
    node_ids['vertex'] = np.arange(0,len(node_ids.index))
    return node_ids


def get_node_ids_2(objects, imagename, attributes):
    """
    Produce dataframe of graph node ids from objects table and specified attributes.
    """
    image_objects = objects[objects['imagename'] == imagename].reset_index()
    image_objects = image_objects[attributes]
    image_objects['vertex'] = np.arange(0,len(image_objects.index))
    return image_objects


def get_node_ids_3(objects, imagename, attributes, nodelabels):
    """
    Produce dataframe of graph node ids from objects table and specified attributes, and the original list of node labels from neigbouRhood out.
    """
    print("OBJECTS:\n",objects)
    print("imagename:\n",imagename)
    print("ATTRIBUTES\n",attributes)
    print("NODELABELS:\n", nodelabels)
    image_objects = objects[objects['imagename'] == imagename].set_index('object').loc[nodelabels, attributes]
    print(image_objects)
    image_objects['vertex'] = np.arange(0,len(image_objects.index)) #image_objects.index #
    node_ids = image_objects.astype(object)
    return node_ids

    
def get_cellType_node_ids(node_ids, phenotyping_level, cellType):
    cellType_node_ids = node_ids[node_ids[phenotyping_level] == cellType]
    return cellType_node_ids


def count_barrier_cells(chain, BARRIER_TYPES):
    if len(set(BARRIER_TYPES).intersection(set(chain))) > 0:
        counts = Counter(chain)
        x = sum([counts[cell] for cell in BARRIER_TYPES])
    else:
        x = 0
    return x


def weighted_barrier_count(chain, BARRIER_TYPES):
    '''Barrier weighting is inversely proportional to its adjacency to the target cell.
    '''

    count = 0
    for i, elem in enumerate(chain[1:-1]):
        weight = 1/(i+1)
        if elem in BARRIER_TYPES:
            count += weight
    return count


def binary_barrier_call(chain, BARRIER_TYPES):
    '''Returns 1 if barrier type is in cell chain'''
    if len(set(BARRIER_TYPES).intersection(set(chain))) > 0:
        barrier = 1
    else:
        barrier = 0
    return barrier


def adjacent_barrier_call(chain, BARRIER_TYPES):
    '''Returns 1 if barrier type is adjacent to target cell (first cell in chain)
    '''
    inner_chain = chain[1:-1]
    
    if len(inner_chain) > 0:
        if inner_chain[0] in BARRIER_TYPES:
            barrier = 1
        else:
            barrier = 0
    else:
        barrier = 0
    return barrier


def degenerate_path_content(closest_epi, shortest_paths, minpath_to_epi, barrier_cells=['Myofibroblasts']):
    
    """
    closest epi = shortest path dataframe with each row an entry corresponding to an epithelial cell at minpath_to_epi
    distance from the source cell for the shortest paths analysis.
    
    shortest_paths = full shortest paths table from the sssp calculation, with assigned cell types for each vertex
    
    minpath_to_epi = minimum distance to an epithelial cell in terms of degree (i.e. hops from source cell)
    
    barrier_cells = list, list of names corresponding to cell types to be considered barriers 
    
    returns:
    myofibroblast_fraction_degen, the fraction of cells along the shortest paths to the N epithelial cells at
    minpath_to_epi distance from the source cell
    """
    
    if minpath_to_epi > 1: #
        
        path_cells = []
        for i in range(len(closest_epi)):
            p = closest_epi.iloc[i,:]['predecessor']
            pdf = shortest_paths[shortest_paths['vertex'] == p]
            chain = []

            for j in range(minpath_to_epi-1):

                if pdf['predecessor'].item() > -1:
                    pdf = shortest_paths[shortest_paths['vertex'] == p]
                    chain.append(pdf)
                    p = pdf['predecessor'].item()
                    path_cells.append(pdf)

            degenerate_path_content = pd.concat(path_cells)

            if set(barrier_cells).intersection(set(degenerate_path_content['cellType'])) != set():
                myofibroblast_fraction_degen = degenerate_path_content['cellType'].value_counts(normalize=True)[barrier_cells].sum()
            else:
                myofibroblast_fraction_degen = 0
    else:
        myofibroblast_fraction_degen = 0
    return myofibroblast_fraction_degen


def degenerate_adjacent_barrier(closest_epi, shortest_paths, minpath_to_epi, barrier_cells=['Myofibroblasts']):
    
    """
    closest epi = shortest path dataframe with each row an entry corresponding to an epithelial cell at minpath_to_epi
    distance from the source cell for the shortest paths analysis.
    
    shortest_paths = full shortest paths table from the sssp calculation, with assigned cell types for each vertex
    
    minpath_to_epi = minimum distance to an epithelial cell in terms of degree (i.e. hops from source cell)
    
    barrier_cells = list, list of names corresponding to cell types to be considered barriers 
    
    returns:
    myofibroblast_fraction_degen, the fraction of cells along the shortest paths to the N epithelial cells at
    minpath_to_epi distance from the source cell
    """
    
    if minpath_to_epi > 1: #
        
        path_cells = []
        for i in range(len(closest_epi)):
            p = closest_epi.iloc[i,:]['predecessor']
            pdf = shortest_paths[shortest_paths['vertex'] == p]
            chain = []

            for j in range(minpath_to_epi-1):

                if pdf['predecessor'].item() > -1:
                    pdf = shortest_paths[shortest_paths['vertex'] == p]
                    chain.append(pdf)
                    p = pdf['predecessor'].item()
                    path_cells.append(pdf)

            degenerate_path_content = pd.concat(path_cells)
            adjacent_cells_df = degenerate_path_content[degenerate_path_content['distance'] == (minpath_to_epi-1)]

            if set(barrier_cells).intersection(set(adjacent_cells_df['cellType'])) != set():                
                mean_adjacent_barrier = adjacent_cells_df['cellType'].value_counts(normalize=True)[barrier_cells].sum() #/ len(adjacent_cells_df)
            else:
                mean_adjacent_barrier = 0

    else:
        mean_adjacent_barrier = 0
    return mean_adjacent_barrier

def permute_phenotypes(objects, in_place=True, region = 'all'):
    '''
    Randomly permute phenotype and positivity information for all images in the objects table.
    
    Args:
    objects = cell objects dataframe with each image as a distinct imagename
    in_place = bool, if True overwrite current cell phenotype data in place, 
            otherwise add '_permuted' columns
    region = str, subset permutation to a give tissue region defined by 'region' column.
            'all' = do not subset, permute all cells.
        
    '''
    if region == 'all':
        for imagename in objects['imagename'].unique():
            rng = np.random.default_rng(seed=123)
            if in_place:
                objects.loc[objects['imagename'] == imagename, 
                            ('cellType', 'majorType', 'positive')] = rng.permutation(objects.loc[objects['imagename'] == imagename,
                                                                                                 ('cellType', 'majorType', 'positive')].values)
            else:
                objects.loc[objects['imagename'] == imagename, 
                            ('cellType_permuted', 'majorType_permuted', 'positive_permuted')] = rng.permutation(objects.loc[objects['imagename'] == imagename,
                                                                                                 ('cellType', 'majorType', 'positive')].values)
    else:
        for imagename in objects['imagename'].unique():
            rng = np.random.default_rng(seed=123)
            if in_place:
                objects.loc[(objects['imagename'] == imagename) & (objects['region'] == region), 
                            ('cellType', 'majorType', 'positive')] = rng.permutation(objects.loc[(objects['imagename'] == imagename) & (objects['region'] == region),
                                                                                                 ('cellType', 'majorType', 'positive')].values)
            else:
                objects.loc[(objects['imagename'] == imagename) & (objects['region'] == region), 
                            (f'cellType_{region}_permuted', f'majorType_{region}_permuted', f'positive_{region}_permuted')] = rng.permutation(objects.loc[(objects['imagename'] == imagename) & (objects['region'] == region),
                                                                                                 ('cellType', 'majorType', 'positive')].values)
                objects.loc[(objects['imagename'] == imagename) & ~(objects['region'] == region), 
                            (f'cellType_{region}_permuted', f'majorType_{region}_permuted', f'positive_{region}_permuted')] = objects.loc[(objects['imagename'] == imagename) & ~(objects['region'] == region),
                                                                                                 ('cellType', 'majorType', 'positive')].values

    return objects