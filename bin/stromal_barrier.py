#!/usr/bin/env python

import pandas as pd
# import cugraph as cg
import networkx as nx
import numpy as np
from tqdm import *
from collections import Counter
import os, sys
import sb_functions as sb 
import argparse
import warnings

def main(args):

    GRAPH_TYPE = args.graph_type # 'nearest_neighbour' or 'spatial_neighbours', 'neighbouRhood'
    NEIGHBOURS = args.neighbours# 10
    RADIUS = args.radius # 5
    NEIGHBOURHOOD_RADIUS = args.neighbourhood_radius #5
    ADJACENCY_DATA_PATH = args.adjacency_data_path # this should be conditional argument; only used if 'neighbouRhood' the specified graph type; graph type should be specified as args input 
    ROOT_OUT = args.root_out
    OBJECTS_PATH = args.objects_path
    OBJECT_SEP = args.objects_sep #','
    PANEL = args.panel # args[3] #'p2'
    imagename = args.imagename

    CALC_CHAIN = args.calc_chain # True
    PERMUTE_PHENOTYPES = args.permute_phenotypes #False
    PERMUTATION_REGION = args.permutation_region #'all'
    EPI_NO_STROMA = False
    LARGE_CLUSTERS_ONLY = False
    SIZE_THRESH_NO_UNCLUSTERED = True
    DOMAIN_SIZE_CUTOFF = args.domain_size_cutoff
    target_cell_type = args.target_cell_type
    calculate_positivity = False
    
    if type(args.barrier_types) == list:
        BARRIER_TYPES = args.barrier_types #['Myofibroblasts']
    else:
        BARRIER_TYPES = [args.barrier_types]

    print(BARRIER_TYPES)

    phenotyping_level = args.phenotyping_level #phenotyping_level # phenotyping_level, 'majorType', 'Positive', 'cellClass'
    
    # Define immune cell subtypes to measure the 'barrier' for:
    cellTypes = [args.source_cell_type] #['Myofibroblasts'] #['Myofibroblasts', 'Stromal', 'Endothelial']
    '''['CD4 T cells', 'CD4 T cells', 'CD57+ CD4 T cells', 'Naive CD4 T cells',
       'Leukocytes - Other', 'CD4+ Myeloid cells', 'Cytotoxic CD8 T cells',
       'CD4 Tcm', 'Tregs', 'CD8 Trm', 'CD8 T cells', 'CD57+ CD8 Trm',
       'CD8 Exhausted TDT', 'Naive CD8 T cells', 'T cells DP',
       'Cytotoxic CD4 T cells', 'T cells - Other']'''


    if GRAPH_TYPE == 'spatial_neighbours':
        RESULTS_DIR = os.path.join(ROOT_OUT, GRAPH_TYPE, f'radius_{RADIUS}', '_'.join(BARRIER_TYPES))
    elif GRAPH_TYPE == 'nearest_neighbour':
        RESULTS_DIR = os.path.join(ROOT_OUT, GRAPH_TYPE, f'n_neigh_{NEIGHBOURS}', '_'.join(BARRIER_TYPES))
    elif GRAPH_TYPE == 'neighbouRhood':
        RESULTS_DIR = os.path.join(ROOT_OUT, GRAPH_TYPE, f'neighbouRhood_radius_{NEIGHBOURHOOD_RADIUS}', '_'.join(BARRIER_TYPES))
    else:
        raise ValueError("GRAPH_TYPE not recognised.")

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR, exist_ok=True)

    # Read in necessary files:
    objects = pd.read_csv(OBJECTS_PATH, sep=OBJECT_SEP, encoding='latin1', quoting=0)
    print('\nThe unique imagenames in the objects table are:\n')
    print(objects["imagename"].unique())
    
    if LARGE_CLUSTERS_ONLY:
        # assign distal stroma epithelial cells to unassigned to test
        objects.loc[(objects['imagename'] == imagename) & (objects[f'{target_cell_type}_cluster_area'] < 2000) & (objects[phenotyping_level] == f'{target_cell_type}'), phenotyping_level] = f'Unclustered {target_cell_type}'
    if EPI_NO_STROMA:
        objects.loc[(objects['domain'] == 'Distal Stroma') & (objects[phenotyping_level] == f'{target_cell_type}'), phenotyping_level] = 'Unassigned'
    
    if SIZE_THRESH_NO_UNCLUSTERED:
        # Alter unclustered
        objects.loc[(objects[f'{target_cell_type}_spatial_cluster_id'] == -1) & (objects[phenotyping_level] == f'{target_cell_type}'), phenotyping_level] = f'Unclustered {target_cell_type}'
        # alter those in small clusters:
        objects.loc[(objects[f'{target_cell_type}_cluster_area'] < DOMAIN_SIZE_CUTOFF) & (objects[phenotyping_level] == f'{target_cell_type}'), phenotyping_level] = f'Unclustered {target_cell_type}'

    ## after filtering only proceed if there are epithelial cells that pass the criteria, else raise warning:   
    if len(objects[objects[phenotyping_level] == f'{target_cell_type}'].index) > 0:

        if imagename in objects['imagename'].unique():

            if PERMUTE_PHENOTYPES == True:
                objects = sb.permute_phenotypes(objects, in_place=True, region=PERMUTATION_REGION)

            # SPECIFY IMC MARKERS:
            if PANEL == 'p1':
                markers = ['CD39', 'CD4', 'LAG3', 'casp3', 'CD31', 'CD27', 'CXCR6', 'PDL1', 'CXCR4', 'TCF1', 'CD103', 'GATA3', 'FAP1', 'ICOS', 'Vimentin', 'CD25', 'PD1', 'CD3', 'GITR', 'alphaSMA', 'CD8a', 'CTLA4', 'B2M', 'TIM3', 'Ki67', 'CD57', 'CD45RA', 'CXCL12', 'panCK', 'CCR7', 'Collagen', 'pSTAT1', 'GZMB', 'CD45', 'FOXP3']
            if PANEL == 'p2':
                markers = ['CD163', 'CLEC9a', 'KIR2DL3', 'VISTA', 'CD56', 'TIM3', 'CD45', 'CD16', 'MHCII', 'panactin', 'TCRd', 'PD1', 'CD79a', 'CD66b', 'CD14', 'MPO', 'CAIX', 'CD38', 'MCT4', 'CD206', 'IDO', 'CD68', 'PDL1', 'CD20', 'CD11c', 'CD4', 'CD8a', 'CD31', 'LAG3', 'CD103', 'CD11b', 'CD3', 'GZMB', 'panCK', 'CD73']


            ## CREATE SPATIAL GRAPH
            if GRAPH_TYPE == 'spatial_neighbours':
                spg = sb.compute_spatial_graph(objects, imagename, markers, phenotyping_level=phenotyping_level, radius=RADIUS)
                G = nx.convert_matrix.from_scipy_sparse_matrix(spg.obsp['spatial_connectivities'])
                node_ids = sb.get_graph_node_ids(spg, phenotyping_level)
                node_label_dict = dict(zip(node_ids['vertex'], node_ids.index))
            elif GRAPH_TYPE == 'nearest_neighbour':
                spg = sb.compute_nn_graph(objects, imagename, markers, n_neighs=NEIGHBOURS, phenotyping_level=phenotyping_level)
                G = nx.convert_matrix.from_scipy_sparse_matrix(spg.obsp['spatial_connectivities'])
                node_ids = sb.get_graph_node_ids(spg, phenotyping_level)
                node_label_dict = dict(zip(node_ids['vertex'], node_ids.index))
            elif GRAPH_TYPE == 'neighbouRhood':

                # specify attributes to attach to nodes:
                node_attr = [phenotyping_level]

                #read in adjacency data to n graph:
                primary_nodes = [] 
                with open(ADJACENCY_DATA_PATH,'r') as f:
                    for line in f:
                        primary_nodes.append(line.split(' ')[0])
                G = nx.read_adjlist(ADJACENCY_DATA_PATH)

                # relabel nodes to seqential integers, but retain dict to keep track of vertexes and celltypes etc
                int_labels =  [int(i) for i in primary_nodes] # [i for i in range(len(primary_nodes))]
                G = nx.relabel.convert_node_labels_to_integers(G, first_label=0, label_attribute='label')

                # create dataframe of node object labels and vertex ids:
                node_ids = sb.get_node_ids_3(objects, imagename, node_attr, int_labels)
                node_label_dict = dict(zip(node_ids['vertex'], node_ids.index))
                ## to do: use node attr variable in spatial_neighbour and nearest_enighbour graph construction for flexibility

            ## Loop through cell types
            barrier_df_list = []

            for cellType in cellTypes:
                minimum_paths = []
                barrier_content = []
                all_cellchains = []
                all_vertexes = []
                all_chain_lengths = []
                all_degenerate_counts = []
                all_degenerate_adjacent = []
                all_vertex_chains = []

                ## node ids of the target cell type:
                cellType_node_ids = sb.get_cellType_node_ids(node_ids, phenotyping_level, cellType)
                if calculate_positivity:
                    positive_dict = dict(zip(node_ids['vertex'], node_ids['positive']))

                print('THIS IS THE CELL TYPE NODE IDS DATAFRAME:')
                print(cellType_node_ids)

                if len(cellType_node_ids.index) > 0:

                    ## loop through all cells of the target type:
                    for v in tqdm(cellType_node_ids.vertex, ascii=True):
                        print(v)

                        if v >= 0: # Starting vertex should be between 0 to number of vertices

                            ## compute shortest paths of vertex to all other nodes in the graph:
                            shortest_paths = sb.compute_shortest_paths(G, node_ids, source=v)

                            print('THIS IS THE SHORTEST PATHS DATAFRAME:')
                            print(shortest_paths)

                            ## what is the minimum path length to an epithelial cell?
                            minpath_to_target = sb.min_path_to_cellType(shortest_paths, phenotyping_level, args.target_cell_type)
                            print("minpath to target:", minpath_to_target)

                            # minpath to target will be call as 1.6..**308 if not connected (skip these)
                            if minpath_to_target < 1000:

                                minimum_paths.append(minpath_to_target)

                                closest_epi = shortest_paths[(shortest_paths['distance'] == minpath_to_target) & (shortest_paths[phenotyping_level] == f'{target_cell_type}')]
                                print('THIS IS THE CLOSEST EPI i.e. target DATAFRAME:')
                                print(closest_epi)
                                
                                degenerate_barrier_fraction = sb.degenerate_path_content(closest_epi, shortest_paths, minpath_to_target, barrier_cells = BARRIER_TYPES, phenotyping_level=phenotyping_level)
                                degenerate_adjacent_count = sb.degenerate_adjacent_barrier(closest_epi, shortest_paths, minpath_to_target, barrier_cells = BARRIER_TYPES, phenotyping_level=phenotyping_level)

                                all_degenerate_counts.append(degenerate_barrier_fraction)
                                all_degenerate_adjacent.append(degenerate_adjacent_count)

                                if CALC_CHAIN == True:
                                    ## compute the chain of the shortest path to the epithelial cell:
                                    ## calculate the 'cell chain' along this path
                                    cellchain, vertexes = sb.followchain(shortest_paths, minpath_to_target, source_cell=cellType, source_vertex=v, target_cell=target_cell_type, phenotyping_level=phenotyping_level)

                                    all_cellchains.append(cellchain)
                                    all_vertex_chains.append(vertexes)
                                    all_chain_lengths.append(len(cellchain))
                                    n_barrier_cells = sb.count_barrier_cells(cellchain, BARRIER_TYPES)
                                    barrier_content.append(n_barrier_cells)
                                    all_vertexes.append(int(v))

                    ## construct dataframe of results:
                    barrier_df = pd.DataFrame(data=list(zip(all_chain_lengths, barrier_content, all_degenerate_counts, all_degenerate_adjacent)), columns=['chainlength', 'barrier_content', 'degenerate_barrier_fraction', 'degenerate_adjacent_fraction'])
                    barrier_df['cell_chain'] = all_cellchains
                    barrier_df['vertex_chain'] = all_vertex_chains
                    barrier_df['internal_chainlength'] = barrier_df['chainlength'] - 2
                    barrier_df['barrier_fraction'] = barrier_df['barrier_content'] / barrier_df['internal_chainlength']
                    barrier_df['imagename'] = imagename
                    barrier_df['source_cell'] = cellType
                    barrier_df['target_cell'] = f'{target_cell_type}'
                    barrier_df['vertex'] = all_vertexes

                    ## relabel vertexes with object IDs from typing tables:
                    barrier_df['object'] = barrier_df['vertex'].map(node_label_dict) ## remap object ids from vertices
                    barrier_df['object_chain'] = barrier_df['vertex_chain'].apply(lambda x: sb.relabel_vertex_chain(x, node_label_dict))

                    ## add marker positivity information:
                    if calculate_positivity == True:
                        barrier_df['positivity_chain'] = barrier_df['vertex_chain'].apply(lambda x: sb.relabel_vertex_chain(x, positive_dict))

                    ## apply lambda function to calculate barrier score with specific cell types:
                    print(barrier_df)
                    barrier_df['weighted_barrier_content'] = barrier_df.apply(lambda x: sb.weighted_barrier_count(x['cell_chain'], BARRIER_TYPES), axis=1)
                    barrier_df['binary_barrier'] = barrier_df.apply(lambda x: sb.binary_barrier_call(x['cell_chain'], BARRIER_TYPES), axis=1)
                    barrier_df['adjacent_barrier'] = barrier_df.apply(lambda x: sb.adjacent_barrier_call(x['cell_chain'], BARRIER_TYPES), axis=1)
                    barrier_df_list.append(barrier_df)

            if len(barrier_df_list) > 0:
                ## Concatenate barrier dataframes for all cell types:
                image_barrier_df = pd.concat(barrier_df_list)

                print('image_barrier_df', image_barrier_df)

                if PERMUTE_PHENOTYPES == True: 
                    spath = os.path.join(RESULTS_DIR, 'permuted', f'{imagename}_{cellType}_to_{target_cell_type}_barrier_results.csv'.replace(' ', '_'))

                elif EPI_NO_STROMA:
                    spath = os.path.join(RESULTS_DIR, 'no_stroma_epithelial', f'{imagename}_{cellType}_to_{target_cell_type}_barrier_results.csv'.replace(' ', '_'))
                elif LARGE_CLUSTERS_ONLY:
                    spath = os.path.join(RESULTS_DIR, 'large_clusters_only', f'{imagename}_{cellType}_to_{target_cell_type}_barrier_results.csv'.replace(' ', '_'))
                elif SIZE_THRESH_NO_UNCLUSTERED:
                    spath = os.path.join(RESULTS_DIR, 'size_threshold_no_unclustered', f'{imagename}_{cellType}_to_{target_cell_type}_barrier_results.csv'.replace(' ', '_'))
                else:
                    spath = os.path.join(RESULTS_DIR, 'unfiltered', f'{imagename}_{cellType}_to_{target_cell_type}_barrier_results.csv'.replace(' ', '_'))

                dirname = os.path.dirname(spath)
                if not os.path.exists(dirname):
                    os.makedirs(dirname)

                image_barrier_df.to_csv(spath, sep=OBJECT_SEP)

        else:
            warnings.warn(f"There are no typed objects for {imagename}.")
    
    else:
        warnings.warn("No Epithelial cell domains larger than the cutoff criteria. Barrier score will not be measured and no output will be produced.")



if __name__ == '__main__':

    # create argument parser
    parser = argparse.ArgumentParser(description = 'Stromal barrier measurement parameters.')
    parser.add_argument('--adjacency_data_path', help='folder containing adjacency list files in csv format for neighbouRhood graphs')
    parser.add_argument('--barrier_types', nargs='+', help='Cell types to assign as barrier cells e.g. Myofibroblasts. Multiple arguments accepted e.g. --barrier_types Myofibroblasts Fibroblasts.')
    parser.add_argument('--calc_chain', type = bool, help='Calculate the chain of cell objects from the starting cell to the end cell.', default=True) # remove this argument as we always want to calculate the chain
    parser.add_argument('--domain_size_cutoff', type = int, help='connectivity type for cell spatial graph construction', default= 2000)
    parser.add_argument('--graph_type', help='connectivity type for cell spatial graph construction')
    parser.add_argument('--imagename', type = str, help='Name of image in cell objects dataframe')
    parser.add_argument('--neighbourhood_radius', type = int, help='Dilation used to determine cell neighbours in neighbouRhood graph')
    parser.add_argument('--neighbours', type = int, help='number of neighbours for nearest neighbour graph')
    parser.add_argument('--objects_path', help='/path/to/cell objects dataframe.')
    parser.add_argument('--objects_sep', help='Objects file delimiter.')
    parser.add_argument('--panel', help='IMC panel name.')
    parser.add_argument('--permutation_region', help='Domain in which to permute cells. e.g. "tumour" or "stroma. Depends on this information being available in the cell objects table under column "region".')
    parser.add_argument('--permute_phenotypes', type = bool, help='Randomly permute cell phenotypes in a given domain.', default=False)
    parser.add_argument('--phenotyping_level', help='Designation of the objects table column to use to determine phenotypes e.g. majorType or cellType, but depends can be other depending on columns in objects.csv')
    parser.add_argument('--radius', type = float, help='radius for spatial neighbours graph')
    parser.add_argument('--root_out', help='Root output directory for saving.')
    parser.add_argument('--source_cell_type', help='source cell type for the shortest path calculation', default='CD8 T cells')
    parser.add_argument('--target_cell_type', help='target cell type for the shortest path calculation', default='Epithelial cells')
    args = parser.parse_args()

    # pass command line arguments to main:
    main(args)