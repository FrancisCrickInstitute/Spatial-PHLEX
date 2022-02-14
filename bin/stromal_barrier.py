#!/usr/bin/env python

import pandas as pd
# import cugraph as cg
import networkx as nx
import numpy as np
from tqdm import *
from collections import Counter
import os, sys
import sb_functions as sb 

def main(args):

    # config:
    GRAPH_TYPE = 'neighbouRhood' # 'nearest_neighbour' or 'spatial_neighbours', 'neighbouRhood'
    NEIGHBOURS = 10
    RADIUS = 5
    NEIGHBOURHOOD_RADIUS = 5
    ADJACENCY_DATA_PATH = args[0]
    ROOT_OUT = args[1]
    OBJECTS_PATH = args[2]
    PANEL = args[3] #'p2'
    imagename = os.path.split(ADJACENCY_DATA_PATH)[1].replace('.txt', '')
    CALC_CHAIN = True
    BARRIER_TYPES = ['Myofibroblasts']
    phenotyping_level = 'cellType' # 'cellType', 'majorType', 'Positive', 'cellClass'
    
    # Define immune cell subtypes to measure the 'barrier' for:
    cellTypes = ['CD4 T cells', 'CD4 T cells', 'CD57+ CD4 T cells', 'Naive CD4 T cells',
       'Leukocytes - Other', 'CD4+ Myeloid cells', 'Cytotoxic CD8 T cells',
       'CD4 Tcm', 'Tregs', 'CD8 Trm', 'CD8 T cells', 'CD57+ CD8 Trm',
       'CD8 Exhausted TDT', 'Naive CD8 T cells', 'T cells DP',
       'Cytotoxic CD4 T cells', 'T cells - Other']


    if GRAPH_TYPE == 'spatial_neighbours':
        RESULTS_DIR = os.path.join(ROOT_OUT, GRAPH_TYPE, f'radius_{RADIUS}', 'barrier', '_'.join(BARRIER_TYPES))
    elif GRAPH_TYPE == 'nearest_neighbour':
        RESULTS_DIR = os.path.join(ROOT_OUT, GRAPH_TYPE, f'n_neigh_{NEIGHBOURS}', 'barrier', '_'.join(BARRIER_TYPES))
    elif GRAPH_TYPE == 'neighbouRhood':
        RESULTS_DIR = os.path.join(ROOT_OUT, GRAPH_TYPE, f'neighbouRhood_radius_{NEIGHBOURHOOD_RADIUS}', 'barrier', '_'.join(BARRIER_TYPES))
    else:
        raise ValueError("GRAPH_TYPE not recognised.")

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR, exist_ok=True)

    # Read in necessary files:
    objects = pd.read_csv(OBJECTS_PATH, sep='\t', encoding='latin1')

    if imagename in objects['imagename'].unique():

        objects = sb.assign_cell_categories(objects, typing='new')
        print(objects)

        # SPECIFY IMC MARKERS:
        if PANEL == 'p1':
            markers = ['CD39', 'CD4', 'LAG3', 'casp3', 'CD31', 'CD27', 'CXCR6', 'PDL1', 'CXCR4', 'TCF1', 'CD103', 'GATA3', 'FAP1', 'ICOS', 'Vimentin', 'CD25', 'PD1', 'CD3', 'GITR', 'alphaSMA', 'CD8a', 'CTLA4', 'B2M', 'TIM3', 'Ki67', 'CD57', 'CD45RA', 'CXCL12', 'panCK', 'CCR7', 'Collagen', 'pSTAT1', 'GZMB', 'CD45', 'FOXP3']
        if PANEL == 'p2':
            markers = ['CD163', 'CLEC9a', 'KIR2DL3', 'VISTA', 'CD56', 'TIM3', 'CD45', 'CD16', 'MHCII', 'panactin', 'TCRd', 'PD1', 'CD79a', 'CD66b', 'CD14', 'MPO', 'CAIX', 'CD38', 'MCT4', 'CD206', 'IDO', 'CD68', 'PDL1', 'CD20', 'CD11c', 'CD4', 'CD8a', 'CD31', 'LAG3', 'CD103', 'CD11b', 'CD3', 'GZMB', 'panCK', 'CD73']


        ## CREATE SPATIAL GRAPH
        if GRAPH_TYPE == 'spatial_neighbours':
            spg = sb.compute_spatial_graph(objects, imagename, markers, radius=RADIUS)
            G = nx.convert_matrix.from_scipy_sparse_matrix(spg.obsp['spatial_connectivities'])
            node_ids = sb.get_graph_node_ids(spg)
        elif GRAPH_TYPE == 'nearest_neighbour':
            spg = sb.compute_nn_graph(objects, imagename, markers, n_neighs=NEIGHBOURS)
            G = nx.convert_matrix.from_scipy_sparse_matrix(spg.obsp['spatial_connectivities'])
            node_ids = sb.get_graph_node_ids(spg)
        elif GRAPH_TYPE == 'neighbouRhood':

            # specify attributes to attach to nodes:
            node_attr = ['majorType', 'cellType', 'positive','cellClass']
            
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
            cellType_node_ids = sb.get_cellType_node_ids(node_ids, 'cellType', cellType)
            
            if len(cellType_node_ids.index) > 0:

                ## loop through all cells of the target type:
                for v in tqdm(cellType_node_ids.vertex, ascii=True):
                    print(v)

                    if v >= 0: # Starting vertex should be between 0 to number of vertices

                        ## compute shortest paths of vertex to all other nodes in the graph:
                        shortest_paths = sb.compute_shortest_paths(G, node_ids, source=v)
                    
                        ## what is the minimum path length to an epithelial cell?
                        minpath_to_epi = sb.min_path_to_epithelial(shortest_paths)
                        print("minpath to epi:", minpath_to_epi)

                        # minpath to epi will be call as 1.6..**308 if not connected (skip these)
                        if minpath_to_epi < 1000:

                            minimum_paths.append(minpath_to_epi)

                            closest_epi = shortest_paths[(shortest_paths['distance'] == minpath_to_epi) & (shortest_paths['cellType'] == 'Epithelial cells')]
                            degenerate_barrier_fraction = sb.degenerate_path_content(closest_epi, shortest_paths, minpath_to_epi, barrier_cells = BARRIER_TYPES)
                            degenerate_adjacent_count = sb.degenerate_adjacent_barrier(closest_epi, shortest_paths, minpath_to_epi, barrier_cells = BARRIER_TYPES)

                            all_degenerate_counts.append(degenerate_barrier_fraction)
                            all_degenerate_adjacent.append(degenerate_adjacent_count)

                            if CALC_CHAIN == True:
                                ## calculate the 'cell chain' along this path
                                cellchain, vertexes = sb.followchain(shortest_paths, minpath_to_epi, source_cell=cellType, source_vertex=v, phenotype_level=phenotyping_level)

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
            barrier_df['target_cell'] = 'Epithelial cells'
            barrier_df['vertex'] = all_vertexes
            barrier_df_list.append(barrier_df)

        ## Concatenate barrier dataframes for all cell types:
        image_barrier_df = pd.concat(barrier_df_list)

        ## relabel vertexes with object IDs from typing tables:
        image_barrier_df['object'] = image_barrier_df['vertex'].map(node_label_dict) ## remap object ids from vertices
        image_barrier_df['object_chain'] = image_barrier_df['vertex_chain'].apply(lambda x: sb.relabel_vertex_chain(x, node_label_dict))

        ## add marker positivity information:
        positive_dict = dict(zip(node_ids['vertex'], node_ids['positive']))
        image_barrier_df['positivity_chain'] = image_barrier_df['vertex_chain'].apply(lambda x: sb.relabel_vertex_chain(x, positive_dict))

        ## apply lambda function to calculate barrier score with specific cell types:
        image_barrier_df['weighted_barrier_content'] = image_barrier_df.apply(lambda x: sb.weighted_barrier_count(x['cell_chain'], BARRIER_TYPES), axis=1)
        image_barrier_df['binary_barrier'] = image_barrier_df.apply(lambda x: sb.binary_barrier_call(x['cell_chain'], BARRIER_TYPES), axis=1)
        image_barrier_df['adjacent_barrier'] = image_barrier_df.apply(lambda x: sb.adjacent_barrier_call(x['cell_chain'], BARRIER_TYPES), axis=1)

        spath = os.path.join(RESULTS_DIR, f'{imagename}_barrier_results.csv')
        image_barrier_df.to_csv(spath)
        
    else:
        print(f"There are no typed objects for {imagename}.")


if __name__ == '__main__':
    main(sys.argv[1:])