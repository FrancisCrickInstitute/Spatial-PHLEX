#!/usr/bin/env python

from centrality_functions import df_select
from centrality_functions import assign_cell_categories, compute_cohort_centralities, compute_spatial_graph, compute_nn_graph
import pandas as pd
import cugraph as cg
import networkx as nx
import numpy as np
from tqdm import *
from collections import Counter
import os, sys

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

def followchain(shortest_paths, minpath_to_target, source_cell, target_cell = 'Epithelial cells', phenotype_level='cellType'):
    '''
    The min_path_to_target is e.g. 3 jumps from the source to the target cell; 
    there may be more than one path of the same length.
    
    Todo: return vertex ids for subsequent referencing of e.g. positive marker expression'''
    
    # print(minpath_to_target)
    if minpath_to_target == 1:
        vertexes = []
        chaincells = [target_cell, source_cell]
        # print(chaincells)
        return chaincells, vertexes
    
    else:
        sorted_shortest_paths = shortest_paths.sort_values(by='distance')
        sorted_shortest_paths = sorted_shortest_paths[(sorted_shortest_paths['distance'] == minpath_to_target) & (sorted_shortest_paths['cellType'] == 'Epithelial cells')]

        if sorted_shortest_paths['predecessor'].max() >= 0: # possible to have a df with target cells that are disconnected
            chaincells = [target_cell] # store cells in chain between starting cell and destination 
            vertexes = []

            predecessors = get_predecessor(sorted_shortest_paths, minpath_to_target) # get list of initial predecessor cells
            # print('predecessor_0: ', predecessors)
            
            p = predecessors[0]

            for j in range(minpath_to_target-1):
            # take first predecessor if more than one:
            
                # print('p:', p)
                if p != -1: ## proceed only if the predecessor is not the source cell
                    predecessor_nodeType = get_predecessor_nodeType(shortest_paths, p, phenotype_level=phenotype_level)
                    chaincells.append(predecessor_nodeType)
                    # print(p, predecessor_nodeType)
                    
                    # print('obtaining chain predecessor for p')
                    p = get_chain_predecessor(shortest_paths, p)
                    predecessor_nodeType = get_predecessor_nodeType(shortest_paths, p, phenotype_level=phenotype_level)
                    # print(p, predecessor_nodeType)
                    chaincells.append(predecessor_nodeType)
                    vertexes.append(p)
                else:
                    chaincells.append(source_cell)
            return chaincells, vertexes
                


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
    # print(df[df['vertex'] == predecessor]) #[phenotype_level].values)
    nodetype = df[df['vertex'] == predecessor].iloc[0][phenotype_level]
#     nodetype = df[df['vertex'] == predecessor][phenotype_level].values.tolist()
    # print('the nodetype of the predecesor: ', nodetype)
    return nodetype

def get_graph_node_ids(spg):
    """
    Produce dataframe of graph node ids from squidpy spatial graph.
    """
#     node_ids = spg.obs['cellType']
    node_ids = pd.DataFrame(spg.obs[['majorType', 'cellType', 'positive', 'cellClass']])
    node_ids['vertex'] = np.arange(0,len(node_ids.index))
    return node_ids

def get_node_ids_2(objects, imagename, attributes):
    """
    Produce dataframe of graph node ids from objects table and specified attributes.
    """
    image_objects = objects[objects['imagename'] == imagename].reset_index()
    image_objects = image_objects[attributes]
#     image_objects = image_objects.rename(columns={'object':'vertex'})
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

## loop over epithelial cells with minpath_to_epi steps in chain:

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
            print('DEGENERATE PATH CONTENT', degenerate_path_content)

            print(set(barrier_cells).intersection(set(degenerate_path_content['cellType'])))
            if set(barrier_cells).intersection(set(degenerate_path_content['cellType'])) != set():
                print('CLOSEST EPI:', closest_epi)
                
                myofibroblast_fraction_degen = degenerate_path_content['cellType'].value_counts(normalize=True)[barrier_cells].sum()
                print('myofibroblast_fraction_degen:', myofibroblast_fraction_degen)
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
            print('DEGENERATE PATH CONTENT', degenerate_path_content)

            print('set intersection:', set(barrier_cells).intersection(set(degenerate_path_content['cellType'])))
            adjacent_cells_df = degenerate_path_content[degenerate_path_content['distance'] == (minpath_to_epi-1)]

            if set(barrier_cells).intersection(set(adjacent_cells_df['cellType'])) != set():
                
                print('ADJACENT CELLS DF:', adjacent_cells_df)
                print('BARRIER CELLS:', barrier_cells)
                mean_adjacent_barrier = adjacent_cells_df['cellType'].value_counts(normalize=True)[barrier_cells].sum() #/ len(adjacent_cells_df)
                print(mean_adjacent_barrier)
            else:
                mean_adjacent_barrier = 0

    else:
        mean_adjacent_barrier = 0
    return mean_adjacent_barrier
## todo: binary degenrate barrier i.e. get 1 on each path if a myofibroblast present; degenerate adjacent


def main(args):

    # config:
    GRAPH_TYPE = 'neighbouRhood' # 'nearest_neighbour' or 'spatial_neighbours', 'neighbouRhood'
    NEIGHBOURS = 10
    RADIUS = 5
    neighbouRhood_radius = 5
    ADJACENCY_DATA_PATH = args[0]
    ROOT_OUT = args[1]
    OBJECTS_PATH = args[2]
    PANEL = args[3] #'p2'

    imagename = os.path.split(ADJACENCY_DATA_PATH)[1].replace('.txt', '')
    
    CALC_CHAIN = True

    BARRIER_TYPES = ['Myofibroblasts']
    phenotyping_level = 'cellType' # 'cellType', 'majorType', 'Positive', 'cellClass'
    
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
        RESULTS_DIR = os.path.join(ROOT_OUT, GRAPH_TYPE, f'neighbouRhood_radius_{neighbouRhood_radius}', 'barrier', '_'.join(BARRIER_TYPES))
    else:
        raise ValueError("GRAPH_TYPE not recognised.")

    if not os.path.exists(RESULTS_DIR):
        os.makedirs(RESULTS_DIR, exist_ok=True)

    # Read in necessary files:
    objects = pd.read_csv(OBJECTS_PATH, sep='\t', encoding='latin1')

    if imagename in objects['imagename'].unique():

        # reference_celldata = pd.read_csv(f'/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/features/{PANEL}_MeanIntensity_final.csv')
        # reference_celldata = reference_celldata.rename(columns={'ObjectNumber': 'object'})
        # objects = pd.merge(objects, reference_celldata, how='left', on=['imagename', 'object'])
        objects = assign_cell_categories(objects, typing='new')
        print(objects)

        # SPECIFY IMC MARKERS:
        if PANEL == 'p1':
            markers = ['CD39', 'CD4', 'LAG3', 'casp3', 'CD31', 'CD27', 'CXCR6', 'PDL1', 'CXCR4', 'TCF1', 'CD103', 'GATA3', 'FAP1', 'ICOS', 'Vimentin', 'CD25', 'PD1', 'CD3', 'GITR', 'alphaSMA', 'CD8a', 'CTLA4', 'B2M', 'TIM3', 'Ki67', 'CD57', 'CD45RA', 'CXCL12', 'panCK', 'CCR7', 'Collagen', 'pSTAT1', 'GZMB', 'CD45', 'FOXP3']
        if PANEL == 'p2':
            markers = ['CD163', 'CLEC9a', 'KIR2DL3', 'VISTA', 'CD56', 'TIM3', 'CD45', 'CD16', 'MHCII', 'panactin', 'TCRd', 'PD1', 'CD79a', 'CD66b', 'CD14', 'MPO', 'CAIX', 'CD38', 'MCT4', 'CD206', 'IDO', 'CD68', 'PDL1', 'CD20', 'CD11c', 'CD4', 'CD8a', 'CD31', 'LAG3', 'CD103', 'CD11b', 'CD3', 'GZMB', 'panCK', 'CD73']


        ## CREATE SPATIAL GRAPH
        if GRAPH_TYPE == 'spatial_neighbours':
            spg = compute_spatial_graph(objects, imagename, markers, radius=RADIUS)
            G = nx.convert_matrix.from_scipy_sparse_matrix(spg.obsp['spatial_connectivities'])
            node_ids = get_graph_node_ids(spg)
        elif GRAPH_TYPE == 'nearest_neighbour':
            spg = compute_nn_graph(objects, imagename, markers, n_neighs=NEIGHBOURS)
            G = nx.convert_matrix.from_scipy_sparse_matrix(spg.obsp['spatial_connectivities'])
            node_ids = get_graph_node_ids(spg)
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
            node_label_dict = nx.get_node_attributes(G,'label')
            node_ids = get_node_ids_3(objects, imagename, node_attr, int_labels)
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

            ## node ids of the target cell type:
            cellType_node_ids = get_cellType_node_ids(node_ids, 'cellType', cellType)
        #     cellType_node_ids.reset_index(inplace=True)
            
            if len(cellType_node_ids.index) > 0:
                print('\ncellType_node_ids:')
                print(cellType_node_ids)

                ## loop through all cells of the target type:
                for v in tqdm(cellType_node_ids.vertex, ascii=True):
                    print(v)

                    if v >= 0: # Starting vertex should be between 0 to number of vertices
        #                 print('Starting vertex:', cellType_node_ids.vertex[i])
                        #print('starting vertex:', v)

                        ## compute shortest paths of vertex to all other nodes in the graph:
                        shortest_paths = compute_shortest_paths(G, node_ids, source=v)
                    
                        ## what is the minimum path length to an epithelial cell?
                        minpath_to_epi = min_path_to_epithelial(shortest_paths)
                        print("minpath to epi:", minpath_to_epi)

                        # minpath to epi will be call as 1.6..**308 if not connected (skip these)
                        if minpath_to_epi < 1000:

                            minimum_paths.append(minpath_to_epi)

                            closest_epi = shortest_paths[(shortest_paths['distance'] == minpath_to_epi) & (shortest_paths['cellType'] == 'Epithelial cells')]
                            degenerate_barrier_fraction = degenerate_path_content(closest_epi, shortest_paths, minpath_to_epi, barrier_cells = BARRIER_TYPES)
                            degenerate_adjacent_count = degenerate_adjacent_barrier(closest_epi, shortest_paths, minpath_to_epi, barrier_cells = BARRIER_TYPES)

                            all_degenerate_counts.append(degenerate_barrier_fraction)
                            all_degenerate_adjacent.append(degenerate_adjacent_count)

                            if CALC_CHAIN == True:
                                ## calculate the 'cell chain' along this path
                                cellchain, vertexes = followchain(shortest_paths, minpath_to_epi, source_cell=cellType, phenotype_level=phenotyping_level)

                                all_cellchains.append(cellchain)
                                all_chain_lengths.append(len(cellchain))
                                n_barrier_cells = count_barrier_cells(cellchain, BARRIER_TYPES)
                                barrier_content.append(n_barrier_cells)
                                
                                all_vertexes.append(int(v))

                if len(all_cellchains) > 0:
                    mod_chains = []
                    for chain in all_cellchains:
                        if chain != None:
                            mod_chains.append(chain[1:-1])

                    mod_chain_df = pd.DataFrame(mod_chains).T

            ## construct dataframe of results:
            barrier_df = pd.DataFrame(data=list(zip(all_chain_lengths, barrier_content, all_degenerate_counts, all_degenerate_adjacent)), columns=['chainlength', 'barrier_content', 'degenerate_barrier_fraction', 'degenerate_adjacent_fraction'])
            barrier_df['cell_chain'] = all_cellchains
            barrier_df['internal_chainlength'] = barrier_df['chainlength'] - 2
            barrier_df['barrier_fraction'] = barrier_df['barrier_content'] / barrier_df['internal_chainlength']
            barrier_df['imagename'] = imagename
            barrier_df['source_cell'] = cellType
            barrier_df['target_cell'] = 'Epithelial cells'
            barrier_df['vertex'] = all_vertexes
            barrier_df_list.append(barrier_df)

        ## Save output:
        image_barrier_df = pd.concat(barrier_df_list)
        image_barrier_df['object'] = image_barrier_df['vertex'].map(node_label_dict)

        ## apply lambda function to calculate barrier score with specific cell types:
        image_barrier_df['weighted_barrier_content'] = image_barrier_df.apply(lambda x: weighted_barrier_count(x['cell_chain'], BARRIER_TYPES), axis=1)
        image_barrier_df['binary_barrier'] = image_barrier_df.apply(lambda x: binary_barrier_call(x['cell_chain'], BARRIER_TYPES), axis=1)
        image_barrier_df['adjacent_barrier'] = image_barrier_df.apply(lambda x: adjacent_barrier_call(x['cell_chain'], BARRIER_TYPES), axis=1)

        spath = os.path.join(RESULTS_DIR, f'{imagename}_barrier_results.csv')
        image_barrier_df.to_csv(spath)
    else:
        print(f"There are no typed objects for {imagename}.")

if __name__ == '__main__':
    main(sys.argv[1:])