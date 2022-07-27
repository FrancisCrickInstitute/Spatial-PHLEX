#!/usr/bin/env python

import pandas as pd
import networkx as nx
import numpy as np
import glob, os, sys

def df_select(df, cat, val):
    return df[df[cat] == val]

def cell_adjacency(object_number, neighbours):
    nlist = neighbours['Second Object Number'].tolist()
    nlist.insert(0, object_number)
    adj_str = to_adj_string(nlist)
    return adj_str

def to_adj_string(adj_list):
    adj_str = ''
    for elem in adj_list:
        adj_str += str(elem) + ' '
    adj_list = adj_str[:-1] # remove trailing space
    return adj_str

def get_imagename(filepath):
    splitpath = filepath.split('/')
    print(splitpath)
    imagename = f'{splitpath[-3]}-{splitpath[-2]}'
    return imagename

def write_adj_list(adj_list, imagename, outdir):
    path = os.path.join(outdir, f'{imagename}.txt')
    with open(path, 'w') as f:
        for item in adj_list:
            f.write("%s\n" % item)

def main(args):
    '''
    Script to convert cellprofiler neighbourhood tables to networkx-compatible adjacency list.
    # args[0] = path to csv neigbourhood file
    # args[1] = Module Number  of the Cellprofiler neighbouRhood module number in the CellProfiler pipeline.
    # 
    # args[3] = root output directory.
    '''
    inputpath = args[0]
    neighbour_module = int(args[1]) 
    root_out = args[2]
    
    outdir = os.path.join(root_out, f'cp_module_{neighbour_module}')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    imagename = get_imagename(inputpath)
    print(imagename)
    cp_nhood = pd.read_csv(inputpath)
    
    # subset 5-micron dilation neighbours from df:
    all_neighbours = df_select(cp_nhood, 'Relationship', 'Neighbors')
    print(all_neighbours)
    neighbours_5 = df_select(all_neighbours, 'Module Number', neighbour_module)
    print(neighbours_5)
    
    # loop through abd extract neighborus for every object id:
    adjacencylist = []
    for objectNumber in neighbours_5['First Object Number'].unique():
        object_neighbours = df_select(neighbours_5, 'First Object Number', objectNumber)
        adj = cell_adjacency(objectNumber, object_neighbours)
        adjacencylist.append(adj)
        
    write_adj_list(adj_list=adjacencylist, imagename=imagename, outdir=outdir)

if __name__ == '__main__':
    main(sys.argv[1:])