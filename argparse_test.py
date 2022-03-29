#! /usr/bin/env python

import argparse
import os

def main(args):

    print(args.graph_type)

    GRAPH_TYPE = args.graph_type
    NEIGHBOURS = args.neighbours

    print(GRAPH_TYPE)
    print('neighbours', NEIGHBOURS)
    print('objects sep', args.objects_sep)

    print(args.neighbours)

    print(args.barrier_types)
    print(type(args.barrier_types))

    if type(args.barrier_types) == list:
        print('its a list')
        BARRIER_TYPES = args.barrier_types #['Myofibroblasts']
    else:
        BARRIER_TYPES = [args.barrier_types]

    print(BARRIER_TYPES)
    print(len(BARRIER_TYPES))

    if GRAPH_TYPE == 'nearest_neighbour':
        RESULTS_DIR = os.path.join('./', GRAPH_TYPE, f'n_neigh_{NEIGHBOURS}', 'barrier', '_'.join(BARRIER_TYPES))
    print(RESULTS_DIR)


if __name__ == '__main__':

    # create argument parser
    parser = argparse.ArgumentParser(description = 'Stromal barrier measurement parameters.')
    parser.add_argument('--graph_type', help='connectivity type for cell spatial graph construction')
    parser.add_argument('--neighbourhood_radius', type = int, help='Dilation used to determine cell neighbours in neighbouRhood graph')
    parser.add_argument('--adjacency_data_path', help='folder containing adjacency list files in csv format for neighbouRhood graphs')
    parser.add_argument('--radius', type = float, help='radius for spatial neighbours graph')
    parser.add_argument('--neighbours', type = int, help='number of neighbours for nearest neighbour graph')
    parser.add_argument('--root_out', help='Root output directory for saving.')
    parser.add_argument('--objects_path', help='/path/to/cell objects dataframe.')
    parser.add_argument('--objects_sep', help='Objects file delimiter.')
    parser.add_argument('--panel', help='IMC panel name.')
    parser.add_argument('--calc_chain', type = bool, help='Calculate the chain of cell objects from the starting cell to the end cell.')
    parser.add_argument('--imagename', help='Name of image in cell objects dataframe')
    parser.add_argument('--permute_phenotypes', type = bool, help='Randomly permute cell phenotypes in a given domain.', default=False)
    parser.add_argument('--permutation_region', help='Domain in which to permute cells. e.g. "tumour" or "stroma. Depends on this information being available in the cell objects table under column "region".')
    parser.add_argument('--barrier_types', type = str, nargs='+', help='Cell types to assign as barrier cells e.g. Myofibroblasts. Multiple arguments accepted e.g. --barrier_types Myofibroblasts Fibroblasts.')
    parser.add_argument('--phenotyping_level', help='Designation of the objects table column to use to determine phenotypes e.g. majorType or cellType, but depends can be other depending on columns in objects.csv')
    args = parser.parse_args()
    # pass command line arguments to main:
    main(args)