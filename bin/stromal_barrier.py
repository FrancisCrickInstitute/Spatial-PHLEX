#!/usr/bin/env python

import pandas as pd

# import cugraph as cg
import networkx as nx
import numpy as np
from tqdm import *
from collections import Counter
import os, sys
import sb_functions as sb
import sb_scores
import argparse
import warnings


def main(args):

    if isinstance(args.barrier_types, list):
        barrier_types = args.barrier_types  # ['Myofibroblasts']
    else:
        barrier_types = [args.barrier_types]

    phenotyping_column = args.phenotyping_column

    # Read in necessary files:
    objects = pd.read_csv(
        args.objects_path, sep=args.objects_sep, encoding="latin1", quoting=0
    )
    objects[args.image_id_col] = objects[args.image_id_col].astype(str)

    # Measure barrier to cell cluster large than an area cutoff:
    if args.clustered_barrier:
        objects = sb.preprocess_clustered_barrier(
            objects,
            args.phenotyping_column,
            args.target_cell_type,
            cluster_size_cutoff=args.cluster_size_cutoff,
        )

    # objects = objects.sample(n=1000, random_state=1)

    ## after filtering only proceed if there are epithelial cells that pass the criteria, else raise warning:
    if len(objects[objects[phenotyping_column] == args.target_cell_type].index) > 0:
        assert (
            args.imagename in objects[args.image_id_col].unique()
        ), f"Image {args.imagename} not found in objects table."

        if args.permute_phenotypes:
            objects = sb.permute_phenotypes(
                objects, in_place=True, region=args.permutation_region
            )

        ## CREATE SPATIAL GRAPH
        if args.graph_type == "spatial_neighbours":
            G, node_ids = sb.make_cell_graph(
                objects,
                args.imagename,
                graph_type="spatial_neighbours",
                x_id=args.x_id,
                y_id=args.y_id,
                phenotyping_column=args.phenotyping_column,
                radius=args.radius,
            )
            results_dir = os.path.join(
                args.root_out,
                args.graph_type,
                f"radius_{args.radius}",
                "_".join(barrier_types),
            )

        elif args.graph_type == "nearest_neighbour":
            G, node_ids = sb.make_cell_graph(
                objects,
                args.imagename,
                graph_type="nearest_neighbour",
                x_id=args.x_id,
                y_id=args.y_id,
                phenotyping_column=args.phenotyping_column,
                n_neighs=args.neighbours,
            )
            results_dir = os.path.join(
                args.root_out,
                args.graph_type,
                f"n_neigh_{args.neighbours}",
                "_".join(barrier_types),
            )

        else:
            raise ValueError("args.graph_type not recognised.")

        if not os.path.exists(results_dir):
            os.makedirs(results_dir, exist_ok=True)

        ## Measure barrier scores for the source cell type
        barrier_df = sb_scores.measure_barrier(
            G, node_ids, phenotyping_column, args.imagename,
            source_cell_type=args.source_cell_type,
            target_cell_type=args.target_cell_type,
            barrier_types=barrier_types,
            image_id_col=args.image_id_col
        )

        spath = os.path.join(
                results_dir,
                f"{args.imagename}_{args.source_cell_type}_to_{args.target_cell_type}_barrier_results.csv".replace(
                    " ", "_"
                ),
        )

        barrier_df['cluster_area_cutoff'] = args.cluster_size_cutoff
        barrier_df.to_csv(spath, index=False)

    else:
        warnings.warn(
            f"No {args.target_cell_type} domains larger than the cutoff criteria. Barrier score will not be measured and no output will be produced."
        )


if __name__ == "__main__":
    # create argument parser
    parser = argparse.ArgumentParser(
        description="Stromal barrier measurement parameters."
    )
    parser.add_argument(
        "--objects_path",
        help="/path/to/cell objects dataframe.",
    )
    parser.add_argument("--objects_sep", help="Objects file delimiter.")
    parser.add_argument(
        "--barrier_types",
        nargs="+",
        help="Cell types to assign as barrier cells e.g. Myofibroblasts. Multiple arguments accepted e.g. --barrier_types Myofibroblasts Fibroblasts.",
    )
    parser.add_argument(
        "--cluster_size_cutoff",
        type=int,
        help="Domain size cutoff in area units to measure barrier scores for.",
        default=2000,
    )
    parser.add_argument(
        "--x_id",
        help="Name of column in cell objects dataframe containing x coordinates.",
        default="x",
    )
    parser.add_argument(
        "--y_id",
        help="Name of column in cell objects dataframe containing y coordinates.",
        default="y",
    )
    parser.add_argument(
        "--image_id_col",
        help="Name of the column in the input data containing different image IDs.",
        default="image_id",
    )
    parser.add_argument(
        "--graph_type", help="connectivity type for cell spatial graph construction"
    )
    parser.add_argument(
        "--imagename", type=str, help="Name of image in cell objects dataframe"
    )
    parser.add_argument(
        "--neighbours",
        type=int, default=6,
        help="number of neighbours for nearest neighbour graph",
    )
    parser.add_argument(
        "--clustered_barrier",
        type=bool,
        help="Calculate barrier score for clustered cells only.",
        default=True,
    )
    parser.add_argument(
        "--permutation_region",
        help='Domain in which to permute cells. e.g. "tumour" or "stroma. Depends on this information being available in the cell objects table under column "region".',
    )
    parser.add_argument(
        "--permute_phenotypes",
        type=bool,
        help="Randomly permute cell phenotypes in a given domain.",
        default=False,
    )
    parser.add_argument(
        "--phenotyping_column",
        help="Designation of the objects table column to use to determine phenotypes e.g. majorType or args.source_cell_type, but depends can be other depending on columns in objects.csv",
    )
    parser.add_argument(
        "--radius", type=float, help="radius for spatial neighbours graph"
    )
    parser.add_argument("--root_out", help="Root output directory for saving.")
    parser.add_argument(
        "--source_cell_type",
        help="source cell type for the shortest path calculation",
        default="CD8 T cells",
    )
    parser.add_argument(
        "--target_cell_type",
        help="target cell type for the shortest path calculation",
        default="Epithelial cells",
    )
    args = parser.parse_args()

    # pass command line arguments to main:
    main(args)
