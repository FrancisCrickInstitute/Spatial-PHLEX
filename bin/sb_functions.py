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
from squidpy.gr._utils import (
    _assert_categorical_obs,
    _assert_non_empty_sequence,
    _get_valid_values,
)
from tqdm import *


def make_legible(string):
    string = string.replace("&", "&\n")
    return string


def df_select(df, cat, val):
    return df[df[cat] == val]

def preprocess_clustered_barrier(objects, 
                            phenotyping_column, 
                            target_cell_type, 
                            cluster_size_cutoff = 2000, 
                            remove_target_cells_in_domain = False,
                            domain_to_restrict = "Distal Stroma"):

    # Alter unclustered
    objects.loc[
        (objects[f"{target_cell_type}_spatial_cluster_id"] == -1)
        & (objects[phenotyping_column] == f"{target_cell_type}"),
        phenotyping_column,
    ] = f"Unclustered {target_cell_type}"
    # alter those in small clusters:
    objects.loc[
        (objects[f"{target_cell_type}_cluster_area"] < cluster_size_cutoff)
        & (objects[phenotyping_column] == f"{target_cell_type}"),
        phenotyping_column,
    ] = f"Unclustered {target_cell_type}"

    if remove_target_cells_in_domain:
        objects.loc[
            (objects["domain"] == domain_to_restrict)
            & (objects[phenotyping_column] == f"{target_cell_type}"),
            phenotyping_column,
        ] = "Unassigned"

    return objects


def compute_cohort_centralities(
    objects,
    markers,
    x_id="centerX",
    y_id="centerY",
    phenotyping_column="cellType",
    radius=35,
):
    imagenames = objects["imagename"].unique()

    cohort_centralities = []
    for i, imagename in enumerate(imagenames):
        # annotate adata with phenotyping info
        image_objects = objects[objects["imagename"] == imagename].astype("category")
        image_cTypes = image_objects[phenotyping_column].dropna().unique()

        if (
            len(image_cTypes) > 0
        ):  # only proceed if there are cells of the specified phenotype level in the image
            features = image_objects[markers].values

            cellTypes = image_objects[phenotyping_column].values
            coords = np.asarray(
                list(zip(image_objects[x_id].values, image_objects[y_id].values))
            )

            # create anndata object:
            adata = AnnData(
                features,
                obsm={"spatial": coords},
                obs={f"{phenotyping_column}": cellTypes},
            )
            # print(adata)

            ## derive spatial neighbours graph:
            sq.gr.spatial_neighbors(adata, radius=radius)

            ## compute centrality scores for image
            sq.gr.centrality_scores(
                adata,
                cluster_key=phenotyping_column,
            )

            ## obtain centrality scores
            image_centralities = adata.uns[f"{phenotyping_column}_centrality_scores"]
            image_centralities["imagename"] = imagename

            cohort_centralities.append(image_centralities)

    cohort_centralities = pd.concat(cohort_centralities)

    return cohort_centralities

def make_cell_graph(
    objects,
    imagename,
    x_id="centerX",
    y_id="centerY",
    phenotyping_column="cellType",
    graph_type="spatial_neighbours",
    radius=35,
    n_neighs=6,
):
    """Create a spatial graph from a dataframe of cell objects with squidpy. Uses random features to create the anndata object.

    Args:
        objects (pd.Dataframe): Cell objects dataframe.
        imagename (str): ID of the image to be processed in the objects dataframe.
        x_id (str, optional): column header of x coordinate. Defaults to "centerX".
        y_id (str, optional): column header og y coordinate. Defaults to "centerY".
        phenotyping_column (str, optional): column header of cell phenotype information. Defaults to "cellType".
        graph_type (str, optional): Method of graph construction to use. Defaults to "spatial_neighbours".
        radius (int, optional): Radius cutoff for spatial neighbours graph construction. Defaults to 35.
        n_neighs (int, optional): Number of nearest neighbours for nearest_neighbours graph construction. Defaults to 6.

    Raises:
        ValueError: Unrecognised graph construction methods.

    Returns:
        Anndata: Anndata object with cell spatial graph produced by graph_type method as attribute.
    """    
    # annotate adata with phenotyping info
    image_objects = objects[objects["imagename"] == imagename].astype("category")
    image_cTypes = image_objects[phenotyping_column].dropna().unique()

    if (
        len(image_cTypes) > 0
    ):  # only proceed if there are cells of the specified phenotype level in the image

        # create dummy features so we can use this method of graph construction
        features = np.random.rand(len(image_objects), 3)
        cellTypes = image_objects[phenotyping_column].values

        coords = np.asarray(
            list(zip(image_objects[x_id].values, image_objects[y_id].values))
        )

        # create anndata object:
        adata = AnnData(
            features, obsm={"spatial": coords}, obs={f"{phenotyping_column}": cellTypes}
        )

        ## derive spatial neighbours graph:
        if graph_type == "spatial_neighbours":
            ## derive spatial neighbours graph:
            sq.gr.spatial_neighbors(adata, radius=radius, coord_type="generic")
        elif graph_type == 'nearest_neighbour':
            ## derive nearest neighbours graph:
            sq.gr.spatial_neighbors(adata, n_neighs=n_neighs)
        else:
            raise ValueError("graph_type not recognised.")

        G = nx.convert_matrix.from_scipy_sparse_matrix(
                adata.obsp["spatial_connectivities"]
            )
        node_ids = get_graph_node_ids(adata, phenotyping_column)
        
    return G, node_ids


def compute_shortest_paths(G, node_ids, source=0):
    '''
    Compute shortest with breadth first search from source node to all other nodes in the graph.
    Args:
    G = networkx graph

    '''
    sps = cg.traversal.bfs(G, start=source)
    sps = merge_node_ids_to_shortest_paths(sps, node_ids)
    return sps


def merge_node_ids_to_shortest_paths(shortest_paths, node_ids):
    # node_ids = pandas dataframe with nodeid in integer format as column called 'vertex'
    shortest = pd.merge(shortest_paths, node_ids, how="left", on="vertex")
    return shortest

def min_path_to_cellType(shortest_paths, phenotyping_column, cType):
    min_path = shortest_paths[shortest_paths[phenotyping_column] == cType][
        "distance"
    ].min()
    return int(min_path)


def get_predecessor(df, minpath_to_target, phenotyping_column, target_cellType):
    # print(df)
    predecessor = df[
        (df["distance"] == minpath_to_target)
        & (df[phenotyping_column] == target_cellType)
    ]["predecessor"].values.tolist()
    # print('initial_predecessor: ', predecessor)
    return predecessor


def get_chain_predecessor(df, predecessor):
    chain_predecessor = df[(df["vertex"] == predecessor)].iloc[0]["predecessor"].item()
    return chain_predecessor


def followchain(
    shortest_paths,
    minpath_to_target,
    source_cell,
    source_vertex,
    target_cell="Epithelial cells",
    phenotyping_column="cellType",
):
    """
    The min_path_to_target is e.g. 3 jumps from the source to the target cell;
    there may be more than one path of the same length.

    Todo: return vertex ids for subsequent referencing of e.g. positive marker expression
    """
    # print("SHORTEST PATHS: \n", shortest_paths)
    sorted_shortest_paths = shortest_paths.sort_values(by="distance")
    sorted_shortest_paths = sorted_shortest_paths[
        (sorted_shortest_paths["distance"] == minpath_to_target)
        & (sorted_shortest_paths[phenotyping_column] == target_cell)
    ]
    # print("SORTED SHORTEST PATHS\n", sorted_shortest_paths)
    target_vertex = sorted_shortest_paths["vertex"].iloc[0]

    # print(minpath_to_target)
    if minpath_to_target == 1:
        vertexes = [target_vertex, source_vertex]
        chaincells = [target_cell, source_cell]
        # print(chaincells)
        return chaincells, vertexes

    else:
        ## calculate internal chain:

        if (
            sorted_shortest_paths["predecessor"].max() >= 0
        ):  # possible to have a df with target cells that are disconnected
            chaincells = [
                target_cell
            ]  # store cells in chain between starting cell and destination
            vertexes = [target_vertex]

            predecessors = get_predecessor(
                sorted_shortest_paths,
                minpath_to_target,
                phenotyping_column,
                target_cellType=target_cell,
            )  # get list of initial predecessor cells
            p = predecessors[0]

            for j in range(minpath_to_target - 1):
                # take first predecessor if more than one:

                if p != -1:  ## proceed only if the predecessor is not the source cell
                    predecessor_nodeType = get_predecessor_nodeType(
                        shortest_paths, p, phenotyping_column=phenotyping_column
                    )
                    chaincells.append(predecessor_nodeType)
                    vertexes.append(p)
                    p = get_chain_predecessor(shortest_paths, p)

            chaincells.append(source_cell)
            vertexes.append(source_vertex)
            return chaincells, vertexes


def relabel_vertex_chain(vertex_chain, relabel_dict):
    relabelled_chain = [relabel_dict[item] for item in vertex_chain]
    return relabelled_chain


def get_predecessor_nodeType(df, predecessor, phenotyping_column="cellType"):
    nodetype = df[df["vertex"] == predecessor].iloc[0][phenotyping_column]
    return nodetype


def get_graph_node_ids(spg, phenotyping_column):
    """
    Produce dataframe of graph node ids from squidpy spatial graph.
    """

    node_ids = pd.DataFrame(spg.obs[[phenotyping_column]])
    node_ids["vertex"] = np.arange(0, len(node_ids.index))
    return node_ids


def get_cellType_node_ids(node_ids, phenotyping_column, cellType):
    cellType_node_ids = node_ids[node_ids[phenotyping_column] == cellType]
    return cellType_node_ids





def permute_phenotypes(
    objects, in_place=True, region="all", phenotyping_column="cellType"
):
    """
    Randomly permute phenotype and positivity information for all images in the objects table.

    Args:
    objects = cell objects dataframe with each image as a distinct imagename
    in_place = bool, if True overwrite current cell phenotype data in place,
            otherwise add '_permuted' columns
    region = str, subset permutation to a give tissue region defined by 'region' column.
            'all' = do not subset, permute all cells.

    """
    if region == "all":
        for imagename in objects["imagename"].unique():
            rng = np.random.default_rng(seed=123)
            if in_place:
                objects.loc[
                    objects["imagename"] == imagename,
                    (phenotyping_column, "majorType", "positive"),
                ] = rng.permutation(
                    objects.loc[
                        objects["imagename"] == imagename,
                        (phenotyping_column, "majorType", "positive"),
                    ].values
                )
            else:
                objects.loc[
                    objects["imagename"] == imagename,
                    ("cellType_permuted", "majorType_permuted", "positive_permuted"),
                ] = rng.permutation(
                    objects.loc[
                        objects["imagename"] == imagename,
                        (phenotyping_column, "majorType", "positive"),
                    ].values
                )
    else:
        for imagename in objects["imagename"].unique():
            rng = np.random.default_rng(seed=123)
            if in_place:
                objects.loc[
                    (objects["imagename"] == imagename) & (objects["region"] == region),
                    (phenotyping_column, "majorType", "positive"),
                ] = rng.permutation(
                    objects.loc[
                        (objects["imagename"] == imagename)
                        & (objects["region"] == region),
                        (phenotyping_column, "majorType", "positive"),
                    ].values
                )
            else:
                objects.loc[
                    (objects["imagename"] == imagename) & (objects["region"] == region),
                    (
                        f"cellType_{region}_permuted",
                        f"majorType_{region}_permuted",
                        f"positive_{region}_permuted",
                    ),
                ] = rng.permutation(
                    objects.loc[
                        (objects["imagename"] == imagename)
                        & (objects["region"] == region),
                        (phenotyping_column, "majorType", "positive"),
                    ].values
                )
                objects.loc[
                    (objects["imagename"] == imagename)
                    & ~(objects["region"] == region),
                    (
                        f"cellType_{region}_permuted",
                        f"majorType_{region}_permuted",
                        f"positive_{region}_permuted",
                    ),
                ] = objects.loc[
                    (objects["imagename"] == imagename)
                    & ~(objects["region"] == region),
                    (phenotyping_column, "majorType", "positive"),
                ].values

    return objects

