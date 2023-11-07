from collections import Counter
import sb_functions as sb
import pandas as pd
from tqdm import tqdm

def measure_barrier(graph, 
                    node_ids, 
                    phenotyping_column, 
                    imagename, 
                    source_cell_type = 'CD8 T cells', 
                    target_cell_type = 'Epithelial cells_tumour', 
                    barrier_types = ['Myofibroblasts'],
                    image_id_col = 'imagename',
                    max_path_length = 1000):
    
    """Measure the barrier presented by a cell type between two other cell types in a spatial graph.

    Args:
        graph (cugraph.Graph): A graph object from cuGraph.
        node_ids (pandas.DataFrame): A dataframe with the node ids of the graph.
        phenotyping_column (str): The column name of the cell type in the node_ids dataframe.
        imagename (str): The name/id of the image being processed.
        source_cell_type (str, optional): The cell type to start the shortest path from. Defaults to 'CD8 T cells'.
        target_cell_type (str, optional): The cell type to end the shortest path at. Defaults to 'Epithelial cells_tumour'.
        barrier_types (list, optional): The cell types that are considered barriers. Defaults to ['Myofibroblasts'].
        max_path_length (int, optional): The maximum path length to consider. Defaults to 1000.

    Returns:
        barrier_df (pd.DataFrame): Pandas dataframe with barrier measurements for every cell of the source cell type present in the spatial graph.
    """    

    ## initialise lists to store results:
    minimum_paths = []
    all_vertexes = []
    all_chain_lengths = []
    all_allpaths_counts = []
    all_allpaths_adjacent = []

    ## node ids of the target cell type:
    cellType_node_ids = sb.get_cellType_node_ids(
        node_ids, phenotyping_column, source_cell_type
    )

    assert(len(cellType_node_ids.index) > 0), "No cells of the source cell type in the graph."

    ## loop through all cells of the target type:
    for v in tqdm(cellType_node_ids.vertex, ascii=True):

        if (
            v >= 0
        ):  # Starting vertex should be between 0 to number of vertices
            ## compute shortest paths of vertex to all other nodes in the graph:
            shortest_paths = sb.compute_shortest_paths(
                graph, node_ids, source=v
            )

            ## what is the minimum path length to an epithelial cell?
            minpath_to_target = sb.min_path_to_cellType(
                shortest_paths,
                phenotyping_column,
                target_cell_type,
            )

            # restrict to paths that are shorter than the max path length
            if minpath_to_target < max_path_length:
                minimum_paths.append(minpath_to_target)

                closest_target = shortest_paths[
                    (shortest_paths["distance"] == minpath_to_target)
                    & (
                        shortest_paths[phenotyping_column]
                        == f"{target_cell_type}"
                    )
                ]

                allpaths_barrier_fraction = (
                    allpaths_path_content(
                        closest_target,
                        shortest_paths,
                        minpath_to_target,
                        barrier_cells=barrier_types,
                        phenotyping_column=phenotyping_column,
                    )
                )
                allpaths_adjacent_count = (
                    allpaths_adjacent_barrier(
                        closest_target,
                        shortest_paths,
                        minpath_to_target,
                        barrier_cells=barrier_types,
                        phenotyping_column=phenotyping_column,
                    )
                )

                all_allpaths_counts.append(
                    allpaths_barrier_fraction
                )
                all_allpaths_adjacent.append(
                    allpaths_adjacent_count
                )

                all_chain_lengths.append(minpath_to_target)
                all_vertexes.append(int(v))

    ## construct dataframe of results:
    barrier_df = pd.DataFrame(
        data=list(
            zip(
                all_chain_lengths,
                all_allpaths_counts,
                all_allpaths_adjacent, 
                all_vertexes
            )
        ),
        columns=[
            "chainlength",
            "all_paths_barrier_fraction",
            "all_paths_adjacent_fraction", 
            "vertex"
        ],
    )

    barrier_df[image_id_col] = imagename
    barrier_df["internal_chainlength"] = barrier_df["chainlength"] - 2
    barrier_df["source_cell"] = source_cell_type
    barrier_df["target_cell"] = f"{target_cell_type}"
    for i, t in enumerate(barrier_types):
        barrier_df[f'barrier_cell_type_{i+1}'] = t

    return barrier_df

def allpaths_path_content(
    closest_epi,
    shortest_paths,
    minpath_to_epi,
    barrier_cells=["Myofibroblasts"],
    phenotyping_column="cellType",
):
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

    if minpath_to_epi > 1:  #
        path_cells = []
        for i in range(len(closest_epi)):
            p = closest_epi.iloc[i, :]["predecessor"]
            pdf = shortest_paths[shortest_paths["vertex"] == p]
            chain = []

            for j in range(minpath_to_epi - 1):
                if pdf["predecessor"].item() > -1:
                    pdf = shortest_paths[shortest_paths["vertex"] == p]
                    chain.append(pdf)
                    p = pdf["predecessor"].item()
                    path_cells.append(pdf)

            allpaths_path_content = pd.concat(path_cells)

            if (
                set(barrier_cells).intersection(
                    set(allpaths_path_content[phenotyping_column])
                )
                != set()
            ):
                myofibroblast_fraction_degen = (
                    allpaths_path_content[phenotyping_column]
                    .value_counts(normalize=True)[barrier_cells]
                    .sum()
                )
            else:
                myofibroblast_fraction_degen = 0
    else:
        myofibroblast_fraction_degen = 0
    return myofibroblast_fraction_degen


def allpaths_adjacent_barrier(
    closest_epi,
    shortest_paths,
    minpath_to_epi,
    barrier_cells=["Myofibroblasts"],
    phenotyping_column="cellType",
):
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

    if minpath_to_epi > 1:  #
        path_cells = []
        for i in range(len(closest_epi)):
            p = closest_epi.iloc[i, :]["predecessor"]
            pdf = shortest_paths[shortest_paths["vertex"] == p]
            chain = []

            for j in range(minpath_to_epi - 1):
                if pdf["predecessor"].item() > -1:
                    pdf = shortest_paths[shortest_paths["vertex"] == p]
                    chain.append(pdf)
                    p = pdf["predecessor"].item()
                    path_cells.append(pdf)

            allpaths_path_content = pd.concat(path_cells)
            adjacent_cells_df = allpaths_path_content[
                allpaths_path_content["distance"] == (minpath_to_epi - 1)
            ]

            if (
                set(barrier_cells).intersection(
                    set(adjacent_cells_df[phenotyping_column])
                )
                != set()
            ):
                mean_adjacent_barrier = (
                    adjacent_cells_df[phenotyping_column]
                    .value_counts(normalize=True)[barrier_cells]
                    .sum()
                )  # / len(adjacent_cells_df)
            else:
                mean_adjacent_barrier = 0

    else:
        mean_adjacent_barrier = 0
    return mean_adjacent_barrier

def count_barrier_cells(chain, BARRIER_TYPES):
    if len(set(BARRIER_TYPES).intersection(set(chain))) > 0:
        counts = Counter(chain)
        x = sum([counts[cell] for cell in BARRIER_TYPES])
    else:
        x = 0
    return x


def weighted_barrier_count(chain, BARRIER_TYPES):
    """Barrier weighting is inversely proportional to its adjacency to the target cell."""

    count = 0
    for i, elem in enumerate(chain[1:-1]):
        weight = 1 / (i + 1)
        if elem in BARRIER_TYPES:
            count += weight
    return count


def binary_barrier_call(chain, BARRIER_TYPES):
    """Returns 1 if barrier type is in cell chain"""
    if len(set(BARRIER_TYPES).intersection(set(chain))) > 0:
        barrier = 1
    else:
        barrier = 0
    return barrier


def adjacent_barrier_call(chain, BARRIER_TYPES):
    """Returns 1 if barrier type is adjacent to target cell (first cell in chain)"""
    inner_chain = chain[1:-1]

    if len(inner_chain) > 0:
        if inner_chain[0] in BARRIER_TYPES:
            barrier = 1
        else:
            barrier = 0
    else:
        barrier = 0
    return barrier