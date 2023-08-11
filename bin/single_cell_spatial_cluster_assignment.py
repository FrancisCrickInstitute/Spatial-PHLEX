#!/usr/bin/env python

import os
import argparse

import numpy as np
import pandas as pd
from tqdm import *
import skimage.io as io
import json

from spclust_util import (
    assemble_cluster_polygons,
    assign_to_spclust,
    do_clustering,
    get_image_shape_from_sampleFile,
    get_shape_from_objects,
    plot_clusters,
    save_cluster_alphashapes,
    pipeline_make_label,
)

from spcluster_properties import intracluster_density

class Config:
    def __init__(self, args):
        self.imagename = args.imagename
        self.phenotyping_column = args.phenotyping_column
        self.phenotype_to_cluster = args.phenotype_to_cluster
        self.eps = args.eps  # eps parameter for dbscan
        self.min_s = args.min_s  # min number of cells in a cluster
        self.objects_filepath = args.objects_filepath
        self.objects_sep = args.objects_sep
        self.sampleFile = args.sampleFile
        self.sampleFile_sep = args.sampleFile_sep
        self.x_coord = args.x_coord
        self.y_coord = args.y_coord
        self.alpha = 0.05  # Alphashape curvature parameter
        self.root_outdir = f"{args.root_outdir}/{args.phenotyping_column}/eps_{args.eps}_min_size_{args.min_s}_alpha_{args.alpha}"
        self.palette = args.plot_palette
        self.image_id_col = args.image_id_col

    def display(self):
        """Display Configuration values."""
        print("\n" + "~" * 140)
        print("CONFIG:")
        for attribute, value in vars(self).items():
            if not attribute.startswith("__") and not callable(value):
                print(f"{attribute:30} {value}")
        print("~" * 140 + "\n")


def process_image(imagename, cType, image_df, imshape):
    """Cluster a given image with DBSCAN for the specified cType.
    Args:
        imagename (str): name of image to be processed
        cType (str): cell type to be clustered
        image_df (pd.DataFrame): dataframe containing cell objects for that image
        imshape (tuple): shape of image

    Returns:
        None

    Creates:
        - output directory for image
        - cluster assignments for each cell
        - cluster polygons
        - cluster alphashapes
        - cluster assignment label image
        - plots of clusters
    """

    out_dir = os.path.join(CONFIG.root_outdir, f"{cType}/{imagename}")
    os.makedirs(out_dir, exist_ok=True)

    clusters = image_df.loc[image_df[CONFIG.phenotyping_column] == cType]

    if len(clusters) > 0:
        cluster_labels = do_clustering(clusters, CONFIG.eps, CONFIG.min_s)

        clusters["dbscan_cluster"] = cluster_labels

        spatial_cluster_polygons = assemble_cluster_polygons(
            clusters, cluster_labels, imagename, cType, x_col_id=CONFIG.x_coord, y_col_id=CONFIG.y_coord, alpha=CONFIG.alpha
        )
        alphashape_df = save_cluster_alphashapes(
            spatial_cluster_polygons, cType, imagename, out_dir, CONFIG.phenotyping_column
        )

        all_cells_list = list(
            zip(image_df[CONFIG.x_coord].values, image_df[CONFIG.y_coord].values)
        )
        (
            spclust_assignments,
            cluster_areas,
            nearest_cluster,
            distances,
        ) = assign_to_spclust(all_cells_list, spatial_cluster_polygons)

        image_df[f"{cType}_spatial_cluster_id"] = spclust_assignments
        image_df[f"{cType}_cluster_area"] = cluster_areas
        image_df[f"nearest_{cType}_cluster_id"] = nearest_cluster
        image_df[f"distance_to_nearest_{cType}_cluster_boundary"] = distances

        alphashape_spath = os.path.join(
            out_dir, f"{imagename}_{cType}_alphashape_polygons_label.png"
        )


    else:
        clusters["dbscan_cluster"] = []
        cols_to_fill = [
            f"{cType}_{col}"
            for col in [
                "spatial_cluster_id",
                "cluster_area",
                "nearest_cluster_id",
                "distance_to_nearest_cluster_boundary",
            ]
        ]
        image_df[cols_to_fill] = -1, np.nan, -1, np.nan

        alphashape_mask = np.zeros(imshape)
        alphashape_spath = os.path.join(
            out_dir, f"{imagename}_{cType}_alphashape_polygons_label.tiff"
        )
        io.imsave(alphashape_spath, alphashape_mask)

    spath = os.path.join(out_dir, f"{imagename}_object_cluster_assignment.csv")
    image_df.to_csv(spath, sep=CONFIG.objects_sep)

    # load palette:

    if CONFIG.palette is not None:
        # load the palette from the config filepath:
        with open(CONFIG.palette, "r") as f:
            plot_palette = json.load(f)

        plot_clusters(image_df, 
                    cluster_id_col = f'{cType}_spatial_cluster_id', 
                    clustering_cell_type = cType, 
                    x_col_id = CONFIG.x_coord,
                    y_col_id = CONFIG.y_coord,
                    phenotyping_column = CONFIG.phenotyping_column,
                    bg_image = None, 
                    image_shape = imshape, 
                    sample_name = imagename,  
                    outdir=out_dir, 
                    alphashape_param = CONFIG.alpha,
                    palette = plot_palette[CONFIG.phenotyping_column])
    else:
        plot_clusters(image_df, 
                    cluster_id_col = f'{cType}_spatial_cluster_id', 
                    clustering_cell_type = cType, 
                    x_col_id = CONFIG.x_coord,
                    y_col_id = CONFIG.y_coord,
                    phenotyping_column = CONFIG.phenotyping_column,
                    bg_image = None, 
                    image_shape = imshape, 
                    sample_name = imagename,  
                    outdir=out_dir, 
                    alphashape_param = CONFIG.alpha,
                    palette = None)
        
    intracluster_density(image_df, CONFIG.objects_sep, CONFIG.phenotyping_column,CONFIG.imagename, out_dir, cType)

    print(f"\n{imagename} done.")


########
# MAIN #
########


def main(CONFIG):
    # display configuration:
    CONFIG.display()

    # define imagename to be processed:
    imagename = CONFIG.imagename

    # read df:
    cell_objects = pd.read_csv(
        CONFIG.objects_filepath, sep=CONFIG.objects_sep, encoding="latin1"
    )

    # read sampleFile
    if CONFIG.sampleFile != None:
        sampleFile = pd.read_csv(
            CONFIG.sampleFile, sep=CONFIG.sampleFile_sep, encoding="latin1"
        )
        imshape = get_image_shape_from_sampleFile(sampleFile, imagename, image_id_col=CONFIG.image_id_col)
    else:
        imshape = get_shape_from_objects(cell_objects, x_coord=CONFIG.x_coord, y_coord=CONFIG.y_coord)

    if CONFIG.phenotype_to_cluster == "all":
        clustering_cell_types = cell_objects[
            CONFIG.phenotyping_column
        ].unique()  # ['Epithelial cells']
    else:
        clustering_cell_types = [CONFIG.phenotype_to_cluster]

    # create base output directory:
    if os.path.exists(CONFIG.root_outdir) != True:
        os.makedirs(CONFIG.root_outdir)

    image_df = cell_objects.loc[cell_objects[CONFIG.image_id_col] == imagename]

    for cType in clustering_cell_types:
        process_image(imagename, cType, image_df, imshape)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--imagename",
        type=str,
        default="test_image",
        help="name of image to be processed",
    )
    parser.add_argument(
        "--objects_filepath",
        type=str,
        default="test_objects.csv",
        help="path to objects file",
    )
    parser.add_argument(
        "--objects_sep", type=str, default=";", help="separator in objects file"
    )
    parser.add_argument(
        "--sampleFile", type=str, default=None, help="path to sampleFile file"
    )
    parser.add_argument(
        "--sampleFile_sep", type=str, default="\t", help="separator in sampleFile file"
    )
    parser.add_argument(
        "--root_outdir", type=str, default=".", help="path to output directory"
    )
    parser.add_argument(
        "--eps", type=float, default=25, help="eps parameter for dbscan"
    )
    parser.add_argument(
        "--min_s", type=int, default=1, help="min number of cells in a cluster"
    )
    parser.add_argument(
        "--alpha", type=float, default=0.05, help="alpha parameter for alphashape"
    )
    parser.add_argument(
        "--phenotyping_column",
        type=str,
        default="cellType",
        help="what level of phenotyping, majorType or cellType",
    )
    parser.add_argument(
        "--phenotype_to_cluster",
        type=str,
        default="all",
        help="which cell types to cluster, all or specific cell type",
    )
    parser.add_argument("--x_coord", type=str, default='centerX', help="x coordinate")
    parser.add_argument("--y_coord", type=str, default='centerY', help="y coordinate")
    parser.add_argument("--image_id_col", type=str, default='imagename', help="image id column")
    parser.add_argument("--plot_palette", type=str, default=None, help="path to json file containing palette for plotting")
    args = parser.parse_args()

    # create configuration based on input args:
    CONFIG = Config(args)

    # pass CONFIG to main():
    main(CONFIG)
