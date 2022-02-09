#!/usr/bin/env python

import pandas as pd
import numpy as np
import scanpy as sc
import squidpy as sq
from anndata import AnnData
from numpy.random import default_rng
from squidpy.gr._utils import (
    _get_valid_values,
    _assert_categorical_obs,
    _assert_non_empty_sequence,
)


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


    
    
    # # classes = df['cellClass'].unique()

    # # df.loc[~df['cellClass'].isin(classes), 'cellClass'] = df.loc[~df['cellClass'].isin(classes), 'cellType']
    # print(df)
    # # df.loc[~df["S"].isin(allowed_vals), "S"] = "None"

    

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