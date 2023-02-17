#!/usr/bin/env python

import numpy as np
import pandas as pd
import glob, os, sys
import re
from matplotlib.path import Path
import multiprocessing as mp

import skimage.io as io
from skimage.util import img_as_uint
## WKT functions

def reconstruct_labels_from_wkt(shapes, imshape):
    """
    Reconstruct a label image from a dataframe containing wkt-encoded polygons.
    """
    # make a canvas with coordinates:
    x, y = np.meshgrid(np.arange(imshape[0]), np.arange(imshape[1]), indexing='ij')
    x, y = x.flatten(), y.flatten()
    points = np.vstack((x,y)).T 
    image = np.zeros(imshape)

    # loop through wkt=encoded shapes in dataframe:
    for i in range(len(shapes.index)):
        
        # extract given wkt:
        wkt = shapes.loc[i, 'well-known_txt']

        # transform wkt string to coords:
        coords = parse_wkt(wkt)

        # check if coords is a list of lists (i.e. a multipolygon):
        if any(isinstance(el, list) for el in coords) == True:

            # loop through polygons if parse_wkt returns list of polygons
            for j in range(len(coords)):
                p = Path(coords[j]) # make a polygon
                grid = p.contains_points(points)
                if j == 0 :
                    mask = grid.reshape(imshape[0],imshape[1]) 
                else:
                    mask = np.add(mask, grid.reshape(imshape[0],imshape[1]))
                
                image = np.where(mask!=0, i+1, image)

        else:
            # just process single polygon:
            p = Path(coords)
            grid = p.contains_points(points)
            mask = grid.reshape(imshape[0],imshape[1]) # now you have a mask with points inside a polygon
            image = np.where(mask!=0, i+1, image)

    
    return image


def parse_wkt(wkt):
    """
    Parse a well known-text format polygon string into a list of coordinates (for polygons), 
    or a list of list of coords for multipolygons.
    """
    # PROCESS POLYGON STRINGS:
    if wkt.split('(')[0] == 'POLYGON ':
        parsed_wkt = wkt.split('(')[2].split(')')[0].split(', ')
        xs = []
        ys = []
        for i in range(len(parsed_wkt)):
            xs.append(float(parsed_wkt[i].split(' ')[0]))
            ys.append(float(parsed_wkt[i].split(' ')[1]))
        coords = list(zip(xs,ys))

    # PROCESS MULTIPOLYGON:
    elif wkt.split('(')[0] == 'MULTIPOLYGON ':
        coords = []
        # print(wkt)
        parsed_wkt = wkt.split(')), ((') 
        for i in range(len(parsed_wkt)):
            parser = parsed_wkt[i] 
            parser = re.sub(r'\(','',parser)
            parser = re.sub(r'\)','',parser)
            parser = re.sub(r'MULTIPOLYGON ','',parser)
            parser = parser.split(', ')
            xs = []
            ys = []
            for j in range(len(parser)):
                xs.append(float(parser[j].split(' ')[0]))
                ys.append(float(parser[j].split(' ')[1]))
            poly_coords = list(zip(xs,ys))
            coords.append(poly_coords)
    else:
        raise ValueError("wkt polygon type not recognised.")

    return coords

def parse_imagename(path):
    return path.split('/')[-2]

def parse_shape(sampleFile, imagename):
    width = sampleFile[sampleFile['imagename'] == imagename]['image_width'].item()
    height = sampleFile[sampleFile['imagename'] == imagename]['image_height'].item()
    shape = (int(height), int(width))
    return shape

def make_label(wktpath, sampleFile, shape):
    wkt = pd.read_csv(wktpath)

    # get imagename and shape:
    imagename = parse_imagename(wktpath)
    print(f'Processing {imagename}.')
    shape = parse_shape(sampleFile, imagename)

    reconstruction_path = wktpath.replace('alphashape_polygons_wkt.csv', 'alphashape_polygons_label.tiff')

    label_image = reconstruct_labels_from_wkt(wkt, shape)
    label_image = label_image / np.amax(label_image)
    label_image = img_as_uint(label_image)
    io.imsave(reconstruction_path, label_image)

def pipeline_make_label(wkts, shape, path):
    '''
    Reconstruct the label image of the dataframe of alpha shapes in wkt format, wkts.
    ''' 
    label_image = reconstruct_labels_from_wkt(wkts, shape)
    label_image = label_image / np.amax(label_image)
    label_image = img_as_uint(label_image)
    io.imsave(path, label_image)


def main(args):

    wkt_file = args[0]
    print(wkt_file)
    # cellType = args[0]
    # panel = args[1]
    # cellType_wkt_list = glob.glob(f'/camp/lab/swantonc/working/Alastair/spatial_analysis/spatial_clustering/results/single_cell_assignment/2021-09-25_release/tx100/{panel}/majorType/dbscan_25/min_size_0/alpha_0.05/{cellType}_clustering/*/*_alphashape_polygons_wkt.csv')
    sampleFile = pd.read_csv(args[1], sep=args[2], encoding='latin1')
    # print(cellType_wkt_list)
    # print(len(cellType_wkt_list))
    # cellType_wkt_list.sort()

    #     # apply
    # for wkt in cellType_wkt_list:

    # make_label(wkt_file, sampleFile)

if __name__== '__main__':
    main(sys.argv[1:]) 