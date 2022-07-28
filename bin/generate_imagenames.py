#!/usr/bin/env python

import pandas as pd
import argparse

def main(args):
        
    objects = pd.read_csv(args.objects, sep=args.delimiter, encoding=args.encoding)
    imagenames = objects['imagename'].unique().tolist()
    for imagename in imagenames:
        print(imagename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--objects', help='path to objects.csv')
    parser.add_argument('--delimiter', default='\t', help='delimiter for objects.csv')
    parser.add_argument('--encoding', default='latin1', help='encoding for objects.csv')
    args = parser.parse_args()
    main(args)