#!/usr/bin/env bash


python ./argparse_test.py --graph_type nearest_neighbour \
        --imagename test354 \
        --neighbours 10 \
        --root_out ./ \
        --objects_path /path.to/objects \
        --objects_sep '\t' \
        --panel p1 \
        --calc_chain True \
        --barrier_types Myofibroblasts \
        --phenotyping_level cellType \

# python ./argparse_test.py --graph_type nearest_neighbour --imagename test354 --neighbours 10 --root_out ./ --objects_path /path.to/objects --objects_sep '\t' --panel p1 --calc_chain True --barrier_types Myofibroblasts --phenotyping_level cellType