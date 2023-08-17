#!/bin/bash

ml purge
ml Nextflow/22.04.0
ml Singularity/3.6.4

export NXF_SINGULARITY_CACHEDIR='./singularity'

nextflow run ./main.nf \
    --objects "./data/mydata.csv"\
    --objects_delimiter ','\
    --image_id_col "Image_ID"\
    --x_id "Location_Center_X"\
    --y_id "Location_Center_Y"\
    --barrier_phenotyping_column "Phenotype" \
    --outdir "../results" \
    --release 'PHLEX_example' \
    --workflow_name 'default' \
    --barrier_source_cell_type "CD8 T cells"\
    --barrier_target_cell_type "Epithelial cells"\
    --barrier_cell_type "aSMA+ Fibroblasts"\
    --singularity_bind_path '/camp,/nemo'\
    --n_neighbours 5\
    -w './scratch'\
    -profile {your_nf-core_profile}\
    -resume