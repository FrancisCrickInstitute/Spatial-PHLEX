#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks 1
#SBATCH --mem 8GB
#SBATCH --time 48:00:0
#SBATCH --job-name SpatialPHLEX

ml purge
ml Nextflow/22.04.0
ml Singularity/3.6.4

export NXF_SINGULARITY_CACHEDIR='./singularity'

nextflow run ./main.nf \
    --workflow_name 'clustered_barrier' \
    --objects "$PWD/data/PHLEX_test_cell_objects.txt"\
    --objects_delimiter "\t" \
    --image_id_col "imagename"\
    --phenotyping_column 'majorType'\
    --phenotype_to_cluster 'Epithelial cells'\
    --x_coord_col "centerX"\
    --y_coord_col "centerY"\
    --barrier_phenotyping_column "majorType" \
    --barrier_source_cell_type "CD8 T cells"\
    --barrier_target_cell_type "Epithelial cells"\
    --barrier_cell_type "aSMA+ cells"\
    --n_neighbours 10\
    --outdir "../results" \
    --release 'PHLEX_testing' \
    --singularity_bind_path '/camp,/nemo'\
    --plot_palette "$PWD/assets/PHLEX_test_palette.json" \
    --dev true\
    --number_of_inputs 1\
    -w "/camp/project/proj-tracerx-lung/txscratch/rubicon/deep_imcyto/work"\
    -profile crick \
    -resume