#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks 1
#SBATCH --mem 8GB
#SBATCH --time 48:00:0
#SBATCH --job-name simbarrier
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=alastair.magness@crick.ac.uk

ml purge
ml Nextflow/22.04.0
ml Singularity/3.6.4
ml CUDA/11.4.1

export NXF_SINGULARITY_CACHEDIR='/camp/project/proj-tracerx-lung/tctProjects/rubicon/inputs/containers/deep-imcyto'

# P1 
nextflow run ./main.nf \
    --objects '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/spatial/simulations/simulated_barrier_objects_2022-03-27.csv'\
    --objects_delimiter ','\
    --phenotyping_level 'cellType' \
    --barrier_phenotyping_level 'cellType' \
    --graph_type 'nearest_neighbour' \
    --outdir '../../results' \
    --release '2022-11-29_full_circular_tumour_simulations' \
    --workflow_name 'default' \
    --barrier_source_cell_type 'CD8 T cells'\
    --barrier_target_cell_type 'Epithelial cells'\
    --barrier_cell_type 'Myofibroblasts'\
    --dev false \
    -resume