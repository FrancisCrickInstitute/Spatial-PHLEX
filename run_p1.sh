#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks 1
#SBATCH --mem 8GB
#SBATCH --time 48:00:0
#SBATCH --job-name p1_pipe
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
    --metadata '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/archive/metadata.tracerx.220815.txt'\
    --metadata_delimiter '\t'\
    --objects '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/spatial/modified_objects_tables/20221205_release/tx100_cell_objects_tx100_publication_p1_reassigned.csv'\
    --objects_delimiter '\t'\
    --phenotyping_level 'majorType' \
    --barrier_phenotyping_level 'majorType' \
    --graph_type 'nearest_neighbour' \
    --outdir '../../results' \
    --release '2022-12-12_normal_epithelial_barrier/p1' \
    --workflow_name 'default' \
    --barrier_source_cell_type 'CD8 T cells'\
    --barrier_target_cell_type 'Epithelial cells_normal'\
    --barrier_cell_type 'Myofibroblasts'\
    --dev false \
    # -resume