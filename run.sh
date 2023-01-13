#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks 1
#SBATCH --mem 8GB
#SBATCH --time 4:00:0
#SBATCH --job-name pipe
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=alastair.magness@crick.ac.uk

ml purge
ml Nextflow/22.04.0

# ml Singularity/2.6.0-foss-2016b


ml Singularity/3.6.4

ml CUDA/11.4.1

# export SINGULARITY_CACHEDIR='/camp/project/proj-tracerx-lung/tctProjects/rubicon/inputs/containers/spatial_phlex/test'
# export NXF_SINGULARITY_CACHEDIR='/camp/project/proj-tracerx-lung/tctProjects/rubicon/inputs/containers/spatial_phlex/test'
export SINGULARITY_CACHEDIR='/camp/project/proj-tracerx-lung/tctProjects/rubicon/inputs/containers/deep-imcyto'
export NXF_SINGULARITY_CACHEDIR='/camp/project/proj-tracerx-lung/tctProjects/rubicon/inputs/containers/deep-imcyto'

nextflow run ./main.nf \
    --metadata '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/archive/metadata.tracerx.220815.txt'\
    --metadata_delimiter '\t'\
    --objects '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/spatial/modified_objects_tables/2022-09-09_tx100_cell_objects_tx100_publication_p2_reassigned_majorType_for_spclustering.csv'\
    --objects_delimiter '\t'\
    --phenotyping_column 'majorType' \
    --barrier_phenotyping_column 'majorType' \
    --graph_type 'nearest_neighbour' \
    --outdir '../../results' \
    --release '2022-09-21_combined_tumour_assignments' \
    --workflow_name 'spatial_clustering' \
    --barrier_source_cell_type 'CD8 T cells'\
    --barrier_target_cell_type 'Epithelial cells_tumour'\
    --barrier_cell_type 'Myofibroblasts'\
    --dev false \
    -resume
    # --dev



    #'cellType_majorType'
    #'spatial_clustering', 'default'
