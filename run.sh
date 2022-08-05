#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks 1
#SBATCH --mem 8GB
#SBATCH --time 3:00:0
#SBATCH --job-name pipe
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=alastair.magness@crick.ac.uk

ml purge
ml Nextflow/22.04.0

# ml Singularity/2.6.0-foss-2016b


ml Singularity/3.6.4

ml CUDA/11.4.1

export SINGULARITY_CACHEDIR='/camp/project/proj-tracerx-lung/tctProjects/rubicon/inputs/containers/spatial_phlex/test'
export NXF_SINGULARITY_CACHEDIR='/camp/project/proj-tracerx-lung/tctProjects/rubicon/inputs/containers/spatial_phlex/test'


nextflow run ./main.nf \
    --metadata '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'\
    --metadata_delimiter '\t'\
    --objects "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_p1.txt" \
    --objects_delimiter '\t'\
    --PANEL 'p1' \
    --phenotyping_level 'cellType' \
    --barrier_phenotyping_level 'cellType' \
    --graph_conda_env "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-22.02" \
    --graph_type 'nearest_neighbour' \
    --md_conda 'Anaconda3' \
    --outdir '../../results' \
    --release '2022-08-05' \
    --spclust_conda_env "/camp/lab/swantonc/working/Alastair/.conda/envs/tf" \
    --workflow_name 'default' \
    --dev \
    -resume



    #'cellType_majorType'
    #'spatial_clustering', 'default'