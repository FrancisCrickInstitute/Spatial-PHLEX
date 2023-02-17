#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks 1
#SBATCH --mem 8GB
#SBATCH --time 48:00:0
#SBATCH --job-name p2_pipe
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=alastair.magness@crick.ac.uk

ml purge
ml Nextflow/22.04.0
ml Singularity/3.6.4
ml CUDA/11.4.1


export NXF_SINGULARITY_CACHEDIR='/camp/project/proj-tracerx-lung/tctProjects/rubicon/inputs/containers/deep-imcyto'


nextflow run ./main.nf \
    --metadata '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/archive/metadata.tracerx.220815.txt'\
    --metadata_delimiter '\t'\
    --objects '/camp/project/proj-tracerx-lung/tctProjects/rubicon/PHLEX/release_testing/Spatial-PHLEX/data/PHLEX_test_data.csv'\
    --objects_delimiter '\t'\
    --phenotyping_column 'majorType' \
    --barrier_phenotyping_column 'majorType' \
    --graph_type 'nearest_neighbour' \
    --outdir '../results' \
    --release 'PHLEX_test' \
    --workflow_name 'default' \
    --barrier_source_cell_type 'CD8 T cells'\
    --barrier_target_cell_type 'Epithelial cells'\
    --barrier_cell_type 'Myofibroblasts'\
    --dev false \
    -w '/camp/project/proj-tracerx-lung/txscratch/rubicon/Spatial-PHLEX/work'\
    # -resume