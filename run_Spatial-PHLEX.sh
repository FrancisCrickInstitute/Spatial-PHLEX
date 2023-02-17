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
    --sampleFile "$PWD/data/sample_data.tracerx.txt"\
    --objects "$PWD/../data/PHLEX_test_data.csv"\
    --phenotyping_column 'majorType' \
    --barrier_phenotyping_column 'majorType' \
    --outdir '../results_2023-02-17' \
    --release 'PHLEX_test' \
    --workflow_name 'clustered_barrier' \
    --barrier_source_cell_type 'CD8 T cells'\
    --barrier_target_cell_type 'Epithelial cells'\
    --barrier_cell_type 'aSMA+ cells'\
    -w '/camp/project/proj-tracerx-lung/txscratch/rubicon/Spatial-PHLEX/work'\
    # -resume