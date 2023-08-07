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
    --objects "/camp/lab/swantonc/working/Alastair/other_crick_imaging/megan_cole/data/dataset2_communityC_MRTX_Tregs_barrier_scoring_input_2.csv"\
    --objects_delimiter ','\
    --image_id_col "Image_ID"\
    --x_id "Location_Center_X"\
    --y_id "Location_Center_Y"\
    --barrier_phenotyping_column "Cell_phenotype" \
    --outdir "../results" \
    --release 'MRTX_barrier' \
    --workflow_name 'barrier_only' \
    --barrier_source_cell_type "T cells CD8"\
    --barrier_target_cell_type "Dendritic cells"\
    --barrier_cell_type "T reg cells"\
    --singularity_bind_path '/camp,/nemo'\
    --n_neighbours 2\
    -w '/camp/project/proj-tracerx-lung/txscratch/rubicon/Spatial-PHLEX/work'\
    -resume