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


nextflow run ./main.nf \
    --BARRIER_DELIMITER ',' \
    --CALCULATE_BARRIER true \
    --METADATA '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'\
    --METADATA_DELIMITER '\t'\
    --OBJECTS "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_p1.txt" \
    --OBJECTS_DELIMITER '\t' \
    --OVERWRITE true \
    --PANEL 'p1' \
    --PHENOTYPING_LEVELS 'cellType' \
    --barrier_phenotyping_level 'cellType' \
    --dev \
    --graph_conda_env "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-22.02" \
    --graph_type 'nearest_neighbour' \
    --md_conda 'Anaconda3' \
    --outdir '../../results' \
    --publish_dir_mode 'copy' \
    --release '2022-08-30_DSL2_dev' \
    --spclust_conda_env "/camp/lab/swantonc/working/Alastair/.conda/envs/tf" \
    --workflow_name 'default' \
    # -resume



    #'cellType_majorType'
    #'spatial_clustering', 'default'