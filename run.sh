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


nextflow run ./test.nf \
    --workflow_name 'default' \ #'spatial_clustering', 'default'
    --spclust_conda_env "/camp/lab/swantonc/working/Alastair/.conda/envs/tf" \
    --OBJECTS "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_p1.txt" \
    --PANEL 'p1' \
    --CALCULATE_BARRIER true \
    --PHENOTYPING_LEVELS 'cellType' \
    --OBJECTS_DELIMITER '\t' \
    --BARRIER_DELIMITER ',' \
    --graph_type 'nearest_neighbour' \
    --dev \
    --graph_conda_env "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-22.02" \
    --md_conda 'Anaconda3' \
    --outdir '../../results' \
    --release '2022-07-29_DSL2_dev_b' \
    --publish_dir_mode 'copy' \
    --OVERWRITE true \
    --METADATA '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'\
    --METADATA_DELIMITER '\t'\
    # -resume



    #'cellType_majorType'