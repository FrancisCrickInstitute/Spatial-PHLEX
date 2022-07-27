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

# nextflow run ./main.nf \
#     --PANEL 'p2' \
#     --CALCULATE_BARRIER false \
#     --OBJECTS '/camp/lab/swantonc/working/leec/nolan/cn_output/31012022/cells_p2_10n.csv' \
#     --PHENOTYPING_LEVELS 'neighborhood10' \
#     --OBJECTS_DELIMITER ',' \

nextflow run ./test.nf \
    --OBJECTS "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_p1.txt" \
    --PANEL 'p1' \
    --CALCULATE_BARRIER true \
    --PHENOTYPING_LEVELS 'cellType_majorType' \
    --OBJECTS_DELIMITER '\t' \
    --BARRIER_DELIMITER "\t" \
    --graph_type 'neighbouRhood' \
    --dev \
    --graph_conda_env "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-22.02" \
    --md_conda 'Anaconda3' \
    --outdir '../../results' \
    --release '2022-07-01_DSL2_dev' \
    --publish_dir_mode 'copy' \
    --OVERWRITE true \
    -resume