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
ml Nextflow/21.04.3

nextflow run ./main.nf \
    --PANEL 'p2' \
    --CALCULATE_BARRIER false \
    --OBJECTS '/camp/lab/swantonc/working/leec/nolan/cn_output/31012022/cells_p2_10n.csv' \
    --PHENOTYPING_LEVELS 'neighborhood10' \
    --OBJECTS_DELIMITER ',' \
