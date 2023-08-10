#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks 1
#SBATCH --mem 8GB
#SBATCH --time 48:00:0
#SBATCH --job-name p2_pipe
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=alastair.magness@crick.ac.uk

#!/bin/bash

ml purge
ml Nextflow/22.04.0
ml Singularity/3.6.4
ml CUDA/11.4.1

export NXF_SINGULARITY_CACHEDIR='./singularity'

nextflow run ./main.nf \
    --objects "/nemo/project/proj-tracerx-lung/tctProjects/rubicon/PHLEX/revision_testing/angelom_typex_mccs/cell_objects_revision_p1.txt"\
    --objects_delimiter "\t" \
    --image_id_col "imagename"\
    --phenotyping_column 'majorType'\
    --phenotype_to_cluster 'Epithelial cells'\
    --x_id "centerX"\
    --y_id "centerY"\
    --barrier_phenotyping_column "majorType" \
    --outdir "../results" \
    --release 'PHLEX_edit_chatgpt_2' \
    --workflow_name 'clustered_barrier' \
    --barrier_source_cell_type "CD8 T cells"\
    --barrier_target_cell_type "Epithelial cells"\
    --barrier_cell_type "aSMA+ cells"\
    --singularity_bind_path '/camp,/nemo'\
    --n_neighbours 5\
    --dev \
    --number_of_inputs 1 \
    -w "/camp/project/proj-tracerx-lung/txscratch/rubicon/deep_imcyto/work"\
    -resume