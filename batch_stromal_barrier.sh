#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --ntasks 1
#SBATCH --gres=gpu:1
#SBATCH --time 3:00:0
#SBATCH --job-name pipe
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=alastair.magness@crick.ac.uk



ml Anaconda3
source activate /camp/home/magnesa/.conda/envs/rapids-0.18

python ./bin/stromal_barrier.py /camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/spatial/src/rubicon-sppipeline/work/38/1063bd953066a2268012c711d54856/cp_module_865/P1_TMA_REC_20190508-roi_16.txt ./ /camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_p1.txt p1