#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*********************************************
* RUBICON NEXTFLOW SPATIAL ANALYSIS PIPELINE *
*********************************************/

/*
COHORT PARAMETERS
*/
params.COHORT = 'tx100'
params.PANEL = 'p1'

/*
 * GRAPH ANALYSIS CONFIG PARAMETERS
 */

params.neighborhood_input = "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/nextflow/${params.PANEL}/publication/*/results/segmentation/*/*/neighbourhood.csv"
params.neighbourhood_module_no = 865
params.md_cuda = "CUDA/10.1.105"
params.md_conda = "Anaconda3" 
params.graph_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-0.18"
params.RELEASE = '2022_02_09_release'
params.CALCULATE_BARRIER = true

/*
 * SPATIAL CLUSTERING CONFIG PARAMETERS
 */

params.OBJECTS = "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_${params.PANEL}.txt"
params.OBJECTS_DELIMITER = '\t'
params.spclust_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/tf"
params.PHENOTYPING_LEVELS = 'cellType'

/*
 * Pipeline execution parameters:
 */

params.dev = false
params.number_of_inputs = 2
params.publish_dir_mode = 'copy'
params.OVERWRITE = true
params.outdir = '../../results'
project_dir = projectDir

/*****************
* BEGIN PIPELINE *
*****************/
workflow {
  
  include { Primary } from './workflows/primary_pipeline.nf'
  Primary ()
  
}