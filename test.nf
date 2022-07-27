#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*********************************************
* RUBICON NEXTFLOW SPATIAL ANALYSIS PIPELINE *
*********************************************/

include {NEIGHBOURHOOD_WF; NEAREST_NEIGHBOUR_WF} from './workflows/barrier.nf'

/*
COHORT PARAMETERS
*/
params.COHORT = 'tx100'
params.PANEL = 'p1'
params.release = "2022-07-01_DSL2_dev" //'2022_02_11_release'

/*
 * SPATIAL CLUSTERING CONFIG PARAMETERS
 */

params.do_spatial_clustering = true
params.OBJECTS = "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_p1.txt" //"/camp/lab/swantonc/working/Alastair/graph_analysis/stromal_barrier/simulations/2022-03-28/from_mask/data/all_P1_TMA_REC_20190508-roi_14_BT_25_MTF_0.01_MAF_0.95_MSF_0.05_simulations.csv" // "/camp/lab/swantonc/working/Alastair/graph_analysis/stromal_barrier/simulations/2022-03-27/data/simulated_barrier_objects_2022-03-27.csv"//
params.OBJECTS_DELIMITER = '\t'
params.spclust_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/tf"
params.PHENOTYPING_LEVELS = 'cellType' //'majorType,cellType'
params.MAKE_SPATIAL_CLUSTER_MASKS = true
params.METADATA = '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'
params.METADATA_DELIMITER = '\t'

/*
 * GRAPH ANALYSIS CONFIG PARAMETERS
 */

params.graph_type = 'neighbouRhood' // 'nearest_neighbour' //
params.neighborhood_input = "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/nextflow/p1/publication/*/results/segmentation/*/*/neighbourhood.csv"
params.neighbourhood_module_no = 865
params.md_cuda = "CUDA/10.1.105"
params.md_conda = "Anaconda3" 
params.graph_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-0.18"
params.CALCULATE_BARRIER = true
params.BARRIER_DELIMITER = ','


/*
 * Pipeline execution parameters:
 */

params.dev = true
params.number_of_inputs = 2
params.publish_dir_mode = 'copy'
params.OVERWRITE = true
params.outdir = '../../results'
project_dir = projectDir

// directly create a list channel for phenotyping levels for combinatoric input of phenotype levels and imagenames
pheno_list = params.PHENOTYPING_LEVELS?.tokenize(',')
ch_phenotyping = Channel.fromList(pheno_list)

// channel for neighbourhood input csv files
if (params.neighborhood_input) {
    Channel
        .fromPath(params.neighborhood_input, checkIfExists: true)
        .take( params.dev ? params.number_of_inputs : -1 )
        .map { it }
        .ifEmpty { exit 1, "Input file not found: ${params.neighborhood_input}" }
        .set { ch_nhood }
} else {
   exit 1, "Neighbourhood input file not specified!"
}


workflow {

    if (params.graph_type == 'neighbouRhood') {
        NEIGHBOURHOOD_WF ( ch_nhood, params.neighbourhood_module_no, params.OBJECTS, params.OBJECTS_DELIMITER)
    }

}