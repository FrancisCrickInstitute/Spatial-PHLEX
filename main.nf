#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*********************************************
* RUBICON NEXTFLOW SPATIAL ANALYSIS PIPELINE *
*********************************************/

include {NEIGHBOURHOOD_WF; NEAREST_NEIGHBOUR_WF} from './workflows/barrier.nf'
include {SPATIAL_CLUSTERING_WF; CLUSTERED_BARRIER_WF } from './workflows/spatial.nf'
include { print_logo; check_params } from './modules/util.nf'

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

ch_phenotyping = ch_phenotyping.first() // this will take only the first phenotyping level even if multiple are specified -- > deprecate multiple input?

workflow {
    
    print_logo()
    check_params()
    
    if (params.workflow_name == 'stromal_barrier_only') {

        if (params.graph_type == 'neighbouRhood') {
            NEIGHBOURHOOD_WF ( ch_nhood, params.neighbourhood_module_no, params.OBJECTS, params.OBJECTS_DELIMITER)
        }
        if (params.graph_type == 'nearest_neighbour') {
            NEAREST_NEIGHBOUR_WF ( ch_nhood, params.neighbourhood_module_no, params.OBJECTS, params.OBJECTS_DELIMITER)
        }
    }

    if (params.workflow_name == 'spatial_clustering') {
        SPATIAL_CLUSTERING_WF ( params.OBJECTS, ch_phenotyping)
    }
    
    if (params.workflow_name == 'default') {
        CLUSTERED_BARRIER_WF ( params.OBJECTS, ch_phenotyping)
        CLUSTERED_BARRIER_WF.GRAPH_BARRIER.out.ch_barrier_results.collectFile('barrier_results.csv', keepHeader: true, skip: 1)
    }

}