#!/usr/bin/env nextflow

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

// directly create a list channel for phenotyping levels for combinatoric input of phenotype levels and imagenames
pheno_list = params.PHENOTYPING_LEVELS?.tokenize(',')
// ch_phenotyping = Channel.fromList(['cellType', 'majorType'])
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

/*****************
* BEGIN PIPELINE *
*****************/

process GENERATE_IMAGENAMES {

    /*
    Generate unique imagenames from the cell objects file. 
    */

    module params.md_conda
    conda params.spclust_conda_env

    output:
    stdout into ch_imagenames
 
    """
    #!/usr/bin/env python
    import pandas as pd
 
    objects = pd.read_csv('${params.OBJECTS}', sep='${params.OBJECTS_DELIMITER}', encoding='latin1')
    imagenames = objects['imagename'].unique().tolist()
    for imagename in imagenames:
        print(imagename)
    """

}

process NEIGHBOURHOOD_GRAPH {
    /*
    * process to create the spatial graph from the cellprofiler neighbourhood output
    */

    // executor "slurm"
	// time "0.25h"
	// clusterOptions "--part=cpu --cpus-per-task=4 --mem=2GB"

    when:
    params.CALCULATE_BARRIER

    module params.md_conda
    conda params.spclust_conda_env

    publishDir "${params.outdir}/${params.RELEASE}/graph/adjacency_lists/neighbourhood", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    val nhood_file from ch_nhood

    output:
    path "*/*.txt" into adj_output_ch

    """
    cp2nx.py $nhood_file ${params.neighbourhood_module_no} ./
    """
}

process GRAPH_BARRIER {
    /*
    Run graph barrier scoring for all cell types.
    */

    when:
    params.CALCULATE_BARRIER

    executor "slurm"
	time "6h"
	clusterOptions "--part=gpu --gres=gpu:1"

    module params.md_conda
    conda params.graph_conda_env

    echo true

    publishDir "${params.outdir}/${params.RELEASE}/graph/barrier", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    val adj_list from adj_output_ch

    output:
    file "**/*.csv" optional true into ch_barrier_results

    """
    stromal_barrier_batch_neighbouRhood.py $adj_list ./ ${params.OBJECTS} ${params.PANEL}
    """

}

// need to concatenate outputs-- use collectFile?
// process CONCAT_BARIER {

//     input:
//     pat
// }

process SPATIAL_CLUSTERING {
    /*
    Perform spatial clustering of cell positions.
    */

    executor "slurm"
	time "0.5h"
	clusterOptions "--part=cpu --cpus-per-task=4 --mem=2GB"

    module params.md_conda
    conda params.spclust_conda_env

    echo true

    publishDir "${params.outdir}/${params.RELEASE}/spatial_clustering", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    tuple imagename, level from ch_imagenames.splitText().map{x -> x.trim()}.combine(ch_phenotyping) //split imagenames and remove trailing newline; create a tuple with channel of phenotyping levels

    output:
    file "**/*cluster_assignment.csv" optional true into ch_spclusters
    file "**/*wkt.csv" optional true into ch_wkts
    file "**/*.png" optional true into ch_cluster_plots

    """
    single_cell_spatial_cluster_assignment.py $imagename ${params.COHORT} ${params.PANEL} $level ${params.OBJECTS} ${params.OBJECTS_DELIMITER}
    """

}

// process DICE {

// }

// process NUCLEAR_MORPH {

// }

// process TOPOLOGICAL_DATA_ANALYSIS {

// }

// process TUMOUR_NORMAL {

// }