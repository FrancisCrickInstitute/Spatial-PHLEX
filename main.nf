#!/usr/bin/env nextflow

/*
 * GRAPH ANALYSIS CONFIG PARAMETERS
 */

 //
params.neighborhood_input = '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/nextflow/p1/publication/*/results/segmentation/*/*/neighbourhood.csv'
params.neighbourhood_module_no = 865
params.md_cuda = "CUDA/10.1.105"
params.md_conda = "Anaconda3" 
params.graph_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-0.18"
params.publish_dir_mode = 'copy'
params.outdir = '../results'
project_dir = projectDir
params.imagenames = './inventory/p2_tumour_communities_imagenames.csv' //redirect to inventory

/*
 * SPATIAL CLUSTERING CONFIG PARAMETERS
 */
// params.METADATA = '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'

params.spclust_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/tf"


// // # EPS density parameter:
// params.EPS = 25

// // # minimum samples for clustering
// params.MIN_S = 0

// // # alphashape curvature parameter:
// params.ALPHA = 0.05  

// params.RELEASE_VERSION = '2021-09-25_release'

// // # base output directory:
// params.ROOT_OUT_DIR = '../../results/single_cell_assignment/{}/{}/{}/{}/dbscan_{}/min_size_{}/alpha_{}'.format(params.RELEASE_VERSION, params.COHORT, params.PANEL, params.PHENOTYPING_LEVEL, params.EPS, params.MIN_S, params.ALPHA)
params.COHORT = 'tx100'
params.PANEL = 'p1'
params.PHENOTYPING_LEVEL = 'neighborhood10' //'cellType'
params.OBJECTS = '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_p1.txt'

params.dev = false
params.number_of_inputs = 2

Channel
    .fromPath(params.imagenames)
    .splitCsv(header:true)
    .map{ row->row.imagename }
    .set { ch_imagenames }

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

ch_panels = Channel.value('p2')
ch_phenotyping = Channel.fromList(['cellType', 'majorType'])

// process GENERATE_IMAGENAMES {

//     module params.md_conda
//     conda params.spclust_conda_env

//     output:
//     stdout into ch_imagenames_2
 
//     """
//     #!/usr/bin/env python
//     import pandas as pd
 
//     objects = pd.read_csv('${params.OBJECTS}', sep='\t', encoding='latin1')
//     imagenames = objects['imagename'].unique().tolist()
//     for imagename in imagenames:
//         print(imagename)
//     """

// }

process NEIGHBOURHOOD_GRAPH {
    /*
    * process to create the spatial graph from the cellprofiler neighbourhood output
    */

    // executor "slurm"
	// time "0.25h"
	// clusterOptions "--part=cpu --cpus-per-task=4 --mem=2GB"

    module params.md_conda
    conda params.spclust_conda_env

    publishDir "${params.outdir}/graph/adjacency_lists/neighbourhood", mode: params.publish_dir_mode

    input:
    val nhood_file from ch_nhood

    output:
    path "*/*.txt" into adj_output_ch

    """
    cp2nx.py $nhood_file ${params.neighbourhood_module_no} ./
    """
}

process GRAPH_BARRIER {

    executor "slurm"
	time "6h"
	clusterOptions "--part=gpu --gres=gpu:1"

    module params.md_conda
    conda params.graph_conda_env

    echo true

    publishDir "${params.outdir}/graph/barrier", mode: params.publish_dir_mode

    input:
    val adj_list from adj_output_ch

    output:
    file "**/*.csv" optional true into ch_barrier_results

    """
    stromal_barrier_batch_neighbouRhood.py $adj_list ./ ${params.OBJECTS} ${params.PANEL}
    """

}

// process SPATIAL_CLUSTERING {

//     executor "slurm"
// 	time "0.5h"
// 	clusterOptions "--part=cpu --cpus-per-task=4 --mem=2GB"

//     module params.md_conda
//     conda params.spclust_conda_env

//     echo true

//     publishDir "${params.outdir}", mode: params.publish_dir_mode

//     input:
//     val imagename from ch_imagenames

//     """
//     single_cell_spatial_cluster_assignment.py $imagename ${params.COHORT} ${params.PANEL} ${params.PHENOTYPING_LEVEL} ${params.OBJECTS}
//     """

// }

// process DICE {

// }

// process NUCLEAR_MORPH {

// }

// process TOPOLOGICAL_DATA_ANALYSIS {

// }