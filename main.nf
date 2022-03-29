
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

params.graph_type = 'nearest_neighbour' //neighbouRhood'
params.neighborhood_input = "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/nextflow/${params.PANEL}/publication/*/results/segmentation/*/*/neighbourhood.csv"
params.neighbourhood_module_no = 865
params.md_cuda = "CUDA/10.1.105"
// params.md_cuda = "CUDA/11.1.1-GCC-10.2.0"
params.md_conda = "Anaconda3" 
params.graph_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-0.18"
params.RELEASE = "development_testing_20220329" //'2022_02_11_release'
params.CALCULATE_BARRIER = true

/*
 * SPATIAL CLUSTERING CONFIG PARAMETERS
 */

params.do_spatial_clustering = false
params.OBJECTS = "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_${params.PANEL}.txt"
params.OBJECTS_DELIMITER = '\t'
params.spclust_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/tf"
params.PHENOTYPING_LEVELS = 'majorType,cellType'
params.MAKE_SPATIAL_CLUSTER_MASKS = true
params.METADATA = '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'
params.METADATA_DELIMITER = '\t'

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

    module params.md_conda
    conda params.graph_conda_env

    publishDir "${params.outdir}/${params.RELEASE}/graph/adjacency_lists/neighbourhood", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    val nhood_file from ch_nhood

    output:
    path "*/*.txt" into adj_output_ch

    when:
    params.CALCULATE_BARRIER

    """
    cp2nx.py $nhood_file ${params.neighbourhood_module_no} ./
    """
}

process GRAPH_BARRIER {
    /*
    Run graph barrier scoring for all cell types.
    */

    executor "slurm"
    time "6h"
    clusterOptions "--part=gpu --gres=gpu:1"

    module params.md_conda
    conda params.graph_conda_env

    publishDir "${params.outdir}/${params.RELEASE}/graph/barrier", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    val adj_list from adj_output_ch
    val imagename from ch_imagenames.splitText().map{x -> x.trim()}

    output:
    file "**/*.csv" optional true into ch_barrier_results

    when:
    params.CALCULATE_BARRIER

    script:
    if (params.graph_type == 'neighbouRhood')

        """
        stromal_barrier.py --graph_type ${params.graph_type}\
        --neighbourhood_radius 5 \
        --adjacency_data_path $adj_list \
        --root_out ./ \
        --objects_path ${params.OBJECTS} \
        --objects_sep '${params.OBJECTS_DELIMITER}' \
        --panel ${params.PANEL} \
        --calc_chain True \
        --barrier_types Myofibroblasts \
        --phenotyping_level cellType \
        """

    else if (params.graph_type == 'nearest_neighbour')

        """
        stromal_barrier.py --graph_type ${params.graph_type} \
        --imagename $imagename \
        --neighbours 10 \
        --root_out ./ \
        --objects_path ${params.OBJECTS} \
        --objects_sep '${params.OBJECTS_DELIMITER}' \
        --panel ${params.PANEL} \
        --calc_chain True \
        --barrier_types Myofibroblasts \
        --phenotyping_level cellType \
        """

    // """
    // stromal_barrier.py $adj_list ./ ${params.OBJECTS} ${params.PANEL}
    // """

}


process SPATIAL_CLUSTERING {
    /*
    Perform spatial clustering of cell positions.
    */

    executor "slurm"
	time "15m"
	clusterOptions "--part=cpu --cpus-per-task=4 --mem=2GB"

    module params.md_conda
    conda params.spclust_conda_env

    publishDir "${params.outdir}/${params.RELEASE}/spatial_clustering", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    tuple imagename, level from ch_imagenames.splitText().map{x -> x.trim()}.combine(ch_phenotyping).take(2) //split imagenames and remove trailing newline; create a tuple with channel of phenotyping levels

    output:
    file "**/*cluster_assignment.csv" optional true into ch_spclusters
    file "**/*wkt.csv" optional true into ch_wkts
    file "**/*.png" optional true into ch_cluster_plots
    file "**/*.tiff" optional true into ch_alpha_labels

    when:
    params.do_spatial_clustering

    """
    single_cell_spatial_cluster_assignment.py $imagename ${params.COHORT} ${params.PANEL} $level ${params.OBJECTS} '${params.OBJECTS_DELIMITER}' ${params.METADATA} '${params.METADATA_DELIMITER}'
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