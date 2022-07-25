
/*********************************************
* RUBICON NEXTFLOW SPATIAL ANALYSIS PIPELINE *
*********************************************/

/*
COHORT PARAMETERS
*/
params.COHORT = 'tx100'
params.PANEL = 'p1'
params.RELEASE = "2022-04-01_tx100_spclustered_into_barrier" //'2022_02_11_release'

/*
 * SPATIAL CLUSTERING CONFIG PARAMETERS
 */

params.do_spatial_clustering = true
params.OBJECTS = "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_${params.PANEL}.txt" //"/camp/lab/swantonc/working/Alastair/graph_analysis/stromal_barrier/simulations/2022-03-28/from_mask/data/all_P1_TMA_REC_20190508-roi_14_BT_25_MTF_0.01_MAF_0.95_MSF_0.05_simulations.csv" // "/camp/lab/swantonc/working/Alastair/graph_analysis/stromal_barrier/simulations/2022-03-27/data/simulated_barrier_objects_2022-03-27.csv"//
params.OBJECTS_DELIMITER = '\t'
params.spclust_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/tf"
params.PHENOTYPING_LEVELS = 'cellType' //'majorType,cellType'
params.MAKE_SPATIAL_CLUSTER_MASKS = true
params.METADATA = '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'
params.METADATA_DELIMITER = '\t'

/*
 * GRAPH ANALYSIS CONFIG PARAMETERS
 */

params.graph_type = 'nearest_neighbour' //neighbouRhood'
params.neighborhood_input = "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/nextflow/${params.PANEL}/publication/*/results/segmentation/*/*/neighbourhood.csv"
params.neighbourhood_module_no = 865
params.md_cuda = "CUDA/10.1.105"
params.md_conda = "Anaconda3" 
params.graph_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-0.18"
params.CALCULATE_BARRIER = true
params.BARRIER_DELIMITER = ','


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
    params.graph_type == 'neighbouRhood'

    """
    cp2nx.py $nhood_file ${params.neighbourhood_module_no} ./
    """
}

process SPATIAL_CLUSTERING {
    /*
    Perform spatial clustering of cell positions.
    */

    executor "slurm"
	time "30m"
	// clusterOptions "--part=cpu --mem=2GB"
    clusterOptions "--part=gpu --gres=gpu:1"


    module params.md_conda
    conda params.spclust_conda_env

    publishDir "${params.outdir}/${params.RELEASE}/spatial_clustering", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    tuple imagename, level from ch_imagenames.splitText().map{x -> x.trim()}.combine(ch_phenotyping).take( params.dev ? params.number_of_inputs : -1 ) //split imagenames and remove trailing newline; create a tuple with channel of phenotyping levels

    output:
    file "**/*cluster_assignment.csv" optional true into ch_spclusters
    file "**/Epithelial cells_clustering/*/*cluster_assignment.csv" optional true into ch_epi_spclusters
    file "**/*wkt.csv" optional true into ch_wkts
    file "**/*.png" optional true into ch_cluster_plots
    file "**/*.tiff" optional true into ch_alpha_labels
    val "$imagename" into ch_imagenames_post_spclust


    when:
    params.do_spatial_clustering

    """
    single_cell_spatial_cluster_assignment.py $imagename ${params.COHORT} ${params.PANEL} $level ${params.OBJECTS} '${params.OBJECTS_DELIMITER}' ${params.METADATA} '${params.METADATA_DELIMITER}'
    """

}

// Duplicate spatial clustering output channels as neighbourhood method requires different input.
ch_imagenames_post_spclust.into {ch_imagenames_post_spclust_b ; ch_imagenames_post_spclust_nb}
ch_epi_spclusters.into {ch_epi_spclusters_b ; ch_epi_spclusters_nb}

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
    val imagename from ch_imagenames_post_spclust_b //ch_imagenames.splitText().map{x -> x.trim()} //
    path spclustered_objects from ch_epi_spclusters_b

    output:
    file "**/*.csv" optional true into ch_barrier_results_b

    when:
    params.CALCULATE_BARRIER

    script:
    // if (params.graph_type == 'neighbouRhood')

    //     """
    //     stromal_barrier.py --graph_type ${params.graph_type}\
    //     --neighbourhood_radius 5 \
    //     --adjacency_data_path $adj_list \
    //     --root_out ./ \
    //     --objects_path ${params.OBJECTS} \
    //     --objects_sep '${params.OBJECTS_DELIMITER}' \
    //     --panel ${params.PANEL} \
    //     --calc_chain True \
    //     --barrier_types Myofibroblasts \
    //     --phenotyping_level cellType \
    //     """

    // else if (params.graph_type == 'nearest_neighbour')

    """
    stromal_barrier.py --graph_type ${params.graph_type} \
    --imagename $imagename \
    --neighbours 10 \
    --root_out ./ \
    --objects_path $spclustered_objects \
    --objects_sep '${params.BARRIER_DELIMITER}' \
    --panel ${params.PANEL} \
    --calc_chain True \
    --barrier_types Myofibroblasts \
    --phenotyping_level cellType \
    """
}

process NEIGHBOURHOOD_BARRIER {
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
    val imagename from ch_imagenames_post_spclust_nb //ch_imagenames.splitText().map{x -> x.trim()} //
    path spclustered_objects from ch_epi_spclusters_nb

    output:
    file "**/*.csv" optional true into ch_barrier_results_nb

    when:
    params.graph_type == 'neighbouRhood'

    script:

    """
    stromal_barrier.py --graph_type ${params.graph_type}\
    --neighbourhood_radius 5 \
    --adjacency_data_path $adj_list \
    --root_out ./ \
    --objects_path $spclustered_objects \
    --objects_sep '${params.BARRIER_DELIMITER}' \
    --panel ${params.PANEL} \
    --calc_chain True \
    --barrier_types Myofibroblasts \
    --phenotyping_level cellType \
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