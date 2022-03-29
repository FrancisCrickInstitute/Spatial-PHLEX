#!/usr/bin/env nextflow

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
params.RELEASE = 'development_testing'
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
params.number_of_inputs = 1
params.publish_dir_mode = 'copy'
params.OVERWRITE = true
params.outdir = '../../results'
project_dir = projectDir

// directly create a list channel for phenotyping levels for combinatoric input of phenotype levels and imagenames
pheno_list = params.PHENOTYPING_LEVELS?.tokenize(',')
// ch_phenotyping = Channel.fromList(['cellType', 'majorType'])
ch_phenotyping = Channel.fromList(pheno_list)

// WKT input

params.wkt_input = "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/spatial/results/2022_02_09_release/spatial_clustering/single_cell_assignment/tx100/p1/*/dbscan_25/min_size_0/alpha_0.05/*/*/*_alphashape_polygons_wkt.csv"
params.MAKE_SPATIAL_CLUSTER_MASKS = true
params.METADATA = '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'
params.METADATA_DELIMITER = '\t'


// // channel for neighbourhood input csv files
// if (params.wkt_input) {
//     Channel
//         .fromPath(params.wkt_input, checkIfExists: true)
//         .take( params.dev ? params.number_of_inputs : -1 )
//         .map { it }
//         .ifEmpty { exit 1, "Input file not found: ${params.wkt_input}" }
//         .set { ch_wkts }
// } else {
//    exit 1, "Neighbourhood input file not specified!"
// }

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

process SPATIAL_CLUSTERING {
    /*
    Perform spatial clustering of cell positions.
    */

    // executor "slurm"
	// time "5m"
	// clusterOptions "--part=cpu --cpus-per-task=4 --mem=2GB"

    module params.md_conda
    conda params.spclust_conda_env

    echo true

    publishDir "${params.outdir}/${params.RELEASE}/spatial_clustering", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    tuple imagename, level from ch_imagenames.splitText().map{x -> x.trim()}.combine(ch_phenotyping).take(2) //split imagenames and remove trailing newline; create a tuple with channel of phenotyping levels

    output:
    file "**/*cluster_assignment.csv" optional true into ch_spclusters
    file "**/*wkt.csv" optional true into ch_wkts
    file "**/*.png" optional true into ch_cluster_plots
    file "**/*.tiff" optional true into ch_alpha_labels

    """
    single_cell_spatial_cluster_assignment.py $imagename ${params.COHORT} ${params.PANEL} $level ${params.OBJECTS} '${params.OBJECTS_DELIMITER}' ${params.METADATA} '${params.METADATA_DELIMITER}'
    """

}

process CONVERT_WKT {
    /*
    * process to create the spatial graph from the cellprofiler neighbourhood output
    */

    // executor "slurm"
	// time "0.25h"
	// clusterOptions "--part=cpu --cpus-per-task=4 --mem=2GB"

    when:
    params.MAKE_SPATIAL_CLUSTER_MASKS

    module params.md_conda
    conda params.spclust_conda_env

    // publishDir "${params.outdir}/${params.RELEASE}/graph/adjacency_lists/neighbourhood", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    val wkt_file from ch_wkts

    // output:
    // stdout
    // path "*/*.txt" into adj_output_ch

    """
    convert_wkt_to_tiff.py '${wkt_file}' ${params.METADATA} '${params.METADATA_DELIMITER}'
    """
}