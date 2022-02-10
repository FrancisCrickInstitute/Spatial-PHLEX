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
    path(nhood_file)

    output:
    path "*/*.txt", emit: adj_output

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
    path(adj_list)

    output:
    file "**/*.csv", emit: into ch_barrier_results

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