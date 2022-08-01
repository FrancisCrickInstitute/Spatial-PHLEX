process NEIGHBOURHOOD_GRAPH {
    /*
    * process to create the spatial graph from the cellprofiler neighbourhood output
    */

    // executor "slurm"
	// time "0.25h"
	// clusterOptions "--part=cpu --mem=2GB"

    module params.md_conda
    conda params.graph_conda_env

    publishDir "${params.outdir}/${params.release}/graph/adjacency_lists/neighbourhood", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
    val nhood_file
    val nhood_module_no

    output:
    path "*/*.txt", emit: adj_output_ch

    """
    cp2nx.py $nhood_file $nhood_module_no .
    """
}


process GRAPH_BARRIER {
    /*
    Run graph barrier scoring for all cell types.
    */
    
    tag "${imagename}"

    executor "slurm"
    time "6h"
    clusterOptions "--part=gpu --gres=gpu:1"

    module params.md_conda
    conda params.graph_conda_env

    publishDir "${params.outdir}/${params.release}/graph/barrier", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
        tuple val(imagename), path(spclustered_objects)
        

    output:
        path "**/*barrier_results.csv", emit: ch_barrier_results, optional: true


    shell:
        '''
        stromal_barrier.py --graph_type !{params.graph_type} \
            --imagename !{imagename} \
            --neighbours 10 \
            --root_out . \
            --objects_path !{spclustered_objects} \
            --objects_sep $'!{params.BARRIER_DELIMITER}' \
            --panel !{params.PANEL} \
            --calc_chain True \
            --barrier_types Myofibroblasts \
            --phenotyping_level !{params.barrier_phenotyping_level} \
        '''
}

process NEIGHBOURHOOD_BARRIER {
    /*
    Run graph barrier scoring for all cell types.
    */

    tag "${imagename}"

    executor "slurm"
    time "6h"
    clusterOptions "--part=gpu --gres=gpu:1"

    module params.md_conda
    conda params.graph_conda_env

    publishDir "${params.outdir}/${params.release}/graph/barrier", mode: params.publish_dir_mode, overwrite: params.OVERWRITE

    input:
        tuple val(imagename), path(adj_list)
        path objects

    output:
        path "**/*.csv" optional true //, emit: ch_barrier_results_nb

    shell:

    '''
    stromal_barrier.py --graph_type !{params.graph_type}\
    --imagename !{imagename} \
    --neighbourhood_radius 5 \
    --adjacency_data_path !{adj_list} \
    --root_out . \
    --objects_path !{objects} \
    --objects_sep $'!{params.OBJECTS_DELIMITER}' \
    --panel !{params.PANEL} \
    --calc_chain True \
    --barrier_types Myofibroblasts \
    --phenotyping_level !{params.barrier_phenotyping_level} \
    '''

}