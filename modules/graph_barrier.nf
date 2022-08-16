process NEIGHBOURHOOD_GRAPH {
    /*
    * process to create the spatial graph from the cellprofiler neighbourhood output
    */

    // executor "slurm"
	// time "0.25h"
	// clusterOptions "--part=cpu --mem=2GB"

    module params.md_conda
    conda params.graph_conda_env

    publishDir "${params.outdir}/${params.release}/graph/adjacency_lists/neighbourhood", mode: params.publish_dir_mode, overwrite: params.overwrite

    input:
        val nhood_file
        val nhood_module_no

    output:
        path "*/*.txt", emit: adj_output_ch
    
    script:
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

    // module params.md_conda
    // conda params.graph_conda_env
    module params.md_cuda

    publishDir "${params.outdir}/${params.release}/graph/clustered_barrier", mode: params.publish_dir_mode, overwrite: params.overwrite

    input:
        tuple val(imagename), path(spclustered_objects)
        

    output:
        path "**/*barrier_results.csv", emit: ch_barrier_results, optional: true


    shell:
        '''
        stromal_barrier.py --graph_type !{params.graph_type} \
            --source_cell_type '!{params.barrier_source_cell_type}' \
            --target_cell_type '!{params.barrier_target_cell_type}' \
            --imagename !{imagename} \
            --neighbours !{params.n_neighbours} \
            --root_out . \
            --objects_path !{spclustered_objects} \
            --objects_sep $'!{params.objects_delimiter}' \
            --panel !{params.PANEL} \
            --barrier_types '!{params.barrier_cell_type}' \
            --phenotyping_level !{params.barrier_phenotyping_level} \
        '''
}


process NN_BARRIER {
    /*
    Run graph barrier scoring without spatial clustering with nearest neighbour graph construction.
    */
    
    tag "${imagename}"

    executor "slurm"
    time "6h"
    clusterOptions "--part=gpu --gres=gpu:1"

    module params.md_conda
    conda params.graph_conda_env

    publishDir "${params.outdir}/${params.release}/graph/unclustered_barrier", mode: params.publish_dir_mode, overwrite: params.overwrite

    input:
        path objects
        val imagename 
        

    output:
        path "**/*barrier_results.csv", emit: ch_barrier_results, optional: true


    shell:
        '''
        stromal_barrier.py --graph_type !{params.graph_type} \
            --source_cell_type '!{params.barrier_source_cell_type}' \
            --target_cell_type '!{params.barrier_target_cell_type}' \
            --imagename !{imagename} \
            --neighbours !{params.n_neighbours} \
            --root_out . \
            --objects_path !{objects} \
            --objects_sep $'!{params.objects_delimiter}' \
            --panel !{params.PANEL} \
            --barrier_types '!{params.barrier_cell_type}' \
            --phenotyping_level !{params.barrier_phenotyping_level} \
        '''
}

process NEIGHBOURHOOD_BARRIER {
    /*
    Run graph barrier scoring without spatial clustering with neighbouRhood graph construction.
    */

    tag "${imagename}"

    executor "slurm"
    time "6h"
    clusterOptions "--part=gpu --gres=gpu:1"

    module params.md_conda
    conda params.graph_conda_env

    publishDir "${params.outdir}/${params.release}/graph/unclustered_barrier", mode: params.publish_dir_mode, overwrite: params.overwrite

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
        --objects_sep $'!{params.objects_delimiter}' \
        --panel !{params.PANEL} \
        --calc_chain True \
        --barrier_types Myofibroblasts \
        --phenotyping_level !{params.barrier_phenotyping_level} \
        '''

}

process AGGREGATE_SCORES {
    /*
    Aggregate the scores from the barrier scoring.
    */
        
    executor "slurm"
    time "6h"
    clusterOptions "--part=gpu --gres=gpu:1"
    
    module params.md_conda
    conda params.graph_conda_env
    
    publishDir "${params.outdir}/${params.release}/graph/aggregated_barrier_scoring", mode: params.publish_dir_mode, overwrite: params.overwrite
    
    input:
        path scores
        
    output:
        path "*barrier_summary.csv", emit: ch_aggregated_barrier_results, optional: true
    
    shell:
        '''
        agg_barrier_scores.py --outdir . \
            --barrier_scores !{scores} \
            --delimiter $'!{params.objects_delimiter}' 
        '''
}

