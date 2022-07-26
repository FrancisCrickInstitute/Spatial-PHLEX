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
    val nhood_file

    output:
    path "*/*.txt", emit: adj_output_ch

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
    val imagename //from ch_imagenames_post_spclust_b //ch_imagenames.splitText().map{x -> x.trim()} //
    path spclustered_objects //from ch_epi_spclusters_b

    output:
    file "**/*.csv" optional true, emit: ch_barrier_results_b

    script:

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
    val adj_list //from adj_output_ch
    val imagename //from ch_imagenames_post_spclust_nb //ch_imagenames.splitText().map{x -> x.trim()} //
    path spclustered_objects //from ch_epi_spclusters_nb

    output:
    file "**/*.csv" optional true, emit: ch_barrier_results_nb

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