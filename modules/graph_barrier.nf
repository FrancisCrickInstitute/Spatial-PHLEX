// process NEIGHBOURHOOD_GRAPH {
//     /*
//     * process to create the spatial graph from the cellprofiler neighbourhood output
//     */
//     executor "slurm"
// 	time "0.5h"
// 	clusterOptions "--part=cpu --cpus-per-task=4 --mem=2GB"

//     module params.md_conda
//     conda params.spclust_conda_env

//     echo true

//     publishDir "${params.outdir}/graph/adjacency_lists", mode: params.publish_dir_mode

//     input:
//     val nhood_file from ch_nhood

//     """
//     cp2nx.py $nhood_file ${params.neighbourhood_module_no} './'
//     """
// }

// process GRAPH_BARRIER {

//     executor "slurm"
// 	time "6h"
// 	clusterOptions "--part=gpu --gres=gpu:1"

//     module params.md_conda
//     conda params.graph_conda_env

//     echo true

//     publishDir "${params.outdir}", mode: params.publish_dir_mode

//     input:
//     val imagename from ch_imagenames

//     // output:
//     // file "*/*/*/*/*/*.csv"
//     // val preprocessdir into ch_preprocess_results


//     // "echo $imagename ${params.outdir}"
    

//     """
//     stromal_barrier_batch.py $imagename ${params.outdir}
//     """

// }