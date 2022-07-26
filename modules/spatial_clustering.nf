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
    tuple val(imagename), val(level)// from ch_imagenames.splitText().map{x -> x.trim()}.combine(ch_phenotyping).take( params.dev ? params.number_of_inputs : -1 ) //split imagenames and remove trailing newline; create a tuple with channel of phenotyping levels

    output:
    file "**/*cluster_assignment.csv" optional true, emit: ch_spclusters
    file "**/Epithelial cells_clustering/*/*cluster_assignment.csv" optional true, emit: ch_epi_spclusters
    file "**/*wkt.csv" optional true, emit: ch_wkts
    file "**/*.png" optional true, emit: ch_cluster_plots
    file "**/*.tiff" optional true, emit: ch_alpha_labels
    val "$imagename", emit: ch_imagenames_post_spclust


    """
    single_cell_spatial_cluster_assignment.py $imagename ${params.COHORT} ${params.PANEL} $level ${params.OBJECTS} '${params.OBJECTS_DELIMITER}' ${params.METADATA} '${params.METADATA_DELIMITER}'
    """

}