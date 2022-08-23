process SPATIAL_CLUSTERING {
    /*
    Perform spatial clustering of cell positions.
    */

    tag "${imagename}"
    executor "slurm"
	time "30m"
	// clusterOptions "--part=cpu --mem=2GB"
    clusterOptions "--part=gpu --gres=gpu:1"

    publishDir "${params.outdir}/${params.release}/spatial_clustering", mode: params.publish_dir_mode, overwrite: params.overwrite

    input:
        val(objects)
        val(imagename)
        val(level)


    output:
        path "**/*cluster_assignment.csv" optional true //, emit: ch_spclusters
        tuple val(imagename), path("**/${params.barrier_target_cell_type}_clustering/*/*cluster_assignment.csv"), emit: ch_target_spclusters, optional: true //MAKE CELL TYPE DERIVED IN NEXTFLOW TO ALLOW FOR RESULTS TO BE PASSED DOWNSTREAM
        path "**/*wkt.csv" optional true //, emit: ch_wkts
        path "**/*.png" optional true//, emit: ch_cluster_plots
        path "**/*.tiff" optional true//, emit: ch_alpha_labels
        val "$imagename"//, emit: ch_imagenames_post_spclust


    shell:

        '''
        single_cell_spatial_cluster_assignment.py --imagename !{imagename} \
            --phenotyping_level !{level} \
            --objects_filepath !{objects} \
            --objects_sep $'!{params.objects_delimiter}' \
            --metadata_filepath !{params.metadata} \
            --metadata_sep $'!{params.metadata_delimiter}' \
            --root_outdir .
        '''

}


