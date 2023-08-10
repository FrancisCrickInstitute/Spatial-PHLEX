process SPATIAL_CLUSTERING {
    /*
    Perform spatial clustering of cell positions.
    */

    tag "${imagename}"
    executor "slurm"
	time "30m"
	// clusterOptions "--part=cpu --mem=2GB"
    clusterOptions "--part=gpu --gres=gpu:1"

    publishDir "${params.outdir}/Spatial-PHLEX/${params.release}/spatial_clustering", mode: params.publish_dir_mode, overwrite: params.overwrite

    input:
        val(objects)
        val(imagename)
        val(level)


    output:
        tuple val(imagename), path("**/*cluster_assignment.csv"), optional: true, emit: ch_all_spclusters
        tuple val(imagename), path("**/${params.barrier_target_cell_type}/*/*cluster_assignment.csv"), emit: ch_target_spclusters, optional: true //MAKE CELL TYPE DERIVED IN NEXTFLOW TO ALLOW FOR RESULTS TO BE PASSED DOWNSTREAM
        path "**/*wkt.csv" optional true //, emit: ch_wkts
        path "**/*.png" optional true//, emit: ch_cluster_plots
        path "**/*.tiff" optional true//, emit: ch_alpha_labels
        val "$imagename"//, emit: ch_imagenames_post_spclust

    

    shell:
        if (params.sampleFile)
            '''
            single_cell_spatial_cluster_assignment.py --imagename !{imagename} \
                --phenotyping_column '!{level}' \
                --phenotype_to_cluster '!{params.phenotype_to_cluster}'\
                --objects_filepath !{objects} \
                --objects_sep $'!{params.objects_delimiter}' \
                --sampleFile_filepath !{params.sampleFile} \
                --sampleFile_sep $'!{params.sampleFile_delimiter}' \
                --root_outdir .
            '''
        else 
            '''
            single_cell_spatial_cluster_assignment.py --imagename !{imagename} \
                --phenotyping_column '!{level}' \
                --phenotype_to_cluster '!{params.phenotype_to_cluster}'\
                --objects_filepath !{objects} \
                --objects_sep $'!{params.objects_delimiter}' \
                --root_outdir .
            '''


}

process INTRACLUSTER_DENSITY {
    /*
    Calculate intracluster densities of each cell type.
    */
        
    executor "slurm"
    time "1h"
    clusterOptions "--part=cpu --mem=4GB"

    publishDir "${params.outdir}/Spatial-PHLEX/${params.release}/spatial_clustering/intracluster_density/${params.phenotyping_column}/${imagename}/${cellType}", mode: params.publish_dir_mode, overwrite: params.overwrite
    
    input:
        tuple val(imagename), val(cellType), path(data)
        
    output:
        path "*intracluster_densities.csv", optional: true
        path "*intracluster_properties.png", optional: true
    
    shell:
        '''
        spcluster_properties.py --outdir . \
            --clustered_data !{data} \
            --phenotyping_column !{params.phenotyping_column} \
            --delimiter $'!{params.objects_delimiter}' \
            --imagename !{imagename}
        '''
}

