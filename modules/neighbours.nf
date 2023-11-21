process NEAREST_NEIGHBOURS {
    /*
    Run graph barrier scoring for all cell types.
    */
    
    label 'Spatial_PHLEX_CPU'

    publishDir "${params.outdir}/Spatial-PHLEX/${params.release}/neighbours/nearest_neighbours", mode: params.publish_dir_mode, overwrite: params.overwrite

    input:
        path(objects)
        

    output:
        path "**/*.csv", emit: ch_nearest_neighbours, optional: true


    shell:
        '''
        spatial_phlex_neighbours.py\
            --measurements !{objects}\
            --measurements_sep $'!{params.objects_delimiter}'\
            --n_neighbours !{params.n_neighbours}\
            --object_identifier !{params.object_id_col}\
            --image_id_col !{params.image_id_col}\
            --x_id_col !{params.x_coord_col}\
            --y_id_col !{params.y_coord_col}\
            --phenotyping_ids !{params.phenotyping_column}\
            --outdir './nearest_neighbours'

        '''
}

process COMMUNITIES {
    /*
    Calculate communities through kmeans clustering.
    */

    label 'Spatial_PHLEX_CPU'

    publishDir "${params.outdir}/Spatial-PHLEX/${params.release}/communities/n=${params.n_communities}", mode: params.publish_dir_mode, overwrite: params.overwrite

    input:
        path(objects)

    output:
        path "**/*.csv", emit: ch_communities, optional: true
        path "**/*.png", optional: true
        path "**/*.pdf", optional: true

    shell:

        '''
        spatial_phlex_communities.py\
            --input !{objects}\
            --communities_type !{params.phenotyping_column}\
            --n_communities !{params.n_communities}\
            --image_id_col !{params.image_id_col}\
            --x_id_col !{params.x_coord_col}\
            --y_id_col !{params.y_coord_col}
        '''
}