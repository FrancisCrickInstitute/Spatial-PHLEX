include { SPATIAL_CLUSTERING; INTRACLUSTER_DENSITY } from '../modules/spatial_clustering.nf'
include { GENERATE_IMAGENAMES; group_channel } from '../modules/util.nf'
include { GRAPH_BARRIER; AGGREGATE_SCORES } from '../modules/graph_barrier.nf'

workflow SPATIAL_CLUSTERING_WF {

    /*
    * Workflow to perform only spatial clustering of cell objects.
    */

    take:
        cell_objects
        phenotyping_column

    
    main:

        // Generate imagenames from cell_objects dataframe:
        GENERATE_IMAGENAMES(cell_objects)

        // convert stdout to a channel of a list of strings
        imagenames = GENERATE_IMAGENAMES.out.imagenames
        imagenames = imagenames.splitText().map{x -> x.trim()}
                                .take( params.dev ? params.number_of_inputs : -1 )

        // Perform spatial clustering:
        SPATIAL_CLUSTERING(cell_objects, imagenames, phenotyping_column)
        all_spclusters = group_channel(SPATIAL_CLUSTERING.out.ch_all_spclusters)
        all_spclusters.view()

        all_spclusters.collectFile( name:"$workDir/collected_spclusters.csv", keepHeader: true, skip:1 ).set {all_spclusts}
        all_spclusts.view()

        // INTRACLUSTER_DENSITY(all_spclusters)

}

workflow CLUSTERED_BARRIER_WF {

    /*
    * Workflow to perform spatial clustering of cell objects and pass epithelial spatial clusters to the barrier module.
    */
    

    take:
        cell_objects
        phenotyping_column

    
    main:

        // Generate imagenames from cell_objects dataframe:
        GENERATE_IMAGENAMES(cell_objects)

        // convert stdout to a channel of a list of strings
        imagenames = GENERATE_IMAGENAMES.out.imagenames
        imagenames = imagenames.splitText().map{x -> x.trim()}
                                .take( params.dev ? params.number_of_inputs : -1 )

        // Perform spatial clustering:
        SPATIAL_CLUSTERING(cell_objects, imagenames, phenotyping_column)
        
        SPATIAL_CLUSTERING.out.ch_target_spclusters.view()
    
        // Pass epithelial spatial clusters to the barrier module:
        GRAPH_BARRIER(SPATIAL_CLUSTERING.out.ch_target_spclusters)

        GRAPH_BARRIER.out.ch_barrier_results.collectFile( name:"$workDir/${params.barrier_source_cell_type}_to_${params.barrier_target_cell_type}_barrier_scores.csv", keepHeader: true, skip:1 ).set {barriers}

        barriers.view()

        AGGREGATE_SCORES(barriers)

        all_spclusters = group_channel(SPATIAL_CLUSTERING.out.ch_all_spclusters)
        // INTRACLUSTER_DENSITY(all_spclusters)


}

workflow CLUSTERED_NEIGHBORHOOD_WF {

    /*
    * Workflow to perform spatial clustering of cell objects and pass epithelial spatial clusters to the barrier module.
    */
    

    take:
        cell_objects
        phenotyping_column
        nhood_file
        nhood_module_no
        cell_objects
        objects_delimiter

    
    main:

        // Generate imagenames from cell_objects dataframe:
        GENERATE_IMAGENAMES(cell_objects)

        // convert stdout to a channel of a list of strings
        imagenames = GENERATE_IMAGENAMES.out.imagenames
        imagenames = imagenames.splitText().map{x -> x.trim()}
                                .take( params.dev ? params.number_of_inputs : -1 )

        // Perform spatial clustering:
        SPATIAL_CLUSTERING(cell_objects, imagenames, phenotyping_column)
        all_spclusters = group_channel(SPATIAL_CLUSTERING.out.ch_all_spclusters)
        INTRACLUSTER_DENSITY(all_spclusters)

    
        // Pass epithelial spatial clusters to the barrier module:
        GRAPH_BARRIER(SPATIAL_CLUSTERING.out.ch_target_spclusters)


}




