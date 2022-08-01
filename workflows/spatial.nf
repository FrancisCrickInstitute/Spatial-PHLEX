include { SPATIAL_CLUSTERING } from '../modules/spatial_clustering.nf'
include { GENERATE_IMAGENAMES } from '../modules/util.nf'
include { GRAPH_BARRIER } from '../modules/graph_barrier.nf'

workflow SPATIAL_CLUSTERING_WF {

    /*
    * Workflow to perform only spatial clustering of cell objects.
    */

    take:
        cell_objects
        phenotyping_level

    
    main:

        // Generate imagenames from cell_objects dataframe:
        GENERATE_IMAGENAMES(cell_objects)

        // convert stdout to a channel of a list of strings
        imagenames = GENERATE_IMAGENAMES.out.imagenames
        imagenames = imagenames.splitText().map{x -> x.trim()}

        // Perform spatial clustering:
        SPATIAL_CLUSTERING(cell_objects, imagenames, phenotyping_level)

}

workflow CLUSTERED_BARRIER_WF {

    /*
    * Workflow to perform spatial clustering of cell objects and pass epithelial spatial clusters to the barrier module.
    */
    

    take:
        cell_objects
        phenotyping_level

    
    main:

        // Generate imagenames from cell_objects dataframe:
        GENERATE_IMAGENAMES(cell_objects)

        // convert stdout to a channel of a list of strings
        imagenames = GENERATE_IMAGENAMES.out.imagenames
        imagenames = imagenames.splitText().map{x -> x.trim()}
                                .take( params.dev ? params.number_of_inputs : -1 )

        // Perform spatial clustering:
        SPATIAL_CLUSTERING(cell_objects, imagenames, phenotyping_level)

        // Pass epithelial spatial clusters to the barrier module:
        SPATIAL_CLUSTERING.out.ch_epi_spclusters.view()

        GRAPH_BARRIER(SPATIAL_CLUSTERING.out.ch_epi_spclusters)

        GRAPH_BARRIER.out.ch_barrier_results.view()

}




