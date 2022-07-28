include { SPATIAL_CLUSTERING } from '../modules/spatial_clustering.nf'
include { GENERATE_IMAGENAMES } from '../modules/util.nf'

workflow SPATIAL_CLUSTERING_WF {

  
    take:
        cell_objects
        phenotyping_level

    

    main:


        GENERATE_IMAGENAMES(cell_objects)
        imagenames = GENERATE_IMAGENAMES.out.imagenames

        imagenames = imagenames.splitText().map{x -> x.trim()}
        imagenames.view()

        SPATIAL_CLUSTERING(cell_objects, imagenames, phenotyping_level)


}




