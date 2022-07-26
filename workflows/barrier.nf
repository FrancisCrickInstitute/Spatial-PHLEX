include { NEIGHBOURHOOD_GRAPH; NEIGHBOURHOOD_BARRIER } from '../modules/graph_barrier.nf'


workflow NEIGHBOURHOOD_WF {

  
    take:
        nhood_file
        imagename
        cell_objects
    

    main:
        
        //Run the graph barrier algorithm with neighrbouRhood graph constructed from the cell objects:
        NEIGHBOURHOOD_GRAPH(nhood_file)

        NEIGHBOURHOOD_BARRIER(NEIGHBOURHOOD_GRAPH.adj_output_ch, imagename, cell_objects)

}