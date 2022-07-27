include { NEIGHBOURHOOD_GRAPH; NEIGHBOURHOOD_BARRIER } from '../modules/graph_barrier.nf'
include { GENERATE_IMAGENAMES } from '../modules/util.nf'

def generate_imagenames(file){
    def new_array = []
    imagename = file.name.lastIndexOf('.').with {it != -1 ? file.name[0..<it] : file.name}
    new_array.add(imagename)
    new_array.add(file)
    return new_array
}

workflow NEIGHBOURHOOD_WF {

  
    take:
        nhood_file
        nhood_module_no
        cell_objects
        objects_delimiter
    

    main:
        
        //Run the graph barrier algorithm with neighrbouRhood graph constructed from the cell objects:

        // GENERATE_IMAGENAMES(cell_objects)

        NEIGHBOURHOOD_GRAPH(nhood_file, nhood_module_no)

        adjacency_mapped = NEIGHBOURHOOD_GRAPH.out.adj_output_ch.map { generate_imagenames(it) }
        adjacency_mapped.view()

        NEIGHBOURHOOD_BARRIER(adjacency_mapped, cell_objects)

}

workflow NEAREST_NEIGHBOUR_WF {

  
    take:
        nhood_file
        imagename
        cell_objects
    

    main:
        
        //Run the graph barrier algorithm with neighrbouRhood graph constructed from the cell objects:
        NEIGHBOURHOOD_GRAPH(nhood_file)

        NEIGHBOURHOOD_BARRIER(NEIGHBOURHOOD_GRAPH.adj_output_ch, imagename, cell_objects)

}