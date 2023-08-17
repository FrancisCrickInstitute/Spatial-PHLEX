include { NEIGHBOURHOOD_GRAPH; NN_BARRIER} from '../modules/graph_barrier.nf'
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

        NEIGHBOURHOOD_BARRIER(adjacency_mapped, cell_objects)

}

workflow NEAREST_NEIGHBOUR_WF {

  
    take:
        cell_objects

   
    main:

        // Generate imagenames from cell_objects dataframe:
        GENERATE_IMAGENAMES(cell_objects)

        // convert stdout to a channel of a list of strings
        imagenames = GENERATE_IMAGENAMES.out.imagenames
        imagenames = imagenames.splitText().map{x -> x.trim()}
                                .take( params.dev ? params.number_of_inputs : -1 )

        // Pass epithelial spatial clusters to the barrier module:
        NN_BARRIER(cell_objects, imagenames)

}