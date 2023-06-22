
process GENERATE_IMAGENAMES {
    /*
    Generate unique imagenames from the cell objects file. 
    */

    input:
        path objects

    output:
        stdout emit: imagenames
    
    shell:
    '''
    generate_imagenames.py --objects '!{objects}' --image_id_col !{params.image_id_col} --delimiter $'!{params.objects_delimiter}'
    '''

}


def print_logo(){

    println("")
    println("  __________  ___   ________________          ")
    println(""" /_  __/ __ \\/   | / ____/ ____/ __ \\_  __    """)
    println("  / / / /_/ / /| |/ /   / __/ / /_/ / |/_/    ")
    println(" / / / _, _/ ___ / /___/ /___/ _, _/>  <      ")
    println("/_/ /_/_|_/_/ _|_\\____/_____/_/_|_/_/|_|      ")
    println("   / __ \\/ / / / /   / ____/ |/ /             ")
    println("  / /_/ / /_/ / /   / __/  |   /              ")
    println(" / ____/ __  / /___/ /___ /   |               ")
    println("/_/   /_/ /_/_____/_____//_/|_|               ")    
    println("")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~ Multiplex Imaging Pipeline ~~~~~~~")
    println("~~~~~~~ Part III: Spatial-PHLEX ~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("")

}

def check_params() {
    println("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println "\nExecuting pipeline with the following parameters:"
    params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
    println("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

}


def get_cellType_tuple(ArrayList channel) {
    //gets [sample_id,cellType,path_to_file] for clustering data
    def sample = channel[0]
    def file_list = channel[1]

    if (file_list.getClass() == java.util.ArrayList) {

        def new_array = []
        for (int i=0; i<file_list.size(); i++) {
            def item = []
            item.add(sample)
            item.add(file_list[i].getParent().getParent().getName())
            item.add(file_list[i])
            new_array.add(item)
        }

        return new_array
    }
    else {
        def new_array = []
        new_array.add(sample)
        new_array.add(file_list.getParent().getParent().getName())
        new_array.add(file_list)
        return new_array
    }
    
}

def group_channel(x){
    grouped = x.map { get_cellType_tuple(it) }
                .flatten()
                .collate(3)
    return grouped
}

// def group_channel(x){
//     grouped = x.map { get_cellType_tuple(it) }
//                 .flatten()
//                 .collate(3)
//                 .groupTuple(by: [0,1])
//                 .map { it -> [ it[0], it[1], it[2].sort() ] }
//     return grouped
// }