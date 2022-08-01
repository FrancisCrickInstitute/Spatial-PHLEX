
process GENERATE_IMAGENAMES {
    /*
    Generate unique imagenames from the cell objects file. 
    */

    module params.md_conda
    conda params.graph_conda_env

    input:
        path objects

    output:
        stdout emit: imagenames
    
    shell:
    '''
    generate_imagenames.py --objects '!{objects}' --delimiter $'!{params.OBJECTS_DELIMITER}'
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

    // if( params.remove('help') ) {
    //     params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
    //     exit 0
    // }
    println "\nExecuting pipeline with the following parameters:"
    params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
    // additional validation here
}
