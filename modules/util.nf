
/*
Generate unique imagenames from the cell objects file. 
*/
process GENERATE_IMAGENAMES {

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



