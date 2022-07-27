
/*
Generate unique imagenames from the cell objects file. 
*/
process GENERATE_IMAGENAMES {

    module params.md_conda
    conda params.graph_conda_env

    input:
        path objects



    
    script:
    """
    echo $objects
    """

}



// #!/usr/bin/env python
// import pandas as pd

// objects = pd.read_csv('$objects', sep='$delim', encoding='latin1')
// imagenames = objects['imagename'].unique().tolist()
// for imagename in imagenames:
//     print(imagename)