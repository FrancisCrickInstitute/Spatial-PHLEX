
params.COHORT = 'tx100'
params.PANEL = 'p1'
params.PHENOTYPING_LEVEL =  'cellType'// 'majorType', 'neighborhood10'
params.OBJECTS = '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_p1.txt'
params.outdir = '../../results'
params.publish_dir_mode = 'copy'
params.md_conda = "Anaconda3" 
params.spclust_conda_env = "/camp/lab/swantonc/working/Alastair/.conda/envs/tf"
params.RELEASE = "development_testing_20220329"

ch_phenotyping = Channel.fromList(['cellType', 'majorType'])

process GENERATE_IMAGENAMES {

    /*
    Generate imagenames from the cell objects file. 
    */

    module params.md_conda
    conda params.spclust_conda_env

    output:
    stdout into ch_imagenames
 
    """
    #!/usr/bin/env python
    import pandas as pd
 
    objects = pd.read_csv('${params.OBJECTS}', sep='\t', encoding='latin1')
    imagenames = objects['imagename'].unique().tolist()
    i = 0
    for imagename in imagenames:
        print(imagename)
        i += 1
        if i > 1:
            break 
    """

}



process SPATIAL_CLUSTERING {

    executor "slurm"
	time "0.5h"
	clusterOptions "--part=cpu --cpus-per-task=4 --mem=2GB"

    module params.md_conda
    conda params.spclust_conda_env

    echo true

    publishDir "${params.outdir}/${params.RELEASE}/spatial_clustering", mode: params.publish_dir_mode

    input:
    tuple imagename, level from ch_imagenames.splitText().map{x -> x.trim()}.combine(ch_phenotyping) //.take(5) //split imagenames and remove trailing newline

    output:
    file "**/*cluster_assignment.csv" optional true into ch_spclusters
    file "**/*wkt.csv" optional true into ch_wkts
    file "**/*.png" optional true into ch_cluster_plots

    """
    single_cell_spatial_cluster_assignment.py $imagename ${params.COHORT} ${params.PANEL} $level ${params.OBJECTS}
    """

}

// process CONCAT_CLUSTERS {

//     /*
//     Generate imagenames from the cell objects file. 
//     */

//     publishDir "${params.outdir}/${params.RELEASE}/spatial_clustering", mode: params.publish_dir_mode

//     input: 
//     file '*.csv' from ch_spclusters.toList()
    
//     output:
//     file 'result.txt'
    
//     """
//     cat *.csv > result.csv
//     """ 

// }