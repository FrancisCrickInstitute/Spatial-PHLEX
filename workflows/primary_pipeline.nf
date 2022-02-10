///////////////////////////////////////////////
/* --        INCLUDE SUBWORKFLOWS        -- */
//////////////////////////////////////////////
//Note i just set this up to read from modules as a demonstration, we can chat
//about subworkflows
include { GENERATE_IMAGENAMES;
          NEIGHBOURHOOD_GRAPH;
          GRAPH_BARRIER;
          SPATIAL_CLUSTERING } from '../modules/reafctor_processes.nf'

workflow Primary {
  
  // directly create a list channel for phenotyping levels for combinatoric input of phenotype levels and imagenames
  pheno_list = params.PHENOTYPING_LEVELS?.tokenize(',')
  // ch_phenotyping = Channel.fromList(['cellType', 'majorType'])
  ch_phenotyping = Channel.fromList(pheno_list)
  
  neighbourhood_input_file = file(params.params.neighborhood_input)
  ch_nhood = Channel.fromPath( neighbourhood_input_file, checkIfExists: true )
                    .ifEmpty { exit 1, "Input file not found" }
    
  GENERATE_IMAGENAMES()
  ch_nhood | NEIGHBOURHOOD_GRAPH() | GRAPH_BARRIER() 
  
}