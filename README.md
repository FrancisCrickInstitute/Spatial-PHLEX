# rubicon-sppipeline
Nextflow pipeline for IMC spatial data methods.

# Inputs

# Workflow options
Spatial-PHLEX can be run in multiple modes:
- 'default'
-- The default pipeline. This performs density-based spatial clustering of all cell types in the objects dataframe. Subsequently, graph-based barrier scoring is executed for a defined, source (e.g. CD8 T cells), target (e.g. Epithelial cells), and barrier (e.g. Myofibroblasts) cell type.
- 'spatial_clustering'
-- Performs only spatial clustering, for all cell types identified in the cell objects dataframe.
- 'barrier_only'
-- Performs graph-based barrier scoring without performing spatial clustering.

# Example Usage

```
nextflow run ./main.nf \
    --BARRIER_DELIMITER '\t' \
    --CALCULATE_BARRIER true \
    --METADATA '/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/master_files/metadata/metadata.tracerx.txt'\
    --METADATA_DELIMITER '\t'\
    --OBJECTS "/camp/project/proj-tracerx-lung/tctProjects/rubicon/tracerx/tx100/imc/outputs/cell_typing/tx100_cell_objects_tx100_publication_p1.txt" \
    --OBJECTS_DELIMITER '\t' \
    --PANEL 'p1' \
    --PHENOTYPING_LEVELS 'cellType' \
    --barrier_phenotyping_level 'cellType' \
    --graph_conda_env "/camp/lab/swantonc/working/Alastair/.conda/envs/rapids-22.02" \
    --graph_type 'nearest_neighbour' \
    --md_conda 'Anaconda3' \
    --outdir '../../results' \
    --publish_dir_mode 'copy' \
    --release '2022-08-30_DSL2_dev' \
    --spclust_conda_env "/camp/lab/swantonc/working/Alastair/.conda/envs/tf" \
    --workflow_name 'default' \
    -resume 
```