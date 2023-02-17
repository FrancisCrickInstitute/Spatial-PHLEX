![TRACERx-PHLEX/Spatial-PHLEX](docs/images/Spatial_PHLEX_logo.png)
# Introduction
Nextflow pipeline for multiplex imaging spatial data methods.

# Inputs
- `cell_objects.csv`
    - A plaintext, delimited file containing single cell-level coordinate data for a set of images, plus their phenotypic identities.
- `metadata.csv`
    - A plaintext, delimited file containing metadata information about the images in `cell_objects.csv`. To run the pipeline this file must contain, for each image, an image identifier (`'imagename'`), and the width and height in pixels for every image as columns with the header `'image_width'` and `'image_height'`.

# Workflow options
Spatial-PHLEX can be run in multiple modes:
- 'default'
    - The default pipeline. This performs density-based spatial clustering of all cell types in the objects dataframe. Subsequently, graph-based barrier scoring is executed for a defined, source (e.g. CD8 T cells), target (e.g. Epithelial cells), and barrier (e.g. Myofibroblasts) cell type.
- 'spatial_clustering'
    - Performs only spatial clustering, for all cell types identified in the cell objects dataframe.
- 'barrier_only'
    - Performs graph-based barrier scoring without performing spatial clustering.

# Example Usage

```bash
nextflow run ./main.nf \
    --sampleFile "/path/to/sample_data.tracerx.txt"\
    --objects "/path/to/PHLEX_test_data.csv"\
    --phenotyping_column 'majorType' \
    --barrier_phenotyping_column 'majorType' \
    --outdir '../results' \
    --release 'PHLEX_test' \
    --workflow_name 'clustered_barrier' \
    --barrier_source_cell_type 'CD8 T cells'\
    --barrier_target_cell_type 'Epithelial cells'\
    --barrier_cell_type 'aSMA+ cells'\
    -w '/path/to/scratch_directory'\
```