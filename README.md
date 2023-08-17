![TRACERx-PHLEX/Spatial-PHLEX](docs/images/Spatial_PHLEX_logo.png)
# Introduction
Nextflow pipeline for multiplex imaging spatial data methods. 

Spatial PHLEX is a constituent component of [TRACERx-PHLEX](), a modular pipeline for key data analysis tasks in multiplexed imaging. For a detailed guide to Spatial PHLEX's capabilities, please see the TRACERx PHLEX [paper]() or its [documentation](https://tracerx-phlex.readthedocs.io/en/main/spatialPHLEX.html).


# Requirements
Spatial PHLEX was developed for HPC systems with both CPU and GPU compute facilities. This is because Spatial PHLEX cellular barrier scoring requires GPU-accelerated [graph algorithms](https://github.com/rapidsai/cugraph) in order to compute the required shortest paths analysis in a reasonable time. If you do not have GPUs on your system then you can still run Spatial PHLEX in `spatial_clustering` mode. See [Workflow Options](#workflow-options).

# Inputs
- `cell_objects.csv` | `cell_objects.txt`
    - A plaintext, delimited file containing single cell-level coordinate data for a set of images, plus their phenotypic identities.
- `metadata.csv` (optional)
    - A plaintext, delimited file containing metadata information about the images in `cell_objects.csv`. To run the pipeline this file must contain, for each image, an image identifier (`'imagename'`), and the width and height in pixThis option is preferred if it is desirable for mask and plot outputs to match the input image shape exactly.

# Workflow options
Spatial-PHLEX can be run in multiple modes:
- 'default', 'clustered_barrier'
    - The default pipeline. This performs density-based spatial clustering of all cell types in the objects dataframe. Subsequently, graph-based barrier scoring is executed for a defined, source (e.g. CD8 T cells), target (e.g. Epithelial cells), and barrier (e.g. Myofibroblasts) cell type.
- 'spatial_clustering'
    - Performs only spatial clustering, for all cell types identified in the cell objects dataframe.
- 'barrier_only'
    - Performs graph-based barrier scoring without performing spatial clustering.

# Example Usage

```bash
nextflow run ./main.nf \
    --objects "./data/mydata.csv"\
    --objects_delimiter ','\
    --image_id_col "Image_ID"\
    --x_id "Location_Center_X"\
    --y_id "Location_Center_Y"\
    --barrier_phenotyping_column "Phenotype" \
    --outdir "../results" \
    --release 'PHLEX_example' \
    --workflow_name 'default' \
    --barrier_source_cell_type "CD8 T cells"\
    --barrier_target_cell_type "Epithelial cells"\
    --barrier_cell_type "aSMA+ Fibroblasts"\
    --singularity_bind_path '/camp,/nemo'\
    --n_neighbours 5\
    -w './scratch'\
    -profile {your_nf-core_profile}\
    -resume
```