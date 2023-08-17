![TRACERx-PHLEX/Spatial-PHLEX](docs/images/Spatial_PHLEX_logo.png)
# Introduction
Nextflow pipeline for multiplex imaging spatial data methods. 

Spatial PHLEX is a constituent component of [TRACERx-PHLEX](https://github.com/FrancisCrickInstitute/TRACERx-PHLEX), a modular pipeline for key data analysis tasks in multiplexed imaging. For a detailed guide to Spatial PHLEX's capabilities, please see the TRACERx PHLEX [paper]() or its [documentation](https://tracerx-phlex.readthedocs.io/en/main/spatialPHLEX.html).


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

# Example data
Example data is provided in the `data` directory. This data is derived from 5 imaging mass cytometry images from the TRACERx non-small cell lung cancer study, which was first processed with [deep-imcyto](https://github.com/FrancisCrickInstitute/deep-imcyto) to generate segmented cells and measure single cell marker intensities, and then the [TYPEx phenotyping pipeline](https://github.com/FrancisCrickInstitute/TYPEx) to arrive at spatially resolved single cells of discrete types. The data is provided in the following file:
- `data/PHLEX_test_cell_objects.csv`


# Example Usage
First, clone the repository. Then, navigate to the `Spatial-PHLEX` directory.
    
```bash
git clone https://github.com/FrancisCrickInstitute/Spatial-PHLEX
cd Spatial-PHLEX
```

Then load Nextflow and Singularity on your system with e.g. (on SLURM):

```bash
ml Nextflow/22.04.0
ml Singularity/3.6.4
```

Set the Singularity cache directory environment variable, to store the container image:

```bash
export NXF_SINGULARITY_CACHEDIR='./singularity'
```

The Spatial PHLEX pipeline can then be run with the following command:

```bash
nextflow run ./main.nf \
    --workflow_name 'clustered_barrier' \
    --objects "./data/PHLEX_test_cell_objects.csv"\
    --objects_delimiter "\t" \
    --image_id_col "imagename"\
    --phenotyping_column 'majorType'\
    --phenotype_to_cluster 'Epithelial cells'\
    --x_coord_col "centerX"\
    --y_coord_col "centerY"\
    --barrier_phenotyping_column "majorType" \
    --barrier_source_cell_type "CD8 T cells"\
    --barrier_target_cell_type "Epithelial cells"\
    --barrier_cell_type "aSMA+ cells"\
    --n_neighbours 10\
    --outdir "../results" \
    --release 'PHLEX_test' \
    --plot_palette "./assets/PHLEX_test_palette.json" \
    -w "./scratch"\
    -profile {your_nf-core_profile} \
    -resume
```

Replace {your_nf-core_profile} with the profile for your organisation, e.g. `crick`. For more information on nf-core profiles, see the [list of available nf-core configs](https://github.com/nf-core/configs) plus other advice on pipeline configuration [here](https://nf-co.re/usage/configuration#profiles).