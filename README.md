# panregulome-analysis
a Snakemake workflow for exploring the panregulome of grasses with model species Brachypodium distachyon


## REQUIREMENTS
Snakemake
assembly-stats
BUSCO
R
get_homologues
mafft
USAGE
To run this workflow, make sure you have all the required tools installed and the input files available. Then, create a Snakefile and copy the contents of the workflow into it. Adjust the paths and parameters as necessary and run Snakemake.

snakemake --cores <num_cores>


## Installation

### From GitHub
```
git clone https://github.com/ammarabdalrahem/panregulome-analysis.git
```
## Installing Snakemake 
```
pip3 install snakemake
```
## Note
input your username and password in SnakeFile to login JGI

## Usage
```
snakemake --cores all
```
## WORKFLOW
Rule: all
The main rule that generates all the output files.

Rule: make_directories
Creates the necessary output directories.

Rule: obtain_data
Downloads genome sequence data using the provided script and file list.

Rule: uncompressed_files
Unzips and organizes the downloaded data.

Rule: assembly_assessment
Performs statistical analysis on the assembly files and creates an assessment table.

Rule: BUSCO_run
Runs BUSCO analysis on each assembly file and generates short summary files.

Rule: quality_table
Extracts the BUSCO completeness scores from the summary files and creates a quality table.

Rule: assemblies_evaluation
Generates a boxplot to visualize the quality of the assemblies.

Rule: excluded_poor_quality
Identifies and excludes poor quality samples based on the boxplot results.

Rule: repeat_masking_annotation
Downloads and installs the necessary tools for repeat masking and annotation.

Rule: repeats_length
Computes the total length of repeats for each ecotype and merges it with the total genome length.

Rule: repeats_visualization
Generates a plot to visualize the distribution of repeat lengths.

Rule: annotation_analysis
Performs TE annotation analysis by extracting TE families and their lengths from the annotation files.

Rule: annotation_visualization
Generates a plot to visualize the distribution of TE orders.

Rule: promoters_extraction
Extracts proximal promoter sequences from the genome files.

Rule: get_homologues
Downloads and installs the get_homologues tool.

Rule: gene_clustering
Performs gene clustering using the get_homologues tool.

Rule: header_extraction
Extracts the headers of the gene clusters.

Rule: promoter_clustering
Performs clustering of the promoter sequences.

Rule: global_alignment
Performs global sequence alignment on the gene clusters and promoter sequences.

Rule: local_alignment
Performs local sequence alignment on the gene clusters and promoter sequences.

Rule: trimal_alignment
This rule trims the alignment files using the Trimal software.
