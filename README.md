# panregulome-analysis
a Snakemake workflow for exploring the panregulome of grasses with model species Brachypodium distachyon


## REQUIREMENTS
Snakemake <br />
assembly-stats <br />
BUSCO <br />
R <br />
get_homologues <br />
mafft <br />



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
#### Rule: all <br />
The main rule that generates all the output files.

#### Rule: make_directories<br />
Creates the necessary output directories.

#### Rule: obtain_data<br />
Downloads genome sequence data using the provided script and file list.

#### Rule: uncompressed_files<br />
Unzips and organizes the downloaded data.

#### Rule: assembly_assessment<br />
Performs statistical analysis on the assembly files and creates an assessment table.

#### Rule: BUSCO_run<br />
Runs BUSCO analysis on each assembly file and generates short summary files.

#### Rule: quality_table<br />
Extracts the BUSCO completeness scores from the summary files and creates a quality table.

#### Rule: assemblies_evaluation<br />
Generates a boxplot to visualize the quality of the assemblies.

#### Rule: excluded_poor_quality<br />
Identifies and excludes poor quality samples based on the boxplot results.

#### Rule: repeat_masking_annotation<br />
Downloads and installs the necessary tools for repeat masking and annotation.

#### Rule: repeats_length<br />
Computes the total length of repeats for each ecotype and merges it with the total genome length.

#### Rule: repeats_visualization<br />
Generates a plot to visualize the distribution of repeat lengths.

#### Rule: annotation_analysis<br />
Performs TE annotation analysis by extracting TE families and their lengths from the annotation files.

#### Rule: annotation_visualization<br />
Generates a plot to visualize the distribution of TE orders.

#### Rule: promoters_extraction<br />
Extracts proximal promoter sequences from the genome files.

#### Rule: get_homologues<br />
Downloads and installs the get_homologues tool.

#### Rule: gene_clustering<br />
Performs gene clustering using the get_homologues tool.

#### Rule: header_extraction<br />
Extracts the headers of the gene clusters.

#### Rule: promoter_clustering<br />
Performs clustering of the promoter sequences.

#### Rule: global_alignment<br />
Performs global sequence alignment on the gene clusters and promoter sequences.

#### Rule: local_alignment<br />
Performs local sequence alignment on the gene clusters and promoter sequences.

#### Rule: trimal_alignment<br />
This rule trims the alignment files using the Trimal software.
