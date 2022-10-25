# panregulome-analysis
Exploring the panregulome of grasses with model species Brachypodium distachyon

a Snakemake workflow for exploring the panregulome of grasses with model species Brachypodium distachyon

## The analysis includes:

1- Obtain assemblies sequences and annotation files

2- Quality control of assemblies

3- Repeat masking and annotation

4- Construct the Pan-gene 

5- Extract of Proximal promoters (-500)

6- Sequence alignment of genomic sequences (Global, Local)

7- Calculating Nucleotide diversity

8- Mash distance estimation

## Inputs
Scripts
list of name assemblies 


## Installation

### From GitHub
```
git clone https://github.com/ammarabdalrahem/panregulome-analysis.git
```
## Installing Snakemake 
```
pip3 install snakemake
```
## Usage
```
cd panregulome-analysis
snakemake --cores all
```
