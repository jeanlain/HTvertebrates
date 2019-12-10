# HTvertebrates

Scripts used in "Horizontal transfer and evolution of transposable elements in vertebrates" by Hua-Hao Zhang, Jean Peccoud, Min-Rui-Xuan Xu, Xiao-Gu Zhang, Clément Gilbert (submitted).

These scripts are publicly available to indicate how parts of the analysis were automated. There is no garanty regarding their use. 
For those who may want to use the pipeline, see below:

## Requirements
- [R](https://cran.r-project.org) 3.4+ 
- [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) 1.0.10
- [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) 4.0.7
- [BUSCO](https://gitlab.com/ezlab/busco) 3.0.1
- [ncbi blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 2.6.0
- [diamond](https://github.com/bbuchfink/diamond) 0.9.19
- [seqtk](https://github.com/lh3/seqtk) 1.2-r94 
- [Slurm Workload Manager](https://slurm.schedmd.com/download.html) 17.11.7

The pipeline was not tested with other versions of the above programs. 

Hardware requirements: a linux cluster with ≥200 CPUs, ≥0.5 TB of system memory,  ≥2 TB of free hard drive space and working internet connection (fiber is recommended). On this hardware, the pipeline should take 1-2 months to complete.

## Installation
Download all the files of this repository into the same directory.

## Usage
Run R scripts whose name start with numbers in the corresponding order, always from the installation directory, which should be set as the working directory.
Adapting this pipeline to other datasets and automating all procedures require modifications to the code. Some parts of the analysis were not automated.

## File description
The R scripts whose names start with numbers performed successive stages of the analysis. The purpose of each script is described by comments at the beginning of the script. 

- HTvFunctions.R and circularPlots.R contain functions required for the other scripts and are sourced automatically.
- The remaining scripts are launched via Rscript (for long, CPU-intensive tasks), from the scripts whose names start with numbers.

The following files are required by the scripts:
- 307.species.info.txt gives general information about the genomes and is used to download genomes sequences from ncbi
- ftp_links.txt contains URL to the genome sequences
- timetree.nwk is the timetree (newick format) used through the analysis
- namedClades.txt is a table of major vertebrate clades in this tree, with their names and color codes used to make some of the paper's figures (these figures are generated with the scripts).
- superF.txt makes the correspondance between repeatModeler family codes (first column), TE class (2nd column) and more common TE superfamily names (3rd column). It is used in steps 15 and 16.
 
## Output
The final output corresponds to results of the publication (please see the publication for their description).
- Figure2.pdf, Figure3.pdf and Figure4.pdf are produced at steps 14, 15 and 16 respectively. They correspond to figures of the main text
- figureS1.pdf is generated at step 5. It corresponds to the supplementary figure 1.
- figureS2.pdf is generated at step 11. It corresponds to the supplementary figure 2.
- tableS1.txt is generated at step 15. It corresponds to the supplementary table 1.
- tableS2.txt is generated at step 16. It corresponds to the supplementary table 2.
- supplementary-data3-TEcomposition_per_species.txt is generated at step 16. 
- supplementary-data4-retained_hits.txt is generated step 12. 
