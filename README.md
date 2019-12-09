# HTvertebrates

© Jean Peccoud, 2019

Scripts used in "Horizontal transfer and evolution of transposable elements in vertebrates" by Hua-Hao Zhang, Jean Peccoud, Min-Rui-Xuan Xu, Xiao-Gu Zhang, Clément Gilbert (submitted).
These scripts perform the bulk of the analyses.

These scripts come with no warranty. They are not a ready-to-use piece of software that could apply the pipeline to any dataset, as assumptions are made regarding the naming of certain files and about the computer hardware used for the analyses (which varied depending on the task performed, due to various constraints). Also, some early parts of the pipeline are not automated and were run "manually" on each species.


# Requirements
The scripts require R 3.4+ with packages stringi, data.table, Biostrings, matrixStats, ape, igraph, seqinr and RColorBrewer.
The following programs are also required: 
RepeatModeler 1.0.10
RepeatMasker 4.0.7
BUSCO 3.0.1
ncbi blast+ 2.6.0
diamond 0.9.19
seqtk 1.2-r94 


# Installation
Dowload all the files of this repository in the same directory.

# Usage
Run R scripts starting with numbers are run in the order corresponding to their numbers, always from the installation directory, which should be set as the working directory.
It is strongly recommended to run these scripts in interactive mode (from an R session), as adapting this pipeline to other datasets should require some modifications. 

# File description
The R scripts whose names start with numbers perform successive stages of the analysis. The purpose of each script is described by comments at the beginning of the script. 

- HTvFunctions.R and circularPlots.R contain functions required for the other scripts and are sourced automatically.
- The remaining scripts are launched via Rscript (for long, CPU-intensive tasks), from the scripts whose names start with numbers.

The following files are required by the scripts:
- 307.species.info.txt gives general information about the genomes and is used to download genomes sequences from ncbi
- ftp_links.txt is also used to download genomes sequences from ncbi
- timetree.nwk is the timetree (newick format) used through the analysis
- namedClades.txt is a table of major vertebrate clades in this tree, with their names and color codes used to make some of the paper's figures (these figures are generated with the scripts).
- superF.txt makes the correspondance between repeatModeler family codes (first column), TE class (2nd column) and more common TE superfamily names (3rd column). It is used in steps 15 and 16.
 
# Ouput
The final output correspond to results of the publication (please see the publication for their description).
- Figure2.pdf, Figure3.pdf and Figure4.pdf are produced at steps 14, 15 and 16 respectively. They correspond to figures of the main text
- figureS1.pdf is generated at step 5. It corresponds to the supplementary figure 1.
- figureS2.pdf is generated at step 11. It corresponds to the supplementary figure 2.
- tableS1.txt is generated at step 15. It corresponds to the supplementary table 1.
- tableS2.txt is generated at step 16. It corresponds to the supplementary table 2.
- supplementary-data3-TEcomposition_per_species.txt is generated at step 16. 
- supplementary-data4-retained_hits.txt is generated step 12. 


