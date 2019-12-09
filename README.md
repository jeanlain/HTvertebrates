# HTvertebrates

© Jean Peccoud, 2019

Scripts used in "Horizontal transfer and evolution of transposable elements in vertebrates" by Hua-Hao Zhang, Jean Peccoud, Min-Rui-Xuan Xu, Xiao-Gu Zhang, Clément Gilbert (submitted)
These scripts perform the bulk of the analyses.

The R scripts require packages stringi, data.table, Biostrings, matrixStats, ape, igraph, seqinr and RColorBrewer.
Other programs are required: RepeatModeler, RepeatMasker, BUSCO, ncbi blast+, diamond and seqtk, with their dependancies.
It is recommended to use the latest versions of the programs and packages.

These scripts are not a ready-to-use piece of software that could apply the pipeline to any dataset, as assumptions are made regarding the naming of certain files and about the computer hardware used for the analyses (which varied depending on the task performed, due to various constraints). Also, some early parts of the pipeline are not automated and were run "manually" on each species.
To apply the pipeline, R scripts starting with numbers are run in the order corresponding to their numbers, always from the same working directory. The purpose of each script is described by comments at the beginning of the script. It is recommended to run these scripts in interactive mode (from an R session), as adapting this pipeline to other datasets should require some modifications.

The other R scripts whose names do not start with numbers are either:
- sourced from the numbered scripts: HTvFunctions.R and circularPlots.R contain functions required for the other scripts
- or launched via Rscript (for long, CPU-intensive tasks), from the numbered scripts.

Other files are used by the scripts (in addition to supplementary dataset 4 provided in with the paper, but not in this folder do to file size limitations)
- 307.species.info.txt gives general information about the genomes and is used to download genomes sequences from ncbi
- ftp_links.txt is also used to download genomes sequences from ncbi
- timetree.nwk is the timetree (newick format) used through the analysis
- namedClades.txt is a table of major vertebrate clades in this tree, with their names and color codes used to make some of the paper's figures (these figures are generated with the scripts).
- superF.txt makes the correspondance between repeatModeler family codes (first column), TE class (2nd column) and more common TE superfamily names (3rd column). It is used in steps 15 and 16.
 
