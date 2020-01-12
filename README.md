# HTvertebrates

Scripts used in "Horizontal transfer and evolution of transposable elements in vertebrates" by Hua-Hao Zhang, Jean Peccoud, Min-Rui-Xuan Xu, Xiao-Gu Zhang, Clément Gilbert (submitted).

These scripts are publicly available to indicate how parts of the analysis were automated. There is no guaranty regarding their use. 
For those who may want to use the pipeline, see below:

## Requirements
- [R](https://cran.r-project.org) 3.4+ (required R packages are installed automatically)
- [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) 1.0.10
- [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) 4.0.7
- [BUSCO](https://gitlab.com/ezlab/busco) 3.0.1
- [ncbi blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) 2.6.0
- [diamond](https://github.com/bbuchfink/diamond) 0.9.19
- [seqtk](https://github.com/lh3/seqtk) 1.2-r94 
- [Slurm Workload Manager](https://slurm.schedmd.com/download.html) 17.11.7

The pipeline was not tested with other versions of the above programs. 

Hardware requirements: a linux cluster with ≥200 CPUs, ≥0.5 TB of system memory,  ≥2 TB of free hard drive space and working internet connection (fiber is recommended). On this hardware, the pipeline should take 1-2 months to complete.

## File description
The R scripts whose names start with numbers performed successive stages of the analysis. The purpose of each script is described by comments at the beginning of the script. 

- HTvFunctions.R contains functions required for the other scripts and is sourced automatically.
- circularPlots.R contains functions used to draw figure 2. 
- The remaining scripts are launched via Rscript (for long, CPU-intensive tasks) from the scripts whose names start with numbers. 

The files in directory "additonal_files" are required by the scripts:
- supplementary-data1-genomes_and_accessions.txt gives general information about the genomes and is used to download genomes sequences from ncbi
- ftp_links.txt contains URL to the genome sequences
- timetree.nwk is the timetree (newick format) used through the analysis
- namedClades.txt is a table of major vertebrate clades in this tree, with their names and color codes used to make some of the paper's figures (these figures are generated with the scripts).
- superF.txt makes the correspondance between repeatModeler family codes (first column), TE class (2nd column) and more common TE superfamily names (3rd column). It is used in steps 15 and 16.
- supplementary-data3-TEcomposition_per_species.txt is generatd by the scripts and is provided with the paper, but we also provide it here if to facilitate the reproduction of the results.

The directory demo_TeKaKs is provided to demo the script TEKaKs.R (see "Demonstration of TEKaKs.R" below), but is not required to run the pipeline.

## Installation
In a bash-compatible terminal, paste
```
git clone https://github.com/jeanlain/HTvertebrates.git
cd HTvertebrates/
```

## Usage
Run R scripts whose name start with numbers in the corresponding order, always from the HTvertebrates/ directory, which should be set as the working directory.

Adapting this pipeline to other datasets, hardware configuration, and automating all procedures require modifications to the code. Some parts of the analysis were not automated.

### Demonstration of TEKaKs.R 
We detail how to run "TEKaKs.R" on a demo dataset, but we remind that this script (as all others) is not intended for use in any other context than the study associated with the paper.

The hardware requirement for this demo is a computer with at last 8GB of RAM, 1GB of free hard drive space, and which is able to execute R 3.4+ in a terminal. The other programs mentioned in the Requirements section need not be installed for this demo.

The demo_TeKaKs/ directory must be immediately within the HTvertebrates/ directory. It contains the following:
- "TEhitFile.txt" is a file of TE-TE HSPs in typical blast tabular format.
- "blastxFile.txt" is a tabular file of TE-protein HSPs. The fields indicate the TE sequence name, start and end coordinates of the HSP on this sequence (where start < end), start coordinate of the HSP on the protein and whether the TE sequence in aligned on the protein in reverse direction.
- "fastaFile.fas" is a fasta file of the TE sequences whose names are in the two previous file.

The nature of these files is also explained in comments in TEKaKs.R

To run the demo, paste the following in the terminal session that you used to install the pipeline:
```
Rscript TEKaKs.R demo_TeKaKs/TEhitFile.txt demo_TeKaKs/blastxFile.txt demo_TeKaKs/fastaFile.fas demo_TeKaKs/output 2
```
where "demo_TeKaKs/output" is the output folder (automatically created) and "2" is the number of CPU to use.

The script should run in less than 5 minutes on a current desktop PC.

Results will be found in "demo_TeKaKs/output".

"allKaKs.txt" is a tabular file that contains the following fields:
- "hit" is an identifier for each HSP, which corresponds to the row index of each HSP in "TEhitFile.txt".
- "ka", "ks", "vka" and "vks" are the results of Ka and Ks computations (see the kaks() function of seqinr), 
- "length" is the length of the alignment on which the above were computed.
- "nMut" is the number of mutations in this alignment.
- "K80distance" and "rawDistance" are molecular distances (according to Kimura 1980 or without any correction) between sequences in the HSP. These are computed before any of the processing required for the Ka Ks computations.


## Output of the pipeline
More than 1TB of intermediate files are generated.
The final output corresponds to results of the publication (please see the publication for their description).
- Figure2.pdf, Figure3.pdf and Figure4.pdf are produced at steps 14, 15 and 16 respectively. They correspond to figures of the main text
- figureS1.pdf is generated at step 5. It corresponds to the supplementary figure 1.
- figureS2.pdf is generated at step 11. It corresponds to the supplementary figure 2.
- figureS5.pdf, figureS6.pdf and figureS7.pdf are generated at step 16. They corresponds to the supplementary figures 5-7.
- tableS1.txt is generated at step 15. It corresponds to the supplementary table 1.
- tableS2.txt is generated at step 16. It corresponds to the supplementary table 2.
- supplementary-data3-TEcomposition_per_species.txt is generated at step 2. 
- supplementary-data4-retained_hits.txt is generated step 12. 

