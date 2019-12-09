

##%######################################################%##
#                                                          #
####           downloading genomes from ncbi            ####
#                                                          #
##%######################################################%##
####### it's probably better to run this script via: Rscript 1-downloadGenomes.R


##### first installing required packages -----------------------------------------------------------
install.packages(c(
  "stringi",
  "data.table",
  "matrixStats",
  "ape",
  "igraph",
  "seqinr",
  "RcolorBrewer"
))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")


#### now downlading genomes -----------------------------------------------------------

source("HTvFunctions.R")

genomeFolder = "genomes/"
dir.create(genomeFolder)
sp = data.table(read.table(
  file = "307.species.info.txt",
  header = F,
  sep = "\t",
  quote = ""
)[, c(1, 9)])	                                        	# fread() doesn't like this file, hence read.table()
links = fread("ftp_links.txt", select = c(1, 20))				#retreived from ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other 
                                                        #and ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian
sp[, V9 := gsub("GCF", "GCA", V9)]
sp[, url := links[chmatch(V9, assembly_accession), ftp_path]]	#retreives genome project urls by matching accession numbers
sp = sp[!is.na(url)]
sp[, folder := stri_c(genomeFolder, gsub(" ", "_", V1))]			#create a folder names where the genome sequences will go, based on species names. Spaces are placed with underscores
sp[, destination := stri_c(folder, "/", basename(url), "_genomic.fna.gz")]	#the final path of each genome sequence, upon completion of the download
sp[, url := stri_c(url, "/", basename(destination))]				#creates complete url of genome compressed fastas
sp[, temp := stri_c(genomeFolder, basename(destination))]		#each genome will be downloaded to this file (and moved after completion)
f = !file.exists(sp$destination)								            #in case this script needs to be re-run, avoids downloading genomes several times
spList = split(sp[f], 1:sum(f))


for (sp in spList) {
  dir.create(sp$folder)
  download.file(sp$url, sp$temp)
  file.rename(sp$temp, sp$destination)					        	#moves genome to destination upon completion
}
