

## %######################################################%##
#                                                          #
##  This stage downloads the species's genomes from ncbi   #
#                                                          #
## %######################################################%##

# it's probably better to run this script via: Rscript 1-downloadGenomes.R



source("HTvFunctions.R")

genomeFolder <- "genomes/"
dir.create(genomeFolder)

# we import the table of species description
# fread() doesn't like this file, hence the use of read.table()
genomesInfo <- fread("supplementary-data1-genomes_and_accessions.txt" )

# we replace spaces in column names
setnames(genomesInfo, gsub(" ","_", names(genomesInfo), fixed = T))

# we import genome urls retrieved from ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other
# and ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian
links <- fread("ftp_links.txt", select = c(1, 20))

genomesInfo[, Assembly_Accession := gsub("GCF", "GCA", Assembly_Accession)]

# we retrieve genome project urls by matching accession numbers
genomesInfo[, url := links[chmatch(Assembly_Accession, assembly_accession), ftp_path]]
genomesInfo <- genomesInfo[!is.na(url)]

# create a folder names where the genome sequences will go, based on species names. Spaces are replaced with underscores
genomesInfo[, folder := stri_c(genomeFolder, gsub(" ", "_", Organism))]


# the final path of each genome sequence, upon completion of the download
genomesInfo[, destination := stri_c(folder, "/", basename(url), "_genomic.fna.gz")]

# creates complete url of genome compressed fastas
genomesInfo[, url := stri_c(url, "/", basename(destination))]

# each genome will be downloaded to this file (and moved after completion)
# this will allow telling whether a genome has been completely downloaded
genomesInfo[, temp := stri_c(genomeFolder, basename(destination))]

# this avoids downloading genomes several times, in case this script needs to be re-run,
f <- !file.exists(genomesInfo$destination)
spList <- split(genomesInfo[f], 1:sum(f))


# we download the compressed genomes iteratively
for (sp in spList) {
    dir.create(sp$folder)
    download.file(sp$url, sp$temp)

    # we move each genome to destination upon completion
    file.rename(sp$temp, sp$destination)
}
