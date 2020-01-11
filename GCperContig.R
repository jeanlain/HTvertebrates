##%######################################################%##
#                                                          #
####        this scripts counts GC and AT bases         ####
####             per sequences (contigs) of             ####
####          different fasta files (genomes)           ####
#                                                          #
##%######################################################%##

# the only argument is the number of CPUs to use
args <- commandArgs(trailingOnly = TRUE)
nCPUs <- as.integer(args[1])

source("HTvFunctions.R")

# we list the genome sequences downloaded in the first step
genomes <- list.files(pattern = ".gz", full.names = T, recursive = T)

# we will process bigger genomes first
genomes <- genomes[order(-file.size(genomes))]

# the fonction that counts GC and AT bases
GCcontent <- function(file) {
  genome <- readDNAStringSet(file)
  
  # the workhorse command
  counts <- letterFrequency(genome, c("GC", "AT"))

  # we put counts in a data table, one row per sequence (contig)
  dt <- data.table(
      contig = splitToColumns(names(genome), " ", 1), 
      counts, 
      sp = extractSpeciesNames(file)
      )
  
  setnames(dt, 2:3, c("GC", "AT"))
  
  cat(".") #progress indicator
  
  dt
}

# we apply the function in parallel for the different genomes
res <- mclapply(genomes, GCcontent, mc.cores = nCPUs, mc.preschedule = F)
res <- rbindlist(res)

writeT(res, "genomes/GCcontents.txt")

