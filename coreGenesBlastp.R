## %######################################################%##
#                                                           #
####        This scripts launches diamond blastp         ####
####         search of BUSCO proteins pairs of           ####
####             species, the list of which              ####
####     is read from a file ("pairsToBlastp.txt")       ####
#                                                           #
## %######################################################%##


library(data.table)
library(parallel)

# the first argument is the job identifier (from 1 to 10)
args <- commandArgs(trailingOnly = TRUE)
job <- as.integer(unlist(strsplit(args[1], ",")))
print(paste("job:", job))

# the second argument in the number of CPUs to use
nCPUs <- as.integer(args[2])

# the pair of species for which we need to do the searches
pairs <- fread("pairsToBlastp.txt", header = T)
pairs <- pairs[batch %in% job]

# our function that does a search
# we only retain one hit per query
blastp <- function(db, aa, out, daa) {
    system(
        paste(
            "diamond blastp --quiet --sensitive -p 1 -k 1 --max-hsps 1 -d",
            db,
            "-q",
            aa,
            "-a",
            daa
        )
    )
    system(paste("diamond view -a", daa, "-o", out)) # converting output file to tabular format (this version of diamond required it)
    file.remove(daa)
    print(paste(out))
}

# we apply the function in parallel for the different species pairs
res <- pairs[!file.exists(out), mcMap(blastp,
    db,
    aa,
    out,
    daa,
    mc.cores = nCPUs,
    mc.preschedule = F
)]
