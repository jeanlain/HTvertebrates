## %######################################################%##
#                                                          #
####  This script launches blast searches of TE copies  ####
####     of a given super family against themselves     ####
#                                                          #
## %######################################################%##

# the only user input is the job identifier (see 8-prepareClustering.R) and the number of CPUs to use

# the output are tabular blastn output files
# these results will be is used to evaluate "criterion 1" for hit clustering
# and also to investigate TE evolution within genome (stage )

library(data.table)
library(parallel)


args <- commandArgs(trailingOnly = TRUE)

# the first argument represent job identifiers (integers) separate by comas
job <- as.integer(unlist(strsplit(args[1], ",")))

# the second argument is the number of CPUs to use
nCPUs <- as.integer(args[2])

print(paste("job:", job))

# we import the table we need to determine the blast searches to launch
searches <- fread("copiesToSelfBlast.txt")

# in case jobs need to be relaunched, we list output files to avoid redoing already completed searches
done <- list.files("TEs/clustering/blastn/done")
toDo <- searches[, batch %in% job & !basename(out) %in% basename(done)]

# we determine the number of threads per blast search (usually 1, but if only a very few searches are done in parallel it can be higher)
nt <- max(1, as.integer(nCPUs / sum(f)))

blastn <- function(query, db, out) {
    done <- gsub("/out/", "/done/", out, fixed = T)
    system(paste(
        "blastn -query", query, "-db", db,
        "-max_target_seqs 100000 -max_hsps 1 -outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore' -num_threads",
        # max_target_seqs is larger than the number of sequences in the db to ensure that all possible hits are returned
        nt,

        # because this is a self blast, we filter results where the query name (an integer) is < subject name,
        "| awk '{if ($4>=100 && $1<$2) print $0}' >",
        # this reduce output file size by two (we only retreive the hit copy1 vs copy2, not the reverse).
        out
    ))


    # moves the output to final destination once finished
    file.rename(out, done)
    return(NULL)
}

# lauchnes the blast search corresponding to the requested job
res <- searches[toDo, mcMap(
    f = blastn,
    query,
    db,
    out,
    mc.cores = min(nCPUs, sum(f)),
    mc.preschedule = F
)]
print("finished")
