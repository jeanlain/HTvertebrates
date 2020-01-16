

## %######################################################%##
#                                                          #
####     performs blastn search between TEs in many     ####
####         pairs of species of a given batch          ####
#                                                          #
## %######################################################%##

# this script is launched in stage 04 and used the table of pairs of species 
# whose TE species are blasted against each other generated in this step


library(data.table)
library(parallel)
args <- commandArgs(trailingOnly = TRUE)

# the batch number
job <- as.integer(unlist(strsplit(args[1], ",")))

print(paste("job:", args[1]))

# number of searches to launch in parallel
nCPUs <- as.integer(args[2])


# the pairs of species
sel <- fread("pairsToBlastn.txt")

# the output file names of jobs that are already done
done <- unique(list.files("TEs/blastn/done"))

# this is to avoid redoing the search, in case the job is relaunched



# selects the species pairs of the current batch to process
sel <- sel[batch %in% job]

# generates output file names
sel[, out := paste("TEs/blastn/out/", sp1, "-on-", sp2, ".blastn.out", sep = "")]

# discard species pairs if output files are present
sel <- sel[!basename(out) %chin% done, ]

# query and db file names
sel[, fas := paste("TEs/copies/", sp1, ".TEs.fasta.gz", sep = "")]
sel[, db := paste("TEs/blastn/db/", sp2, ".TEs.fasta", sep = "")]
sel <- sel[!out %chin% done, .(fas, db, out, size)]

# we will launch the larger blast first
setorder(sel, -size)
sel[, job := rep(1:(.N / nCPUs), nCPUs, length.out = .N)]
setorder(sel, job)
nt <- max(1, as.integer(nCPUs / nrow(sel)))

# our function to launch a blast search
blastn <- function(fas, db, out) {
    system(
        paste(
            "gunzip -c",
            fas,
            "| blastn -query - -db",
            db,
            "-max_target_seqs 1 -outfmt 6 -num_threads",
            nt,
            # we do not retain HSP with pID < 75 and score < 200 and alignement shorter than 300 pb.
            # Filtering these hits from the beginning was required, else we might have run out of disk space
            "| awk '{if ($3>=75 && $4>=300 && $12>=200) print $0}' >",
            out
        )
    )

    file.rename(out, gsub("out/", "done/", out, fixed = T))


    # compresses the output as they take a lot of space
    system(paste("gzip", out))

    return(NULL)
}

# we launch the search for the current batch of species pairs with the required number of CPUs ---------------------------------------------
sel[, mcMap(
    f = blastn,
    fas = fas,
    db = db,
    out = out,
    mc.cores = nCPUs,
    mc.preschedule = F
)]

print("finished")
