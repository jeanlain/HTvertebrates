## %######################################################%##
#                                                          #
####            this script blasts TE copies            ####
####           in retained hit groups against           ####
####                all TEs of the host                 ####
####        species, to check the reliability of        ####
####            transfers (by the number of             ####
####          copies they may involve), so as           ####
####            to test whether apparent HTT            ####
####       may in fact result from contamination        ####
#                                                          #
## %######################################################%##

# this script is launched at stage 11-hitGroupEvaluation.R

# the script reads fastas of copy sequences generated at stage 11
# and uses the blast databases of TE copies generated at stage 04-blastTEs.R
# the output is a blastn tabular output (hits)

source("HTvFunctions.R")

# the only user input is the number of CPUs to use
args <- commandArgs(trailingOnly = TRUE)
nCPUs <- as.integer(args[1])

# we list the file names of queries to blast
fasFiles <- list.files(
    path = "TEs/clustering/testConta",
    pattern = ".fas",
    full.names = T
)

# we extract species names from file names
sp <- extractSpeciesNames(basename(fasFiles))

# we also list blastn databases (only .nin files) of TEs for all species,
# these are generated in stage 04-blastTEs.R
dbs <- list.files(
    path = "TEs/blastn/db",
    pattern = ".nin",
    full.names = T
)

# we remove the extension name from the file
dbs <- gsub(
    pattern = ".nin",
    replacement = "",
    x = dbs
)

# we extract species names from database names
dbSp <- extractSpeciesNames(dbs)

# so that we can order db files to match query files
dbs <- dbs[match(sp, dbSp)]

dir.create("TEs/clustering/testConta/blastn/done", recursive = T) # were output files will go

# names of future output files during the work
outFiles <- stri_c(
    "TEs/clustering/testConta/blastn/",
    sp,
    ".selectedCopies.blastn.out"
)

# and of final output files
doneFiles <- stri_c(
    "TEs/clustering/testConta/blastn/done/",
    sp,
    ".selectedCopies.blastn.out"
)

f <- !file.exists(doneFiles) # this is to avoid redoing completed blast searches

# the function that blasts TEs
blastAgainstTEs <- function(fasFile, db, out, done) {
    system(paste(
        "blastn -query",
        fasFile,
        "-db",
        db,
        "-outfmt 6 -max_target_seqs 20 -max_hsps 1 -out",
        out
    ))

    # moves result to other folder once finished
    file.rename(out, done)
    NULL
}

# we apply the blast in parallel
m <- mcMap(
    f = blastAgainstTEs,
    fasFile = fasFiles[f],
    db = dbs[f],
    out = outFiles[f],
    done = doneFiles[f],
    mc.cores = nCPUs,
    mc.preschedule = F
)

# we concatenate the output
system("cat TEs/clustering/testConta/blastn/done/*.out > TEs/clustering/testConta/blastn/all.out")

print("finished")
