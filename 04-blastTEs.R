## %######################################################%##
#                                                          #
####          This stage prepares and launches          ####
####    the pair-wise reciprocal blastn searches        ####
####            to find homologous TE copies            ####
####    between species (one genome against another)    ####
#                                                          #
## %######################################################%##


# this script uses the gzipped fasta of TEs generated at stage 2.
# as well as the timetree of the species
# the output are the tabular blastn output between the TEs of every relevant pair of species

source("HTvFunctions.R")


# STEP ONE, we define the relevant species pairs ------------------------------------------------
# we will not search for HTT between species that are too closely related

# hence, we need to timetree of the species
tree <- read.tree("additional_files/timetree.nwk")

# we obtain the divergence times between all species (as a matrix).
distMat <- cophenetic(tree)

# we turn it into a 3-column data.table
pairs = setNames(data.table(as.table(distMat)), c("sp1","sp2","divTime"))

# we only retain species that have diverged before the last 40 My
selectedSpeciesPairs <- pairs[divTime > 80]

# we will launch the longer blast searches first (those involving the bigger fasta files), to better use the CPUs
selectedSpeciesPairs[, size := file.size(stri_c("TEs/copies/", sp1, ".TEs.fasta.gz")) 
                            + file.size(stri_c("TEs/copies/", sp2, ".TEs.fasta.gz"))]

setorder(selectedSpeciesPairs, -size)

# pairs involving Ambystoma mexicanum (which has much more TEs than the others) will constitute  
# a dedicated batch as we wanted to avoid issue with jobs exceeding the allowed runtime at genotoul
axolotl <- selectedSpeciesPairs[, sp1 == "Ambystoma_mexicanum" | sp2 == "Ambystoma_mexicanum"]

# pairs are assigned to "batches" of blast searches to launch on our cluster
selectedSpeciesPairs[axolotl, batch := 0L]

# we create batches that contains approximately 4000 searches (species pairs). 
# This was determined given the number of CPUs per node on the server (30 cores) and the maximal duration of jobs (4 days).
# Using an array of jobs via sbatch --array would have been better, but we didn't know how to do it by then)
selectedSpeciesPairs[!axolotl, batch := rep(1:(.N / 4000), length.out = .N)]

writeT(selectedSpeciesPairs, "pairsToBlastn.txt")



# STEP TWO, we make blastn database for the TE copies of each species ----------------------------------------------

# we list compressed fasta files of TE copies generated in stage 02-TEextractionAndComposition.R
files <- list.files("TEs/copies", pattern = ".gz", full.names = T)

# where the databases will go
dir.create("TEs/blastn/db", recursive = T)

# where the output of the blast search will go
dir.create("TEs/blastn/done")

# where the output of blast go when the search is not finished
# this helps to determine completion
dir.create("TEs/blastn/out")

db <- gsub("copies", "blastn/db", gsub(".gz", "", files))
f <- !file.exists(paste(db, ".nin", sep = ""))

# this function makes blastn databases from gzipped fastas of TE copies :
makedb <- function(fasta, db) {
    system(
        paste(
            "gunzip -c",
            fasta,
            "| makeblastdb -in - -dbtype nucl -title",
            db,
            "-out",
            db
        ),
        ignore.stderr = T
    )
    return(NULL)
}

# we apply the above function in parallel on 10 CPUs
res <- mcMap( 
    f = makedb,
    fasta = files[f],
    db = db[f],
    mc.cores = 10,
    mc.preschedule = F
)




# STEP THREE, we launch the blastn on the selected species pairs (several jobs sent for different batches) --------------------
# jobs were launched manually 
# this is the command for batch 1 with 20 CPUs
system('sbatch --mail-type=BEGIN,END,FAIL --cpus-per-task=20 --mem=60G --wrap="Rscript pairwiseSpeciesBlastn.R 1 20"')


