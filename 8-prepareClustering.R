

## %######################################################%##
#                                                          #
####       This stage prepares the clustering of        ####
####            TE hits to count HTT events.            ####
#                                                          #
## %######################################################%##


# In this step, we estimate the sequence identiy between all TE copies within
# all super families with blastn all vs all (referred to as "self blast"
# afterwards). This will be used to determine of two hits represent the same HTT
# (identity at the hits is lower than identity between copies in each clade,
# refered to as "criterion 1" in the paper), and also for the analysis of TE
# evolution within genomes

# this script uses
# - the TE-TE hits tabular file we retained in the previous step
# - the fasta file fo TE copies in the hits, from the previous step

# the output is the blast hits (tabular format) of all TE copies within each super family

source("HTvFunctions.R")

# where output file will go
dir.create("TEs/clustering")

# selected TE-TE hits from the previous step
httHits <- fread("occ2000ds05.txt")


# STEP ONE, we reduce the number of hits ---------------------------------------------------------
# there are too many hits for this clutering, we only retain 200 of hits
# per single-linkage cluster, as hits from the same cluster are very likely
# to represent the same HTT

# we favor hits that involve the longest protein coding regions
setorder(httHits, -length.aa)
httHits <- httHits[occ <= 200L, ]

# and write them to file
writeT(httHits, "occ200ds05.txt")

# STEP TWO, we extract TE copies sequences in the retained hits ----------------------------------

# we add species and classification information to copy names, as specified in the fasta file of TE copies, so we can extract them
copies <- httHits[, union(
    x = stri_c(sp1, ":", query, ":", f1, "#", superF),
    y = stri_c(sp2, ":", subject, ":", f2, "#", superF)
)]


# writes them to temporary file used by seqtk
write(copies, "temp.bed")

# and we extract these sequences into a new fasta file
seqtk(
    fas = "TEKs/selectedCopiesAA300cl2000.fas",
    bed = "temp",
    out = "TEs/clustering/selectedCopiesKs05occ200.fas"
)
file.remove("temp.bed")

# STEP THREE, we prepare and laucnh the self blastn within each super family of the copies extracted above ------------------------------------------
# we will rename copies with integers, to reduce the size of blastn files, and improve speed in further stages

seqs <- readDNAStringSet("TEs/clustering/selectedCopiesKs05occ200.fas")
copies <- names(seqs)

# we extract super family names and replaces slashes with periods as file names will containes super family names
superF <- gsub("/", ".", stri_extract_last(copies, regex = "[^#]+"), fixed = T)

# attributes integer number to copies, which we will used instead of long copy names
copyID <- data.table(copy = copies, id = 1:length(copies))

# we write this correspondence to disk as it will be reused many times afterwards
writeT(copyID, "TEs/clustering/selectedCopiesKs05occ200.IDs.txt")

# we replace copy names with the integers in the sequecnes
names(seqs) <- copyID[chmatch(names(seqs), copy), id]

dir.create("TEs/clustering/blastn/db", recursive = T)

# we split copy sequences by super families, as blastn searches will be done within super families
all.db <- split(seqs, superF)

# we generate output file names
fasNames <- stri_c("TEs/clustering/blastn/db/", names(all.db), ".fas")

m <- Map(writeXStringSet, all.db, fasNames)
  
# we build blast databases
m <- lapply(fasNames, function(x) {
    system(paste("makeblastdb -dbtype nucl -in", x))
})

# we will split query files for better CPU usage (num_threads of blastn isn't very efficient so we will leave it to 1)
# and for better management of jobs + outputs (which are very big)

# where split query fastas will go
dir.create("TEs/clustering/blastn/split")

# temporary blatn output files with go there
dir.create("TEs/clustering/blastn/out")

# and the output file will go there
dir.create("TEs/clustering/blastn/done")


# the number of sequences per super family, which we used to decide how much we split the queries
nD <- as.numeric(sapply(all.db, length))

# the criterion is the size of query * size of db
nComp <- data.table(
    query = names(all.db),
    nc = nD * as.numeric(sapply(all.db, length))
) 

# we ensure that a query is not split in more than 200 parts
nComp[, nt := round(nc / (max(nc) / 200))]
nComp[nt == 0, nt := 1L]

# we create an integer index corresponding to the sub-query (max 200), for each sequences
indices <- unlist(Map(function(x, y) {
    rep(1:x, length.out = y)
}, nComp$nt, nD))

# this index is used to attribute each sequence to a sub-query (part)
queryParts <- stri_c(
    "TEs/clustering/blastn/split/",
    rep(nComp$query, nD),
    "_",
    indices,
    ".fas"
)


# so we split the copy sequences in these sub queries
splitSequences <- split(unlist(all.db, use.names = F), queryParts)

# and save them to fastas
l <- Map(writeXStringSet, splitSequences, names(splitSequences))

# we create the databases corresponding to the split queries
dbs <- gsub(
    pattern = "split",
    replacement = "db",
    x = stri_extract_last(
        str = names(splitSequences),
        regex = "[A-Z]+.[^_]+"
    )
)


# creates a table listing the blast searches to do
searches <- data.table(query = names(splitSequences), db = stri_c(dbs, ".fas"))
searches[, out := gsub("split", "out", query)]
searches[, out := gsub(".fas", ".out", out)]

# we attribute batches of blastn searches to be launch from a given R instance, on a node (using sarray would probably have been better)
searches[, batch := rep(1:(.N / 30), length.out = .N)]
writeT(searches, "copiesToSelfBlast.txt")

# we now run the script for the self-blastn searches. This is done by batches.
# These use a lot of ram and take quite a bit of time, as we don't limit the number of hits

# example for batch 1
system('sbatch --mail-type=BEGIN,END,FAIL --cpus-per-task=4 --mem=120G --wrap="Rscript TEselfBlastn.R 1 4"')
