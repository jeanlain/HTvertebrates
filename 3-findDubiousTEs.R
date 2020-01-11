## %######################################################%##
#                                                          #
####          finds TE consensuses (families)           ####
####            that may cover non-TE genes             ####
# 				to exclude them as aterfactual			   #
#                                                          #
## %######################################################%##



# We perform a blastx search of all TE consensus sequences generated
# by RepeatModeler against the non-redundant “nr” protein database of
# NCBI, using Diamond
# We do a similar search against a database of known TE proteins. 
# We discard a TE consensus (hence family) if (i) it has homology to an nr
# protein over at least 90 amino acids in a region that did not show homology
# with a RepBase protein and if (ii) this nr protein does not show a homology of
# at least 35% over >=100 amino acids with a RepBase protein. Homology between
# proteins is determined by a Diamond blastp search of nr proteins fulfilling
# criterion (i) against TE proteins. For all searches, we ignore alignments of
# e-value < 10−3.

# the script uses
# - the fasta files of repeatmodeler family consensuses
# - the repeat masker database of repeated proteins

# the output is a text file listing names of family consensuses to exclude

source("HTvFunctions.R")


# folder where intermediate and results file will be save
workFolder <- ("TEs/findDubious/")
dir.create(workFolder)

# STEP ONE, we blastx consensus of repeat elements against the ncbi non-redundant database of proteins (with diamond) -----------------------------------

# we thus collect all the species TE consensuses, and put them in a single fasta file
# (we prefer doing a single blast for all species, and we need to add species names to the sequences to keep track of them)

consensusFiles <- list.files(
    pattern = "families.fa$",
    full.names = T,
    recursive = T
)

# we determine the species for each consensus fasta. As for the previous step, species names must be present in file paths
species <- extractSpeciesNames(consensusFiles)

# this function  imports a single fasta file
import <- function(file, sp) {
    seqs <- readDNAStringSet(filepath = file)

    # we add the species name to each sequence name
    setNames(seqs, stri_c(sp, "-", names(seqs)))
}

# we import fasta sequence in parallel with 10 CPUs
seqs <- mcMap(
    f = import,
    file = consensusFiles,
    sp = species,
    mc.cores = 10,
    mc.preschedule = F
)


# we combine the sequences into a single DNAStringSet
seqs <- do.call(c, seqs)

setwd(workFolder)

# write all consensus sequences to a single fasta
writeXStringSet(seqs, "allTEconsensuses.fas")


# we obtain the nr database (not recommended if you already have it somewhere)
download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz", "nr.gz")

# we make a diamond database fro it
system("diamond makedb --in nr.gz -d nr -p 10")

# file name of the output of the blastx search
out <- "consensusesOn_nr.out"

# column names of the output (used throughout)
fields <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle")

# we launch the diamond blastx search:
system(
    command = paste(
        "diamond blastx -p 1 --more-sensitive -q allTEconsensuses.fas -d nr.dmnd --quiet -k 20 -f 6",
        paste(fields, collapse = " "),
        "-o",
        out
    )
)


# processing the hits -----------------------------------------------------
blast_nr <- fread(
    input = out,
    header = F,
    sep = "\t",
    col.names = fields
)


# we replace subject id by subject title (description)
blast_nr[, c("sseqid", "stitle") := .(stitle, NULL)]

# remove hits of low evalue
blast_nr <- blast_nr[evalue <= 0.001]
setorder(blast_nr, qseqid, qstart, -qend)

# remove hits nested within others (in the query), as these hits are of lower quality (shorter)
blast_nr[, set := regionSets(data.frame(qseqid, qstart, qend))]
blast_nr <- blast_nr[!duplicated(set)]

# we combine relatively adjacent hits (with a lot of margin). This function also renames columns
combined <- combineHits(
    blast = blast_nr,
    maxDist = 300,
    maxOverlap = 10000,
    maxDiff = 300,
    blastX = T
)


# puts best hits on top, per query
setorder(combined, query, -score)

# attributes a number to every hit related to the same query (starting at 1)
combined[, hit := occurrences(query)]




# STEP TWO, we blast TE consensuses against repeat_pep proteins -----------------------------
# we use RepeatPeps.lib = the repeat modeler classification fasta file "RepeatPeps.lib" 
# it is part of the RepeatModeler package, not installed with this script
# this file should be put in current working directory

# we make a diamond database from it
system("diamond makedb --in RepeatPeps.lib -d repeatPeps")
out <- "consensusesOnRepeatPeps.out"

# we launch the diamond search
system(
    paste(
        "diamond blastx -p 1 --more-sensitive -q allTEconsensuses.fas -d repeatPeps.dmnd --quiet -k 20 -f 6",
        paste(fields, collapse = " "),
        "-o",
        out
    )
)

# we process the hits in the same fashion as the previous step: -------------------------------------------------
blast_repBase <- fread(
    input = out,
    header = F,
    sep = "\t",
    col.names = fields
)

blast_repBase <- blast_repBase[evalue <= 0.001]
setorder(blast_repBase, qseqid, qstart, -qend)
blast_repBase[, set := regionSets(data.frame(qseqid, qstart, qend))]
blast_repBase <- blast_repBase[!duplicated(set)]

# replaces subject name by this generic term
blast_repBase[, sseqid := "repProtein"]

# here we combine all hsps of each consensus, regardless of the protein 
combinedRep <- combineHits(
    blast = blast_repBase,
    maxDist = 100000,
    maxOverlap = 100000,
    maxDiff = 100000,
    blastX = T
)


# STEP THREE ----------------------------------------------------------------------------------------------------
# to find consensuses that may not be TEs, we compare their hits on repbase proteins to hits on nr.
# But since there may be several hit on NR per consensus, we do it on a per-hit (not per-consensus) basis.
# The first hit per consensus (query) is processed first, then the second, etc.

# the function below confronts the hit on repbase to the hit on nr for each consensus:
dubious <- function(occ) {
    merged <- merge(
        x = combined[hit == occ],
        y = combinedRep,
        by = "query",
        suffixes = c("", ".r"),
        all = T
    )

    # we compute the length of the part that aligns on the nr protein before it aligns on the rep protein
    merged[, leading := pmin(qStart.r, qEnd.r) - pmin(qStart, qEnd)]

    # the length pf the part that aligns on the nr protein after it aligns on the rep protein
    merged[, trailing := pmax(qStart, qEnd) - pmax(qStart.r, qEnd.r)]

    # we retain pairs of hits where there is no alignment on repbase,
    m <- merged[is.na(subject.r) |
        
        # or where there are parts aligning on nr protein not covered by alignment on rep proteins (of at least 90 bp)
        leading > 90 | trailing > 90]
    m
}

# we apply the above fonctions to hits involving every consensus
dubiousConsensusHits <- lapply(unique(combined$hit), dubious)
dubiousConsensusHits <- rbindlist(dubiousConsensusHits)

# we only retain hits on nr proteins that are not annotated as repeat proteins :
dubiousConsensusHits <- dubiousConsensusHits[!grepl(
    pattern = "transpo|retro|reverse|integrase|gag|pol-|pol |rna-dependent|polyprot|mobile|jockey|jerky|setmar|copia|recombinase|crypton|mariner|tcmar|tc1|gypsy|helitron|harbi|piggy|maveri|polinton|academ|ltr|cmc|envelop",
    # the keyword list above has been determined after carefull inspection of protein names in the hits
    x = subject,
    ignore.case = T
)
                                            # and with alignment length >100 
                                            & abs(qStart - qEnd + 1) > 100, ]

# we put accession numbers of nr proteins in this new column
dubiousConsensusHits[, acc := splitToColumns(subject, " ", 1)]




# STEP FOUR, we blastp nr proteins similar to TE consensuses against repeat proteins -----------------------------------

# we put their names in a file for seqtk 
# (we use all proteins, not just those in the dubiousConsensusHits table. We can afford it)
writeLines(unique(blast_nr$sseqid), "hitted_nrProteins.bed")

# we extract the protein sequences
seqtk(nr.gz, bed, "hitted_nrProteins.fas")

# the futre output of the blastp
out <- "nrProtsOnRepeatPeps.out"

# we launch the search
system(
    paste(
        "diamond blastp -p 1 -q hitted_nrProteins.fas --more-sensitive -d repeatPeps.dmnd --quiet -k 20 -f 6",
        paste(fields, collapse = " "),
        "-o",
        out
    )
)

# processing the hits --------------------------------------
blastp <- fread(out,
    header = F,
    sep = "\t",
    col.names = fields
)

blastp <- blastp[evalue < 0.001]

# again, combining HSPs on the same query
combinedP <- combineHits(
    blast = blastp,
    maxDist = 1000,
    maxOverlap = 100000,
    maxDiff = 1000
)
    

# we get accession number of nr proteins that have an acceptable hit with a rebpase proteins.
# We will consider them as legit TE proteins
repProtAccessions <- combinedP[abs(qEnd - qStart) + 1 >= 100 & identity >= 35, 
                               unique(splitToColumns(query, " ", 1))]


# so we can extract consensuses that hit to other proteins than these
dubiousFamilies <- dubiousConsensusHits[!acc %in% repProtAccessions, unique(query)]

# we write these family names to disk. 
write(dubiousFamilies, "familiesToIgnore.txt")

# We still blast these TEs to find HTT.
# This leaves us the possibility to remove these TEs afterwards

