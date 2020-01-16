## %######################################################%##
#                                                          #
####            this script processes "raw"             ####
####              blast output betwee TEs               ####
####              of different species to               ####
####          select hits with sufficient pID           ####
#      and removes hits involving artefactual TEs          #
#                   identified at stage 3                   #
#                                                          #
## %######################################################%##

# This script is run at stage 06-filterTEhits.R


library(data.table)
library(parallel)
library(stringi)

# the only argument is the number of CPUs to use
args <- commandArgs(trailingOnly = TRUE)
nCPUs <- as.integer(args[1])


# where the filtered blast output files will go
dir.create("TEs/blastn/filtered/")
pairs <- fread("pairsToFilter.txt")

# output files have been moved and compressed in pairwiseSpeciesBlastn.R
# we therefore change file names accordingly
pairs[, out := gsub("/out/", "/done/", stri_c(out, ".gz"), fixed = T)]


# we will process forward and reverse blast searches (sp1 on sp2, sp2 on sp1) simultaneously,
# so as to select the best hit among the two reciprocal ones

# we therefore extract the forward search
forw <- pairs[sp1 < sp2]

# and the reverse searches
rev <- pairs[sp1 > sp2]

# we swap the two species in the reverse search so that it matches the forward search
setnames(rev, 1:3, c("sp2", "sp1", "rev"))
pairs <- merge(forw, rev)

# paths of filtered output files to be generated
pairs[, filtered := gsub(".gz", "", stri_c("TEs/blastn/filtered/", basename(out)))]


# we import the names of dubious TE families that may not represent true TEs generated in stage 03-findDubiousTEs.R
toRemove <- readLines("TEs/findDubious/familiesToIgnore.txt")


# this function filter the hits for two reciprocal searches:
filterHits <- function(out, rev, filtered, q) {

    # we create an empty data table to avoid a bug in fread() if a blast output file is empty (no hit)
    blast <- data.table()

    # to report the number of hits
    nr <- 0
    if (file.size(out) > 0) {
        # imports "forward" blast search between two species
        # we do not name columns as they wont be retained in the output
        blast <- fread(
            paste("gunzip -c", out),
            sep = "\t",
            header = F,
            drop = c(5, 6, 11)
        )
        nr <- nr + nrow(blast)

        # removes hits of insufficient pID
        blast <- blast[V3 > q, ]
    }

    # we do the same for the "reverse" search
    reverse <- data.table()
    if (file.size(rev) > 0) {
        reverse <- fread(
            paste("gunzip -c", rev),
            sep = "\t",
            header = F,
            drop = c(5, 6, 11)
        )
        nr <- nr + nrow(reverse)
        reverse <- reverse[V3 > q, ]

        # reversing query and subject fields so we can concatenate forward and reverse hits, and remove reciprocal hits (amongst other things)
        reverse[, c("V1", "V2", "V7", "V8", "V9", "V10") := .(V2, V1, V9, V10, V7, V8)]
    }

    # if there is no hit, we exit
    if (nrow(blast) == 0 & nrow(reverse) == 0) {
        return(NULL)
    }
    
    # we stack the forward and reverse hits
    blast <- rbind(blast, reverse)
    rm(reverse)

    # we put the best hits on top (column 12 is the score)
    setorder(blast, -V12)

    # we retain the best hsp per pair of copies, which also removes reciprocal hits
    blast <- blast[!duplicated(data.table(V1, V2))]

    # now we remove dubious TE families and create columns that will be used through our pipeline------------------
    
    # splits copy names into matrices of several columns (fields are separated by colons)
    mat1 <- stri_split(blast$V1, fixed = ":", simplify = T)
    mat2 <- stri_split(blast$V2, fixed = ":", simplify = T)
    f1 <- stri_c(mat1[, 1], "-", mat1[, 6])

    # creates TE family names from species names (first field) and family names
    # so we can filter-out hits involving dubious TEs families
    f2 <- stri_c(mat2[, 1], "-", mat2[, 6])
    f <- !f1 %chin% toRemove & !f2 %chin% toRemove
   
     if (sum(f) == 0) {
        # if there is no retained hit, we exit
        return(NULL)
     }
       
    # we use rbind() makes sure that these will be matrices and not vectors, 
    # even if there is only one hit (to avoid a dimension problem)
    mat1 <- rbind(mat1[f, ])
    mat2 <- rbind(mat2[f, ])

    # we separate TE family name form superfamily name
    subf1 <- stri_split(mat1[, 6], fixed = "#", simplify = T)
    subf2 <- stri_split(mat2[, 6], fixed = "#", simplify = T)

    # we create  new copy names (not containing species names, as these will be in a separate column)
    query <- stri_c(mat1[, 2], mat1[, 3], mat1[, 4], mat1[, 5], sep = ":")
    subject <- stri_c(mat2[, 2], mat2[, 3], mat2[, 4], mat2[, 5], sep = ":")
    
    nr2 <- nrow(blast)
    
    #and we create the final file of filter hits
    blast <- data.table(query, subject, mat1[, 1], mat2[, 1], subf1, subf2, blast[f, -c("V1", "V2")])
    
    # we write file of filtered hits (not returned by the function, as we are limited in the memory we can use here)
    fwrite(
        x = blast,
        file = filtered,
        col.names = F,
        row.names = F,
        quote = F,
        sep = "\t"

    )

    # instead we return the initial number of hits and number of retained hits
    c(basename(filtered), nr, nr2, nrow(blast))
}

# applies the function with the requested number of CPUs
hitNumbers <- pairs[!file.exists(filtered), mcMap(
    f = filterHits,
    out,
    rev,
    filtered,
    q,
    mc.cores = nCPUs,
    mc.preschedule = F
)]

# we stack the reports into a single table
hitNumbers <- data.table(do.call(rbind, hitNumbers))

# and write them to file
fwrite(
    x = hitNumbers,
    file = "filteredStats.txt",
    col.names = F,
    row.names = F,
    quote = F,
    sep = "\t"
)

# we concatenate the output files of filtered hits
system("cat TEs/blastn/filtered/*.out > all.quantile005score200.blastn.out")
