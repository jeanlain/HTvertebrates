## %######################################################%##
#                                                           #
####        This script computes pair-wise Ka and        ####
####         Ks between species at BUSCO genes           ####
#                                                           #
## %######################################################%##

# this script is run at stage 05-coreGeneKs.R

# the only user input is the number of CPUs to use
args <- commandArgs(trailingOnly = TRUE)
nCPUs <- as.integer(args[1])

source("HTvFunctions.R")

dir.create("Ks/KaKs")

# we import all protein and CDS sequences at once
allAA <- readAAStringSet(list.files("CDS", pattern = ".aa.fas", full.names = T))
allCDS <- readDNAStringSet(list.files("CDS", pattern = ".CDS.fas", full.names = T))

# we list results files of diamond blastp searches (tabular)
blastpFiles <- list.files(
    "Ks/blastp",
    pattern = ".out",
    recursive = T,
    full.names = T
) 

# to determined best reciprocal hits, we need to identify reciprocal searches
# we first determine the species we have from sequence names 
sp <- unique(extractSpeciesNames(names(allAA)))

# we generate all possible pairs of species, no reciprocity
pairs <- allPairs(sp, reciprocal = F)

# we generate all possible forward searches (output file names)
allForward <- stri_c("Ks/blastp/", pairs$V1, "/", pairs$V2, ".on.", pairs$V1, ".out")

# all reverse searches
allReverse <- stri_c("Ks/blastp/", pairs$V2, "/", pairs$V1, ".on.", pairs$V2, ".out")

# and we make a table of the output files we actually have for these searches
pairs <- data.table(
    forward = intersect(allForward, blastpFiles),
    reverse = intersect(allReverse, blastpFiles)
) 

# a batch will process approximately 50 species pairs. Bigger batches may cause memory issues
pairs[, batch := rep(1:(.N / 50), length.out = .N)]

pairs <- split(pairs, pairs$batch)

# our function to compute Ka and Ks of pairs of genes based on reciprocal blastp results (one species vs another)
compute_Ks <- function(pairs) {
    
    require(seqinr)
    
    # imports forward blastp search
    
    # the fields present in the blast outputs
    fields <- c("query", "subject", "identity", "length", "mismatches", "indels",
                "qStart", "qEnd",  "sStart","sEnd",  "eValue", "score"  )
    
    forward <- lapply(pairs$forward, function(x) {
        suppressWarnings(fread(x, header = F, col.names = fields, sep = "\t"))
    })
    
    reverse <- lapply(pairs$reverse, function(x) {
        suppressWarnings(fread(x, header = F, col.names = fields, sep = "\t"))
    })
    
    forward <- rbindlist(forward)
    reverse <- rbindlist(reverse)
    
    # we apply the reciprocal best hit criterion to retain hits between putative orthologs
    # for this, we create query-subject identifiers
    forwardHits <- forward[, stri_join(query, subject, sep = "-")]
    reverseHits <- reverse[, stri_join(subject, query, sep = "-")]
    
    # we now retain the reciprocal best hits
    RBH <- forward[forwardHits %in% reverseHits, ]
    batch <- pairs$batch[1]
    
    # retain only the best alignment (HSP) per hit (because there can be overlapping alignments)
    setorder(RBH, -score)
    RBH <- RBH[!duplicated(data.table(query, subject)), ]
    
    # because there appears to be a bug with diamond where rare hits have coordinates longer than the query or subject
    RBH <- RBH[qEnd <= nchar(allAA[query]) & sEnd <= nchar(allAA[subject])]
   
     print(paste(nrow(RBH), "RBH selected for batch", batch))
    
    #returns empty results if there are no hits (almost impossible)
    if (nrow(RBH) == 0) return (data.table()) 
    
    # we extract protein regions involved in hits. 
    # If the alignment doesn't work (it appears that Biostring sometimes has issues with aligning AAStringSets), 
    # it could be useful to replace subseq() below by stri_sub(), which returns character vectors
    sp1Seqs <- RBH[, subseq(allAA[query], qStart, qEnd)]
    sp2Seqs <- RBH[, subseq(allAA[subject], sStart, sEnd)]
    
    # aligns these regions
    aln <- alignWithEndGaps(sp1Seqs, sp2Seqs)
    print(paste("alignment done for batch", batch))
    
    # extracts CDS regions corresponding to aligned protein regions
    sp1Nuc <- RBH[, subseq(allCDS[query], qStart * 3 - 2, qEnd * 3)]
    sp2Nuc <- RBH[, subseq(allCDS[subject], sStart * 3 - 2, sEnd * 3)]
    
    # converts protein alignments into nucleotide alignments
    sp1Nuc <- aaToCDS(aln$pattern, sp1Nuc)
    sp2Nuc <- aaToCDS(aln$subject, sp2Nuc)
    
    # makes a list of seqinrAlignment objects for Ka/Ks computation
    seqinrAlns <- apply(cbind(sp1Nuc, sp2Nuc), 1, seqinrAlignment)
    KaKs <- lapply(seqinrAlns, kaks)
    
    # we stack the results into a data.table
    KaKs <- rbindlist(lapply(
        X = KaKs,
        FUN = as.data.table
    ))
    
    res <- data.table(RBH[, c(1, 2, 7:10), with = F], 
                      KaKs[,.(Ks = ks, Ka = ka)],
                      alnLength = nchar(gsub("-", "", aln$pattern)))
    
    # we write results for this batch, for safety
    writeT(res, stri_join("Ks/KaKs/batch", batch, ".txt"), col.names = F)
    
    res
}


# we process batches in parallel
res <- mclapply(pairs,
    FUN = compute_Ks,
    mc.cores = nCPUs,
    mc.preschedule = F
)

writeT(rbindlist(res), "Ks/KaKs/all.KaKs.txt")

