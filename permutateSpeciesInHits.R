## %######################################################%##
#                                                          #
####          This scripts shuffle sspecies to          ####
####       generate random species pairs involved       ####
####            in transfers, for each super            ####
####       family, avoiding to generate "illegal        ####
####   pairs" for which no HTT can be infered, due to   ####
####         our requirement that species must          ####
####             be divergent by 2*120 My.              ####
#                                                          #
## %######################################################%##


source("HTvFunctions.R")


args <- commandArgs(trailingOnly = TRUE)

# number of CPUs to use
nCPUs <- as.integer(args[1])

# an second optional argument is a node number to permute only species within this clade. 
# Defaults to the basal node to permute all species
tree <- read.tree("timetree.nwk")
node <- ifelse(is.na(args[2]), which.max(node.depth(tree)), as.integer(args[2]))

# a third optional argument is the number of replicates we want, defaults to 1000
n <- ifelse(is.na(args[3]), 1000L, as.integer(args[3]))

# we compute the number of legal permutation we need per job launched in parallel
maxn <- ceiling(n / nCPUs)




# STEP ONE, we select the hits and delineate the scope of permutations -------------------------------------------

# we import the hit representing htt
retainedHits <- fread("supplementary-data4-retained_hits.txt")

# as several species are involved in a hit group even though this would just represent one transfer, 
# we only retain the two species for the best hit (best pID) per hit group. We also replace their 
# names with tip numbers of the tree, as integer numbers speed things up and save memory considerably

# we put the best hits on top 
setorder(retainedHits, -pID)

# so we extract the one we one (note that we only consider "independent" transfers)
# and we replace species names with tip numbers at the same time
HTThits <- retainedHits[independent == T, .(
    sp1 = chmatch(sp1[1], tree$tip.label),
    sp2 = chmatch(sp2[1], tree$tip.label),
    repl = 0L
    # repl is a replicate number (= 0 for original, non-permuted species)
), by = .(hitgroup, superfamily)]


# we determine the scope of permutations. Ideally, we permute species involved in htt within a TE super family, 
# but some contain too few HTT for this to be meaningful, 
# and others contain too many, which makes impossible to obtain only "legal permutations"
# to make our choices, we compute the number of independent transfers per super family, 
nTr <- HTThits[, .N, by = superfamily]

# we add a column to denote the TE class we pool super families (within TE classes) that are involved in less than 20 transfers
nTr[, class := ifelse(grepl("CMC|hAT|Mariner|Maverick|Merlin|PIF|PiggyBac", superfamily),
    "DNA",
    "RNA"
)]

nTr[, combined := ifelse(test = N < 20, 
                         yes = paste("other", class, sep = " "),  #for pooled superfamilies, we use "other" + the TE class
                         no = superfamily)]

# we replace names of underrepresented super families with the combined names
HTThits[, superfamily2 := nTr[chmatch(HTThits$superfamily, superfamily), combined]]

# we recompute the numbers of transfers with the pooled super families
nTr <- nTr[, .(N = sum(N)), by = combined]

# as there may be too many hits per superfamily to obtain only "legal" permutations, 
# we will split certain super families in batches of hits (when they encompass more than 120 transfers)
# we add a column for the number of time a superfamily has been seen in a transfer
HTThits[, occ := occurrences(superfamily)]

# we compute the max number of transfers we allow per superfamily batch
nTr[, maxi := N / ceiling(N / 121)]

# which we transfer to the hits table
HTThits[, maxi := nTr[match(HTThits$superfamily, combined), maxi]]


# for mariners, we need even smaller batches of ≤ 61 transfers 
# (for some reason, it is harder to obtain legal permutations in these TEs)
HTThits[superfamily == "Tc1/Mariner", maxi := 61L]
HTThits[, batch := ceiling(occ / maxi)]


# we split the hits by superfamily and batch
# a "hit" will be a pair of species associated with a permutation number
hitList <- split(x = HTThits[, .(sp1, sp2, repl)],
                 f = list(HTThits$superfamily, HTThits$batch),
                 drop = T) 

# we get the total number of species
nSpecies <- length(tree$tip.label)

# this matrix will contain the permuted "replacement species" during the work, one set of permuted species per column.
newsp <- matrix(rep(1:nSpecies, 10^5), ncol = 10^5)
# Note that we anticipate 10^5 permutations per batch, although we keep much fewer. 
# This is because many permutations may be "illegal".
# The row indices of this matrix = the original species (integer ids)


# we determine which simulated transfers (those with permuted species) are "legal"-------------------------------
# we create the matrix of divergence times
divTimeMat <- cophenetic(tree)

# this vector will be TRUE for species pairs that are "too close" to be involved in HTT 
# (a species here is a row/column index of the logical matrix)
tooClose <- divTimeMat < 240

# we retrieve the species we will permute (tip numbers of the tree, for the clade/node we focus on)
toShuffle <- tipsForNode(tree, node)


# STEP TWO, we permute species for htts (hits) of a batch ----------------------------------------------------------
shuffleSpecies <- function(hits, superfamily) {
    # we print progress, which is the only use of the superfamily argument
    print(superfamily)
    
    # to speeds things up, we generate 10^5 permutations in a row as the vast
    # majority lead to invalid transfers. This allows taking advantage of
    # vectorised functions after that.
    # for that, we need the number of hits
    nHits <- nrow(hits)
 
    # we replicate the hits 10^5 times, but incrementing species ids at each replicate (by the total number of species = l)
    sp1 <- rep(0:(10^5 - 1), each = nHits) * nSpecies + hits$sp1
    sp2 <- rep(0:(10^5 - 1), each = nHits) * nSpecies + hits$sp2
    
    # this function performs the permutations in parallel. job is a simple integer identifier
    shuffleWork <- function(job) {

        # we prepare a matrix of "legal permutations" 
        # it is the same format as the newsp matrix, where row numbers = species ids, 
        # columns are different permutations and values are replacement species
        legalPermutations <- NULL

        # indicator to tell when to stop
        g <- 0
         
        # until we obtain maxn legal permutations:
        repeat {
           
            # this creates a matrix of permuted species, with 10^5 columns (each is a vector of shuffled species)
            # this is the workhorse function, all the rest is result handling
            newsp[toShuffle, ] <- replicate(10^5, sample(toShuffle))

            # we replace original species by the sampled ones in the transfers
            # the left-hand column is the original species, the right-hand one the replacement one
            newPairs <- cbind(newsp[sp1], newsp[sp2])

            # we compute the number of illegal transfers per permutation (we use a matrix to quickly count them via colSums)
            illegal <- colSums(matrix(tooClose[newPairs], nrow = nHits))
            
            if (any(illegal == 0L)) {

                # we extract the columns corresponding to permutations with no illegal transfer and add them to the retained permutations
                legalPermutations <- cbind(legalPermutations, newsp[, illegal == 0])
                g <- ncol(legalPermutations)

                # progress indicator
                cat(".")
            }
            
            # we exit once we have the number of legal permutations we want
            if (g >= maxn) {
                break
            }
        }

    
        # we unrolls the matrix of legal permutations (useful to replace original species with the sampled ones),
        newSp <- as.vector(legalPermutations[, 1:maxn])
        #  the index in this vector will be the original species identifier, and its value is the replacement species.
        # We do not retain more than maxn legal permutations (there may actually be up to 10^5 if there were few hits).

        # we generate a permutation id number
        repl <- rep(1:maxn, each = nHits)
        
        # and we make a table of hits involving the permuted species
        # we replicate the transfers maxn times but incrementing species numbers, 
            
        randomHits <- data.table(
            sp1 = rep(hits$sp1, maxn) + (repl - 1L) * nSpecies,
            sp2 = rep(hits$sp2, maxn) + (repl - 1) * nSpecies,
            repl

        )
        
        # we can now easily replace original species with the sampled ones.
        # This is similar to what we did to create the newPairs matrix, except this time we use a data table
        # we also use the job id to generate final permutation identifiers, which will have to differ between jobs
        randomHits[, c("sp1", "sp2", "repl") := .(newSp[sp1], newSp[sp2], repl + job * maxn)]
        
        randomHits
    }
    
    # we apply the above function to batches of hits in parallel
    randomHits <- mclapply(1:nCPUs - 1L,
        shuffleWork,
        mc.cores = nCPUs,
        mc.preschedule = F
    )
    
    #and we stack the results in a single table
    randomHits = rbindlist(randomHits)
    
    # we replace replication numbers with smaller number that do not exceed the number
    # of permutation that was specified
    randomHits[, repl := match(repl, unique(repl))]
    
    randomHits
}


# we apply the permutations to superfamilies successively
randomHTTs <- Map(shuffleSpecies, hitList, names(hitList))

# we stack randomized hits with read ones (that we can differentiate since their replication number is 0)
randomHTTs = Map(rbind, hitList, randomHTTs)

dir.create("permutations")
saveRDS(randomHTTs, file = stri_c("permutations/allPermutations_Node.", node, ".RDS"))
