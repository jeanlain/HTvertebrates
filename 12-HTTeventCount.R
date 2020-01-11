## %######################################################%##
#                                                           #
####            This stage counts HTT events             ####
####         required to explain the transfers          ####
#                                                           #
## %######################################################%##

source("HTvFunctions.R")

# this script uses
# -the table of HTT hits with transfer information generated at step 10
httHits <- fread("oc200HitGroup.txt")

# - the hit group statistics generated at step 11 (imported later)
# - the "self-blastn" of TE copies against themselves, generated at step 8

# The output is the table of HTT hits with a logical column indicating
# whether each transfer (hit group) can be "explained" by others (see below)

# PRINCIPLE:
# similar TEs may have been brought in two vertebrate clades A and B by several transfers ("hit groups")
# this may give the impression of a direct transfer between those clades
# hence this focal transfer may be "explained" by others
# is that is the case, the copies involved in the focal transfer are similar to other transfers
# and  the species involved in the focal transfer are the same as, or related to, 
# those of the explanatory transfers

# In the pseudocode section of the paper, a main function "requirement_1" evaluates these conditions
# but a naive implementation of this function would be vastly inefficient
# given how frequently this function will be called (thousands of times)
# In practice, most of the repetitive work underlying the pseudocode function
# is done upstream



# STEP ONE, we find TE copies that are similar between transfers --------------------------------------

# For each TE copy from each "focal transfer", we retrieve the highest sequence similarity
# that this copy has with any other copy from every other transfer.
# We do not yet control whether the two transfers involve related species
# Doing so allows changing the criterion for species relatedness independently
# of this first step, should we want to.

# Similarities are only evaluated within TE superfamilies
# since we never looked for homology between copies from different super families.
# So we split the hit table per super family, retaining only the columns we need 
# and replacing slashes with periods to avoid issues with file paths in future files
hitList <- split(
    x = httHits[, .(q, s,             # copy integer IDs
        hitGroup                      # the hitGroup identifiers
    )], 
    f = gsub("/", ".", httHits$superF, fixed = T)
) 

# Note that all hit groups are taken, even those that we considered
# as unreliable in 11-hitGroupEvaluation.R. Doing so allow changing our
# selection of reliable transfers without having to re-run this step.

# we save this list of tables for the script we run below
saveRDS(hitList, file = "TEs/clustering/forHighestSimilarity.RDS")

system("Rscript findHighestSimilarity?R 5")


# We import the results
highestSimilarity <- fread(
  input = "TEs/clustering/networks/all.hitGroupSimilarities.txt",
  col.names = c("copy", "pID", "explanatoryTransfer", "focalTransfer")
  )

# from now on, we will refer to a hit group as a "transfer"
# Using this word helps to understand the concepts we use. 
# This required renaming some columns of the imported table
# in this table, the "copy" belongs to the "focalTransfer".
# Column ("pID") denotes the sequence similarity that the "copy" has with the
# most similar copy of "explanatoryTransfer" 
# (that other copy is not named in the table).
# the explanatoryTransfer may therefore "explain" the focalTransfer.
# But for this, explanatoryTransfer has to fulfil some conditions...




# STEP TWO  ---------------------------------------------------------------------
# we assess whether the explanatoryTransfer involves species 
# that are related to the host species of the "copy" from the focalTransfer

# We determine if the clade involved in the focalTransfer to which the copy belongs 
# is, nested in, or encompasses, either of the clades involved in the explanatoryTransfer. 
# Note that this doesn't require that species are shared between these transfers

# defining the clades requires the species tree
tree <- read.tree("timetree.nwk")

# We first retreive the host species of each copy of the focal transfer,
# we encode species as tip numbers of the timetree
spForCopy <- integer(httHits[,max(q, s)])
spForCopy[httHits[,c(q, s)]] <- httHits[,chmatch(c(sp1, sp2), tree$tip.label)]

# we add the host species in a new column
highestSimilarity[, speciesA := spForCopy[copy]]


# We define two clades involved in each transfer by the mrcas of the species involved
# which are node numbers in the species tree
clades <- httHits[, .(
    cladeA = MRCA(tree, unique(sp1)), # the mrca of the left-clade species (sp1)
    cladeB = MRCA(tree, unique(sp2))  # same for the right-clade species (sp2)
), 
by = .(transfer = hitGroup)
]

# we make integer vectors to quickly retrieve the clades (A or B) of each transfer
cladeA <- cladeB <- integer(max(clades$transfer))
cladeA[clades$transfer] <- clades$cladeA
cladeB[clades$transfer] <- clades$cladeB


# We now attribute each copy involved of the focalTransfer to a clade.
# We start with "query" copies, which belong to clade A
queries <- httHits[,.(copy = unique(q)), by =  hitGroup]
queries[,clade := cladeA[hitGroup]]

# we do the same for "subject" copies, those belonging to clade B
subjects <- httHits[,.(copy = unique(s)), by = hitGroup]
subjects[,clade := cladeB[hitGroup]]

cladeForCopies = rbind(queries, subjects)

# we create a unique integer identifier for each copy in each transfer
cladeForCopies[, copyTransfer := copy * 10000 + hitGroup]

# we do the same for the highestSimilarity table
highestSimilarity[, copyfocalTransfer := copy * 10000 + focalTransfer]

# We now add a colum denoting the clade to which the copy belongs in the focalTransfer
# We call this column "clade_A" to match the name we use in the pseudocode
highestSimilarity[,clade_A := cladeForCopies[match(
  copyfocalTransfer,
  copyTransfer
), clade]]

# we add similar columns indicating the clades involved in the explanatoryTransfer
# we use "cladeC" to match the name used in supplementary figure 4
highestSimilarity[, c("clade_C1", "clade_C2") := .(
  cladeA[explanatoryTransfer],
  cladeB[explanatoryTransfer]
)]


#----------
# we create a matrix that is TRUE if a clade (column/row index) is nested in,
# or includes, another clade (row/column index). The diagonal is also TRUE
# see nestedClades() function in HTvFunctions.R
nested <- nestedClades(tree)

# we add a logical column that tells if cladeA is
# nested, or includes, either clade of explanatoryTransfer
highestSimilarity[, nestedOrIncludes := nested[cbind(clade_A, clade_C1)] | 
                        nested[cbind(clade_A, clade_C2)]]




# STEP THREE ------------------------------------------------------------------------
# We prepare our criterion that imposes that the best identity that a copy 
# of the focalTransfer has with any copy of the explanatoryTransfer is higher than
# the lowest identity the copy has within the focalTransfer

# for this, we obtain the min pID that each copy has within each transfer
minCopyIDs <- rbind(
    httHits[, .(lowestSimilarity = min(pID)), by = .(hitGroup, copy = q)], # for "query" copies
    httHits[, .(lowestSimilarity = min(pID)), by = .(hitGroup, copy = s)]  # and for "subject" copies
) 

# as we did previously, we use an integer to identify the copy-transfer pair
minCopyIDs[, copyTransfer := copy * 10000 + hitGroup]

# we can now add the column denoting the lowest_similarity
# we convert this pID to integer, to compare it 
# to the pID column
highestSimilarity[, lowestSimilarity := minCopyIDs[
    match(copyfocalTransfer, copyTransfer),
    as.integer(lowestSimilarity * 1000)
]] 



# we now select pairs of homologous copies fulfilling our conditions -----------------------------------

# We will also discard the transfers that were considered unreliable in the previous step.
# Doing it only now allows changing the retained transfers 
# without having to re-run the whole script

# We thus import the statistics we generated in the previous script
hitGroupStats <- fread("TEs/clustering/hitGroupStats.txt")
retainedTransfers <- hitGroupStats[retained == T, hitGroup]

selectedCopies <- highestSimilarity[, 
                                   
    # the best identity it has with a copy of the explanatory transfer
    # must be higher than the the lowest similarity it has within the focal transfer
    pID > lowestSimilarity & 
  
    # its host clade must be nested in, or encompass, either clade of
    # the explanatory transfer
    nestedOrIncludes == TRUE &

    # the transfers must be among those considered as reliable
    explanatoryTransfer %in% retainedTransfers & focalTransfer %in% retainedTransfers
    ]


# STEP FIVE ----------------------------------------------------------------------------
# we evaluate whether explanatory transfers could have brought TE copies
# in all the species composing the focalTransfer


# We make several objects (lists) to optimise the procedure.
# We first list the species that carry copies that are
# similar to those of other transfers harbouring related species

# at this stage, we no longer care about TE copies, 
# so we simplify the table with unique()
explainedSpecies <- unique(highestSimilarity[selectedCopies, .(
  focalTransfer, 
  explanatoryTransfer, 
  speciesA
  )])

# we split speciesA of this table by explanatoryTransfer then by focalTransfer
# we use a modification of the split() function that adds recursiveness
explainedSpecies <- Split(
    x = explainedSpecies$speciesA,
    f = explainedSpecies[, list(explanatoryTransfer, focalTransfer)],
    drop = T,
    recursive = T
)

# explainedSpecies[[x]][[y]] returns the species that carry the TEs 
# involved in focal transfer "y" that may have been brought 
# by explanatory transfer "x"
# names(explainedSpecies[[x]]) gives the ids of transfers that are  
# partly explained by transfer x

explainedSpecies <- reList(explainedSpecies, max(retainedTransfers))
# now x can be an integer, which speeds-up access to first-level elements


# As we will check whether all species of the focal transfers have TEs that
# may have been brought by others transfers, we need to list all species per transfer
# we encode species as tree tip numbers
spForTransfer <- httHits[, unique(chmatch(c(sp1, sp2), tree$tip.label)), 
    by = hitGroup]

# We again make a list for quick access. spForTransfer[[x]] will return all species
# involved in transfer "x", x being an integer
spForTransfer <- reList(split(
    x = spForTransfer$V1,
    f = spForTransfer$hitGroup
))


# These lists were generated to optimise the speed of the function below
# This is the function that tells if a focal transfer can be explained by others
# i.e., if "requirement 1" as defined in the methods and pseudocode is passed

requirement1_passed <- function(transfer) {

    # we retrieve all species whose TE copies may have been brought by explanatory transfers
    explSpecies <- unlist(
      explainedSpecies[explanatoryTransfers], 
      recursive = F)
    
    # we of course only consider the species from the transfer we investigate
    explSpecies <- explSpecies[names(explSpecies) == transfer]
    # this selection is not very efficient, but a faster version of 
    # the function was more complex to understand and distant from the 
    # pseudocode
    
    # we return wether these species constitue all species of the focal transfer
    # and if there are at least 2 contributing transfers
    all(spForTransfer[[as.integer(transfer)]] %in% 
          unlist(explSpecies, use.names = F)) & 
      length(explSpecies) >= 2L
}



# we evaluate requirement 1 on transfers iteratively, sorted by "reliability" score ---------------------------------

# to determine if a transfer can be explained by others, 
# we will inspect the less "reliable" transfers first
# these are considered less likely to represent a 
# "direct" transfer event between clades

# we put the best htt hits on top, for each transfer
setorder(httHits, hitGroup, -pID)

# the reliability score of a transfer is based on 
# the pID of best hits of copies involved (see Methods)
# the sum of the best pIDs for over "query" copies (cladeA)
pIDQ <- httHits[!duplicated(data.table(hitGroup, q)), .(sumID = sum(pID)), by = hitGroup]

# and for the subject copies (cladeB)
pIDS <- httHits[!duplicated(data.table(hitGroup, s)), .(sumID = sum(pID)), by = hitGroup]

# we add the "score" as a new column (the lowest of the two sum of pIDs)
hitGroupStats[, score := pmin(pIDQ$sumID, pIDS$sumID)]

# we extract retained transfers, ordered by reliability score, as an integer vector
orderedTransfers <- hitGroupStats[retained == T, hitGroup[order(score)]]


# we are now ready to iterate over transfers -----------------------
# we initialise two vectors:
# this one will contain the identifiers of transfers that are explained by others
explained <- NULL

# and this one contains the transfers than may explain others
# initially, all transfers are allowed to explain others
# but if one is explained, hence "indirect", it does not
# correspond to a "movement" of TEs, and therefore it cannot 
# explain other transfers
explanatoryTransfers <- retainedTransfers

for (focalTransfer in orderedTransfers) {
    if (requirement1_passed(focalTransfer)) {
        # if the transfer can be explained by others,
        # we may remove it from explanatory transfers
        explanatoryTransfers <- setdiff(explanatoryTransfers, focalTransfer)
      
        # we investigate if transfers that were explained by the focal one
        # can still be explained without it
      
        # we define the transfers to investigate
        toInvestigate <- intersect(
          explained, 
          names(explainedSpecies[[focalTransfer]])
          )
        
        # we check if all these transfers to can 
        # still be explained without the focal transfer
        stillExplained <- sapply(toInvestigate, requirement1_passed)
        
        if (any(!stillExplained)) {
          
          # if any of these transfers can no longer be explained 
          # after removal of the focal transfer we restore the 
          # focal transfer into the list of explanatory transfers
          explanatoryTransfers <- c(explanatoryTransfers, focalTransfer)
          
          cat(".") # progress indicator
          
        } else {
          # else we consider the focal transfer as explained
          explained = c(explained, focalTransfer)
        }
    }
}

# WE ARE NOW FINISHED WITH THE COUNT OF HTT EVENTS -------------------------------------



# we only retain the columns we need in the htt hit table
httHits <- httHits[, .(
  query, subject, q, s, sp1, sp2, f1, f2, 
  superF, mrca, pID, length, qStart, qEnd, 
  sStart, sEnd, ka, ks, length.aa, com, hitGroup
  )]

# we add two logical columns:
# - "retained" is TRUE for transfers we retained (based on the contamination filter , sufficient divergent time and low Ks)
# - "independent" is TRUE is a transfer is not explained by other 
# i.e. a htt event that may be seen as "independent" or "direct"
httHits[, c("retained", "independent") := .(hitGroupStats$retained[hitGroup], ! hitGroup %in% explained)]

# we save this whole table, including transfers that are not retained
writeT(httHits, "HTThitsAssessed.txt") 



# We make supplementary dataset 4 of retained hits --------------------------

# we import the correspondance between repeat modeler super families and more common super family names,
# which we will use from now on. This file is provided with the scripts
corres <- fread(
    input = "superF.txt",
    header = F,
    col.names = c("superFam", "subClass", "newName")
)

# we add more common super family names to the hits
httHits[, superFName := corres[chmatch(superF, superFam), newName]]

# we remove transfers and columns we don't retain
retainedHits <- httHits[retained == T, -c("retained", "mrca", "superF")]

# we replace column names with more user-friendly ones
setnames(
    x = retainedHits,
    old = c("query", "subject", "f1", "f2", "com", "superFName"),
    new = c("copy1", "copy2", "consensus1", "consensus2", "community", "superfamily")
)


writeT(retainedHits, "supplementary-data4-retained_hits.txt")
