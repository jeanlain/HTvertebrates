

## %######################################################%##
#                                                          #
####   This stage applies the 2nd round of clustering   ####
#                                                          #
## %######################################################%##
# We cluster hits from different communities into "hit groups" (see the paper methods for details)


# this script uses
# - the tabular file of selected TE-TE hit ("HTT hits") from step 7-TEKsAndHTTfilter.R
# - the self-blastn output generated at step 8


# ouptut is
# - a tabular file of HTT hits with a new column attributing a hit-group number to each, "ocCD00_comm.txt"
# - a tabular file of statistics for each pair of hit community (much too long to describe in a comment)

source("HTvFunctions.R")
httHits <- fread("ocCD00_comm.txt") # HTT hits with community IDs from the previous step


# STEP ONE, we apply criterion 1 to hits from different communities---------------------------------------------------------------
# this is done for communities of of hits between TEs of the same super samily and MRCA
# Indeed, two hits involving pairs of species that do not have the
# same MRCA cannot results from the same HTT (see paper for details)
# The procedure is similar to what we did in the previous step

# here a "group" of hits to cluster is defined by the TE super family and MRCA
# we therefore add a column indicating the "group" of eac hit

httHits[, group := stri_c(gsub("/", ".", superF, fixed = T), "_", mrca)]

# like in the previous stage, split the hits by group and select the column we need for the clustering
hitList <- split(
    x = httHits[, .(
        hit = 1:.N,
        query = q,
        subject = s,
        pID = as.integer(pID * 1000L), # pID is again converted to integer for speed
        group,
        com
    )],
    f = httHits$group
)

# determines the number of hits per group. No clustering need be done if there is just one community of hits per group
nCom <- sapply(
    X = hitList,
    FUN = function(x) length(unique(x$com))
)

hitList <- hitList[nCom > 1L] # only select the hit of groups where clustering has to be done

# and we sort the the hit list by deceasing number of hits (for better CPU usage)
hitList <- hitList[order(sapply(
    X = hitList,
    FUN = nrow
),
decreasing = T
)]

dir.create("TEs/clustering/round2") # results of this round of clustering will go there

# we save the list of hits that will be clustered in a dedicated script
saveRDS(
    object = hitList,
    file = "TEs/clustering/round2/hitsFor2ndClustering.RDS"
)

system("Rscript iterativeSecondClustering.R 20") # applies criterion 1 with 20 CPUs.

# we collect results. We now evaluate the connexions between hit communities that may represent the same HTTs
# results are statistics for every relevant pair of hit commmunities (See iterativeSecondClustering.R)
stats <- fread("TEs/clustering/round2/all.txt")
stats[, links := tot - allDiff] # the number of hit pairs that pass criterion 1 between two communities




# STEP 2, applying "criterion 2" -------------------------------------------------------------------------------------
## we evaluate whether 2 communities can represent the same HTT, based on
# "citerion 2", which relies on the inferred age of the transfer. To reflect the
# same HTT the transfers (represented by the 2 clusters of hits) should not be
# more recent that both clades and we infer the maximal times of transfers from
# the mean Ks between copies at hits of the same community we compare it do the
# Ks of buscos of corresponding clades so we determine the clades involved. For 2
# clusters of hits between young clades A-B (community 1) and C-D (community 2)
# both coalescing to the same MRCA (mrca of A and B is the same that of C and D),
# 2 clades are involved. One is composed of A-C and the order of species B-D, if
# A and C belong to one of the 2 subclades diverging from the MRCA. Remember that
# we took care to place species of one subclade in column sp1 and species of the
# other subclade in column sp2, for each mrca. (done at step 8-prepareClustering.R)

# for each community, we compute the mean Ks of TE-TE hits
httHits[ks < 0, ks := 0] # for some reasons rare Ks values are very sightly < 0, we fix that
ksStats <- httHits[, .(meanKs = mean(ks)), by = .(com, mrca)]

# we get the core gene Ks threshold we used for to filter hits between species of a given MRCA
# we therefore import BUSCO gene Ks (200 AA alignments and one score per BUSCO gene per pair of clades) generated in step 5-coreGeneKs.R
Ks <- fread("gunzip -c Ks200AAnoRedundancy.txt.gz")

# for each pair of subclade diverging from an MRCA ("clade"), we get the Ks threshold mentionned above
KsThresholds <- Ks[, .(q05 = quantile(Ks, 0.005)), by = clade]

ksStats[, threshold := KsThresholds[match(mrca, clade), q05]] # we add the Ks thresolf as a new column

# we retrieve the 2 clades (A-C and B-D in the explanation above) that are involved in every pair of communities
tree <- read.tree("timetree.nwk")
mrca <- mrca(tree) # the matrix of MRCA for all species of the tree

# We identify the clades involved in the HTT of the communities.
# we first retreive two species for each community of a pair, those involved in the first hit
species <- cbind(
    httHits[match(stats$com1, com), cbind(sp1, sp2)],
    httHits[match(stats$com2, com), cbind(sp1, sp2)]
)

# since all sp1 and sp2 from a hit community are from different "young" clades of <40 My,
# their MRCA define the clades we want (mrca of clades A-B and of clades C-D)
# we add these MRCA to the table
stats[, c("AB", "CD") := .(
    mrca[species[, c(1, 3)]], # AB defines the MRCA of clades A and B
    mrca[species[, c(2, 4)]]
)]

# we add a last row for a dummy clade with Ks 0
KsThresholds <- rbind(
    KsThresholds,
    data.table(clade = 0L, q05 = 0)
)

# the "nomatch" argument below below is set to the last row, hence the q05 value retreived would be 0
# hence, Ks would be 0 for MRCA for which we did not compute Ks
stats[, c("KsAB", "KsCD") := .(
    KsThresholds[match(AB, clade, nomatch = .N), q05],
    KsThresholds[match(CD, clade, nomatch = .N), q05]
)]

# and retreive the mean Ks of the hit communities involved
stats[, c("Ks1", "Ks2") :=
    data.table(
        ksStats[match(com1, com), meanKs],
        ksStats[match(com2, com), meanKs]
    )]

# we are now read to apply criterion 2 by adding a logical column to this effect
stats[, crit2 := (pmin(Ks1, Ks2) >= pmin(KsAB, KsCD) &
    pmax(Ks1, Ks2) >= pmax(KsAB, KsCD))]

# and criterion 1 needs to be passed by more than 5% of hit pairs
# to consider the communities as representing the same HTT.
stats[, crit1 := links / tot > 0.05]
# But we are not done with criterion 1 yet (see next)


# STEP 3 ------------------------------------------------------------------------------------------
# to check whether different clusters of hits could represent the
# retention of different parts of a TE transferred once, we compare
# their protein regions.
# the approach is that used in Peccoud et al. 2017 PNAS

# For this we reuse the blastx results of copies against repeat proteins (done in step 6-filterTEhits.R)
blastx <- fread("TEs/blastx/all.copies.successiveBlastx.out")

copyIDs <- fread("TEs/clustering/selectedCopiesKs05ocCD00.IDs.txt") # we retreive copy integer IDs, to speed up some functions
blastx[, query := copyIDs[chmatch(query, copy), id]]

# to reduce the workload discard hits not those involving copies that we did not retain in the pipeline
blastx <- blastx[!is.na(query)]

# we blast proteins similar to copies against themselves. This will be used to
# determine whether two copies have some homologies, that is, is they have
# homology to protein regions that are themselves homologous. This should be
# quite sentitive.

prot <- blastx[, unique(subject)] # we get the protein that are similar to retained copies
protSeqs <- readAAStringSet("TEs/blastx/RepeatPeps.lib") # import protein database used by repeat modeler (retreive from repeatmasker)
names(protSeqs) <- splitToColumns(names(protSeqs), " ", 1) # we have to discard protein descriptions in protein names

dir.create("TEs/blastp") # results of this step will go there

# we create a fasta of all proteins that are similar to TEs
writeXStringSet(
    x = protSeqs[prot],
    filepath = "TEs/blastp/involvedProtOCC200Ks05.fas"
)

# we make a blastp databse
system("makeblastdb -in TEs/blastp/involvedProtOCC200Ks05.fas -dbtype prot")

# and launch the blastp search
system(paste(
    "blastp -query TEs/blastp/involvedProtOCC200Ks05.fas -db TEs/blastp/involvedProtOCC200Ks05.fas -max_target_seqs",
    length(prot), # we set max_target_seqs as the query size to ensure that all possible hits are reported
    "-outfmt 6 -num_threads 5 -evalue 1e-4 -out TEs/blastp/involvedProtOCC200Ks05.self.out"
))


# for the following, we generate these tables:
# all copies involved in each community, we also need the mrca column in this table
dt <- httHits[, .(copy = unique(c(q, s))), by = .(com, mrca)]
# this table is only need to generata the following one:

# blastx hits for copies in each community. (This takes a bit of time)
protHits <- dt[, blastx[query %in% copy], by = .(com, mrca)]

# we count the number of different communities that have copies
# matching proteins for each MRCA (within which hits were clustered).
comPerMRCA <- protHits[, length(unique(com)), by = mrca]

# if there is just one such community, there is nothing to do so we can ignore these
protHits <- protHits[mrca %in% comPerMRCA[region1 > 1L, mrca]]

# imports results of the self blastp launched above
blastp <- fread(
    input = "TEs/blastp/involvedProtOCC200Ks05.self.out",
    header = F,
    col.names = c(
        "query", "subject", "pID", "length", "mismatches", "gapOpen",
        "qStart", "qEnd", "sStart", "sEnd", "evalue", "score"
    )
)

# since we did a basted protein against themselves, we may remove reciprocal hits.
# We will retain the ones with best score (reciprocal hits should have the same score, but we never know)
# we make every hit reciprocal (forward and reverse, this will be useful later).
# In theory, all hit should have a reciprocla one since we did a self-blast.
# but some reciprocal hits are not found for some reasons, and some may have
# different score than "forward" hits

# we select the best HSP for each subject-query pair (including reciprocal hits)
blastp <- removeReciprocal(blastp, removeSelf = F)

# we make an equivalent table with query and subject reversed
blastpR <- copy(blastp[query != subject])
blastpR[, c("query", "subject", "qStart", "qEnd", "sStart", "sEnd") :=
    .(subject, query, sStart, sEnd, qStart, qEnd)]

# we rbind this table with the original hits
blastp <- rbind(blastp, blastpR)
rm(blastpR)

# we creae a unique identifier for the query-subject pair
blastp[, pair := paste(query, subject)]

# for a community, copies may have several protein-coding regions that are
# overlapping or adjacent on a given protein. We can combine these regions, so we
# get larger, non redundant, regions that should more likely represent longer
# ancestral elements

# we combine protein regions that are distant by 10 amino
# acids or less on the same protein in a given community
protRegions <- combineRegions(protHits[, .(
    stri_c(com, subject, sep = " "), # hence, we paste the community number to the protein name
    sStart, sEnd
)],
distance = 10L
)

protRegions <- combineHomologous(protRegions, blastp, protSeqs) # see function for details
writeT(
    data = protRegions,
    path = "TEs/clustering/protRegionsInCommunities.txt"
)


# now comparing communities of TEs in respect at protein regions ----------------------------------------------------
setorder(protRegions, commID)

# we will generate every possible pair of rows from this table in order to
# compare protein regions with each others and assess how they overlap
# protein regions are identified by row indices in the protRegion table
nr <- nrow(protRegions)
rows <- 1:(nr - 1L)

# we generate all possible pairs of protein regions, which requires quite a bit of RAM
pairs <- data.table(
    region1 = rep(rows, nr - rows),
    region2 = unlist(lapply(
        X = 2:nr,
        FUN = function(x) x:nr
    ))
)

# we retreive the communities of these protein regions
pairs[, c("com1", "com2") := data.table(
    protRegions[region1, commID],
    protRegions[region2, commID]
)]

# we will retain only the region pairs for community pairs that are relevant
# we this create an identifier for community pairs
pairs[, comPair := com1 * 10^5 + com2]

# and similar identifiers for our table containing statistics on community pairs.
# Note that com1 < com2 is both cases
stats[, comPair := com1 * 10^5 + com2]

# and remove pairs of regions we no longer need
pairs <- pairs[comPair %in% stats$comPair]
gc() # reclaims some memory

# we create a table containing the coordinate of the protein
# regions between which we want to determine the homology
regionComp <- pairs[, data.table(
    comPair,
    protRegions[region1, .(prot, start, end)],
    protRegions[region2, .(prot2 = prot, start2 = start, end2 = end)]
)]

# and we determine the homology at each pair of protein regions
regionComp <- proteinOverlap(
    regionComp,
    blastp,
    protSeqs
)

# for each pair of communities, we obtain statistics about protein homology
comPairStats <- regionComp[, data.table(
    maxOverlap = max(interAll), # length of the longest hsp between protein regions of the 2 communities
    sumOverlap = sum(interAll), # the total length of hsp between protein regions of the 2 communities
    maxID = max(pID, na.rm = T)
), # the max pID of the HSP
by = comPair
]

# we merge these statistics with our previous statistics on community pairs
mStats <- merge(stats, comPairStats, by = "comPair", all = T)

# if a community was not even present in the protRegions, the values are NAs.
# We put zeros (we leave NAs for com1 == com2 as there was nothing to compare)
mStats[
    com1 != com2 & is.na(maxOverlap),
    c("maxOverlap", "sumOverlap", "maxID") := 0L
]


# we also obtain useful information per community regarding the proteins they cover
# this first requires retreiving total protein lengths from the sequences
protRegions[, protLength := nchar(protSeqs[prot])]

comStats <- protRegions[, .(
    nProt = length(unique(prot)), # number of different proteins (remaining after combineHomologous()),
    protLength = max(protLength)
), # length of longest protein,
longestRegion = max(end - start + 1), # length of longest protein region that is covered by TE copies
by = commID
]

# we add these statistic to our table of community pairs
mStats[, c(
    "nProt1", "longestRegion1", "protLength1", # these columns will apply to com1
    "nProt2", "longestRegion2", "protLength2"
) # same for com2
:= data.table(
        comStats[match(com1, commID), .(nProt, longestRegion, protLength)],
        comStats[match(com2, commID), .(nProt, longestRegion, protLength)]
    )]

# we replace NAs by zeros
mStats[is.na(nProt1), c("nProt1", "longestRegion1", "protLength1") := 0L]
mStats[is.na(nProt2), c("nProt2", "longestRegion2", "protLength2") := 0L]

# we now apply our criterion
mStats[
    maxOverlap < 100 & # if the longest aligned protein part ("maxOverlap") between copies of the 2 transfers is less than 100 aa.
        valid == 0, # and if here is no nucleotide identity between copies within clade (i.e, no self blastn hit),
    crit1 := T
] # criterion 1 is passed

# this leaves the possility that communities with no nucleotide homology and insufficient protein homology
# be different parts of the same TEs that were retained by different clades after a single HTT.



# STAGE THREE, we apply complete linkage clustering to cluster communities in "hit groups"----------------------------------------------------
# the approach is explained in Peccoud et al. 2017 PNAS
# we link communities that pass criteria 1 and 2 as being possibly the result of
# the same HTT. An HTT can only include communities that are all linked with each
# other (= complete-linkage clustering). Because linkage doesn't mean 2
# communities MUST reflect the same HTT, it could just mean a lack of information
# to distinguish HTTs (while two communities not passing the criteria cannot
# possibly be included in the same HTT). So we have to take decisions as to which
# group to from. Within each clade pair and super family, we want to group first:
# communities that cover very similar (or the same) proteins (high maxID)

# we thus place these pairs of communities on top, since we will regroup from top to bottom
mStats <- mStats[order(-maxID)]

# and apply our complete-linkage clustering
hitGroups <- mStats[, cLinkFromPairs(
    V1 = com1,
    V2 = com2,
    linked = crit2 & crit1
)]

# and we add the hitGroup indentfiers (integers) to our community pair table
mStats[, c("hitGroup1", "hitGroup2") :=
    data.table(
        hitGroups[as.character(com1)],
        hitGroups[as.character(com2)]
    )]

# we check that all pairs of communities belonging to the same group passed both criteria
mStats[hitGroup1 == hitGroup2, all(crit2 & crit1)] # this should return TRUE

# we save these satistics for later use
writeT(mStats, "TEs/clustering/comPairStats.Ks05.OCC200.txt")

# we assign each HTT hit to a hit group (new column added)
httHits[, hitGroup := hitGroups[as.character(com)]]

# any hit that is not in a hit group constitutes its own hit group
# (but these small hit groups will be removed by our contamination filter in the next step)
httHits[is.na(hitGroup), hitGroup := 1:.N + max(hitGroups)]

writeT(
    data = httHits,
    path = "oc200HitGroup.txt"
)
