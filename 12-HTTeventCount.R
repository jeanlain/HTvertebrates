## %######################################################%##
#                                                          #
####            This stage counts HTT events            ####
####         required to explain the hit groups         ####
#                                                          #
## %######################################################%##

source("HTvFunctions.R")

# this script uses
# -the table of HTT hits generated with hit group information at step 10
httHits <- fread("oc200HitGroup.txt")

# - the hit groups statistics generated at step 11 (imported later)
# - and the "self-blastn" of TE copies against themselves, generated at step 8
# The output is the table of HTT hits with a logical column indicating
# whether each hit group can be "exlpained" by others (see below)

# PRINCIPLE:
# similar TEs may have been brought in two animal clades A and B by several transfers ("hit groups")
# this may give the impression of a direct transfer between those clades
# hence this focal transfer (focal hit group) may be "explained" by others
# is so, this requires that the copies involved in the transfer are similar to other transfers
# and that the species involved in the focal transers are the same as, or related to, those in the other transfers
# (plus other subtleties that are evaluated below)



# STAGE ONE, we find copies that are similar between hit groups -------------------------------------------------
# at this tage, we do not control whether hit groupes involve related species
# this allows more flexibility to change this criterion of relatdness in the later stages

# comparison between copies is only done within TE superfamilies
# since we never look for homology between copies from different super families
# so we split the hit table per super family, retaining only the columns we need (copy ids)

hitList <- split(
    x = httHits[, .(q, s, # we need copy IDs
        pID = as.integer(pID * 1000), # blast pIDs, conveted to integer for speed
        hitGroup
    )], # the hitgroup identifiers
    f = gsub("/", ".", httHits$superF, fixed = T)
) # per super family, replacing slashes with periods
# to avoid issues with file paths in future files

# we save this list of table for the script below
saveRDS(hitList, file = "TEs/clustering/forHitGroupPairs.RDS")

# this script finds homologous TE copies between any 2 hit groups, using 5 CPUs
system("Rscript hitGroupPairHomologousCopies.R 5")


# retreiving results from hitGroupPairHomologousCopies.R.
homologousCopies <- fread("TEs/clustering/networks/all.hitGroupPairSharedCopies.txt")
# in this table, the "focalHitGroup" column is a hit group with a TE "copy" (column) having some degree of
# similarity ("pID") with a (unnammed) copy of "explanatoryHitGroup".
# explanatoryHitGroup may therefore "explain" focalHitGroup.
# But for this, explanatoryHitGroup has to fulfil some conditions...



# STAGE TWO, we assess whether explanatoryHitGroup involves species related to those of the focalHitGroup ---------------------------------------------------------------------
# We assess if the clades of focalHitGroup are the same as, nested in, or
# encompassing, clades of explanatoryHitGroup. Note that this doesn't require
# that species are shared.
# "clade A" is be the clade on the left of table, (hence containing
# the "query" copies in the httHits table) in the hit group hits
# the other clade is "clade B"

# defining clades requires the species tree
tree <- read.tree("timetree.nwk")

# we generate the pairs of clades involved in each hit group (transfer), encoded as node numbers in the species tree
clades <- httHits[, .(
    cladeA = MRCA(tree, unique(sp1)), # the mrca (tree node) of the left-clade species (sp1)
    cladeB = MRCA(tree, unique(sp2))
), # same for the right-clade species
by = hitGroup
]

# this vector will make the correspondance between a hit group (its index) and the mrca of the species involved (its value)
mrcas <- httHits[match(1:max(hitGroup), hitGroup), mrca]

# same principle as above for the 2 mrcas composing the hit groups.
# these vector speed up the command below as this avoid performing a match on the huge homologousCopies table
cladeA <- clades[match(1:max(hitGroup), hitGroup), cladeA]
cladeB <- clades[match(1:max(hitGroup), hitGroup), cladeB]

# we add these elements as new columns
homologousCopies[, c(
    "cladeA", # clade A of the focalHitGroup
    "cladeB", # clade B
    "cladeA_2", "cladeB_2", # equivalent for explanatoryHitGroup
    "mrcaFocal"
) # mrca of the species in focalHitGroup (both clades combined)
:= .(
        cladeA[focalHitGroup], cladeB[focalHitGroup],
        cladeA[explanatoryHitGroup], cladeB[explanatoryHitGroup], mrcas[focalHitGroup]
    )]


#----------
# we create a matrix that is TRUE if a clade (column/row index) is nested in,
# or includes, another clade (row/column index). The diagonal is also TRUE
nested <- nestedClades(tree)
nested2 <- nestedClades(tree, F)

# we add logical colums:
homologousCopies[, c(
    "nestedA", # TRUE if cladeA of focalHitGroup is nested, or indludes, a clade of explanatoryHitGroup
    "nestedB", # same for cladeB of focalHitGroup
    "nestedMRCA"
) # TRUE if the mrca of focalHitGroup is nested within a clade of explanatoryHitGroup
# (in which case we will consider that explanatoryHitGroup cannot exlpain focalHitGroup)
:= .(
        nested[cbind(cladeA, cladeA_2)] | nested[cbind(cladeA, cladeB_2)],
        nested[cbind(cladeB, cladeA_2)] | nested[cbind(cladeB, cladeB_2)],
        nested2[cbind(cladeA_2, mrcaFocal)] | nested2[cbind(cladeB_2, mrcaFocal)]
    )]

# we remove non-longer needed columns (as the table is really big)
homologousCopies[, c("cladeA", "cladeA_2", "cladeB", "cladeB_2", "mrcaFocal") := NULL]

# we also retreive the host species of each copy of hit group 1,
temp <- httHits[, .(
    cop = c(q, s),
    sp = chmatch(c(sp1, sp2), tree$tip.label)
)] # this encodes the species as integer

# we add the host species in a new column
homologousCopies[, sp := temp[match(copy, cop), sp]]
rm(temp)


# STAGE THREE ------------------------------------------------------------------------
# We prepare our criterion that imposes that the best identity that a copy of focalHitGroup has with any
# copy of explanatoryHitGroup is not lower that the lowest identity the copy has within focalHitGroup

# for this we obtain the min pID each copy has within each hit group
minCopyIDs <- rbind(
    httHits[, (minID <- min(pID)), by = .(hitGroup, copy = q)], # for "query" copies
    httHits[, (minID <- min(pID)), by = .(hitGroup, copy = s)]
) # and for "subject" copies

# we add to the homologousCopies table an integer column to identify the copy-hitgroup pair
# this uses the fact that both are integers, and hitGroup is < 10000
minCopyIDs[, copyHitGroup := copy * 10000 + hitGroup]

# we do the same for the homologousCopies table
homologousCopies[, copyfocalHitGroup := copy * 10000 + focalHitGroup]

# now we can add the column denoting the min pID
homologousCopies[, minID := minCopyIDs[
    match(copyfocalHitGroup, copyHitGroup),
    as.integer(V1 * 1000)
]] # we take the opportunity to convert pID to integer (for memory reason)

# we add the following logical columns
homologousCopies[, c(
    "inCladeA", # TRUE if the copy of focalHitGroup belongs to cladeA
    "inCladeB"
) # TRUE if it belongs to cladeB
:= .(
        copyfocalHitGroup %in% httHits[, q * 10000 + hitGroup],
        copyfocalHitGroup %in% httHits[, s * 10000 + hitGroup]
    )]



#---------------------------------------------------------------------------------
# we now determine if focalHitGroup and explanatoryHitGroup were not considered independent
# based on the criterion 1 of hit clustering that separated hit groups (not
# counting the absence of homology at protein regions). If that is the case, we
# will decide that explanatoryHitGroup cannot explain focalHitGroup.
mStats <- fread("TEs/clustering/comPairStats.Ks05.OCC200.txt") # file generated in step 10-hitClusteringRound2.R
crit1 <- mStats[, sum(links) / sum(tot), by = .(hitGroup1, hitGroup2)]
crit1[, pair := ifelse(test = hitGroup1 < hitGroup2,
    yes = hitGroup1 * 10000 + hitGroup2,
    no = hitGroup2 * 10000 + hitGroup1
)]

homologousCopies[, hitGroupPair := ifelse(focalHitGroup < explanatoryHitGroup,
    focalHitGroup * 10000 + explanatoryHitGroup,
    explanatoryHitGroup * 10000 + focalHitGroup
)]
homologousCopies[, indep := hitGroupPair %in% crit1[V1 < 0.05, pair]]
rm(mStats, crit1)



# STAGE FOUR, we order hit groups by "reliability" score (see methods) ---------------------------------------------------------------
# to determine if a hit group can be explained by others, we will inspect the less "reliable" hit groups first
# these are considered less likely to represent a "direct" transfer event between clades A and B

# we put the best htt hits on top, for each hit group
setorder(httHits, hitGroup, -pID)

# the reliability score of a hit group is based on the pID of best hits of copies involved (see Methods)
# the sum of the best pIDs for over "query" copies (cladeA)
pIDQ <- httHits[!duplicated(data.table(hitGroup, q)), .(sumID = sum(pID)), by = hitGroup]

# and for the subject copies (cladeB)
pIDS <- httHits[!duplicated(data.table(hitGroup, s)), .(sumID = sum(pID)), by = hitGroup]

# as we not consider the hit groups that we considered unreliable in the previous step,
# we imports the statistics on hit groups we generated
hitGroupStats <- fread("TEs/clustering/hitGroupStats.txt")

# and adds the "scoer" as a new column (the lowest of the two sum of pIDs)
hitGroupStats[, score := pmin(pIDQ$sumID, pIDS$sumID)]

# we extract retained hitGroups, ordered by reliability score, as an integer vector
orderedHitGroups <- hitGroupStats[keep == T, hitGroup[order(score)]]


# we now select homologous copies fulfilling our conditions ----------------------------------------------------------
# some could have been evaluated earlier, but we prefer applying filters at the last moment


selectedCopies <- homologousCopies[, pID >= minID & # the best identity with it has with a copy of the explanatory hit group
    # must not be lower than its lowest id within the focal hit group
    indep == F &

    # the species of the explanatory hit group must be related to its host species
    ((inCladeA == T & nestedA == T) | (inCladeB == T & nestedB == T)) &

    # the mrca of focalHitGroup must not belong to a clade of explanatoryHitGroup
    nestedMRCA == F &

    # we also ignore any hitgroup that is not among those considered are reliable
    explanatoryHitGroup %in% orderedHitGroups & focalHitGroup %in% orderedHitGroups]



# STAGE FIVE --------------------------------------------------
# we evaluate whether explanatory hitGroups could have brought TE copies
# in all the species composing the focalHitGroup

# to do so,  we list the species that have copies that are
# similar to those of other hit groups harboring related species.

# we make a temporary table of the columns we need for the copies that passed our criteria
temp <- unique(homologousCopies[selectedCopies, .(focalHitGroup, explanatoryHitGroup, sp)])

# we split the species of this table by pair of focal hit group - explanatory hit group
spList <- split(
    x = temp$sp,
    f = temp[, list(focalHitGroup, explanatoryHitGroup)],
    drop = T
)

# we extact the two hit groups (focal and expanatory) from the list's names
hitGroups <- splitToColumns(names(spList), ".")

# we use the explanatory hit groups (column 2) as names of the list
names(spList) <- hitGroups[, 2]

# we re-split the list by focal hitgroup
spList <- split(spList, hitGroups[, 1])

# and we reorder the list so that explainedSpecies[[x]][[y]] below gives the species
# of expanatory hit group "y" that carry TEs which are similar to those carried
# by species of focal hit group "x"
# x is an integer, hence indices of the first-level elements.
# y is still a character, and represent names of the second-level elements
explainedSpecies <- reList(
    spl = spList,
    len = max(orderedHitGroups)
)

# we will also need to know the hit groups that are partly explained by each each group
# this is more or less the reciprocal of the above list, except that we do not include species
# for this we just need every pair of focal hitgroup and explanatory hit group
temp <- unique(temp[, .(focalHitGroup, explanatoryHitGroup)])

# explainedHitgroups[[x]] below returns the focal hit groups whose TEs may have been "explained"
# by expanatory hit group x (because these TEs are similar to those of hit group x)
# x is an integer, hence the indices of the elements of the list
explainedHitgroups <- reList(
    spl = split(
        x = temp$focalHitGroup,
        f = temp$explanatoryHitGroup
    ),
    len = max(orderedHitGroups)
)

# as we will have to check whether all species of the focal hit groups have TEs that
# may have been brought by explanatory hit groups, we need to list all species per hit group
spForHitgroup <- httHits[, unique(chmatch(c(sp1, sp2), tree$tip.label)), # we encode species as tree tip numbers
    by = hitGroup
]

# we again make a list for quick acess. spForHitgroup[[x]] will return all species
# for hitGroup "x", x being an integer
spForHitgroup <- reList(spl = split(
    x = spForHitgroup$V1,
    f = spForHitgroup$hitGroup
))


# our function that returns a logical to tell if a focal hit group (integer number) can be explained
isExplained <- function(hitGroup, explainedSpecies) {

    # we retreive the species whose copies may have been transferred by other hit groups
    explSpecies <- explainedSpecies[[hitGroup]]

    # and we get the other species
    nonExplSpecies <- setdiff(
        x = spForHitgroup[[hitGroup]],
        y = unlist(explSpecies)
    )

    # we return wether all species should have have homologous copies with other hit groups
    # and if there are at least 2 such explanatory hit groups
    length(nonExplSpecies) == 0 & length(explSpecies) >= 2L
}


# if we consider a hit group as "explained" by others, it cannot itself explain others
# so we have another function to remove a hit group from those than may have grought TEs
# into species from focal hit groups
removeHitGroupFrom <- function(hitGroup, explainedSpecies) {
    lapply(
        X = explainedSpecies,
        FUN = function(x) x[!names(x) %chin% as.character(hitGroup)]
    )
}

# we now processes focal hit groups in increasing "reliability" score -----------------------------
# to determine if each can be is explained by others or "required" to exlpain others

# exlpained[x] below will be TRUE if hitgroup x is explained by others
explained <- logical(nrow(hitGroupStats))

# another logical vector that will tell whether a hit group is "required" to explain others
required <- logical(nrow(hitGroupStats))

for (focalHitgroup in orderedHitGroups) {
    if (isExplained(focalHitgroup, explainedSpecies)) {
        # if the hit group can be explained by others,
        # we investigate whether we can remove it from the explanatory hit groups
        # we duplicate the explainedSpecies object

        explainedSpecies_without <- explainedSpecies

        # we retrive the hit groups that the focal hit group explains
        explainedByFocal <- explainedHitgroups[[focalHitgroup]]

        # we remove the hit group from the explanatory ones
        explainedSpecies_without[explainedByFocal] <- removeHitGroupFrom(focalHitgroup, explainedSpecies[explainedByFocal])

        # but we may not remove the hit group if doing so makes is so that previous
        # focal hit groups that are considered explained can no longer be explained

        # the previous hit groups to investigate are those that were explained
        # by the focal hit group
        toInvestigate <- intersect(
            x = which(explained | required),
            y = explainedByFocal
        )

        if (any(toInvestigate)) {

            # we check if all these hit groups to can still be explained without the focal hit group
            stillExplained <- sapply(
                X = toInvestigate,
                FUN = isExplained,
                explainedSpecies = explainedSpecies_without
            )

            if (!all(stillExplained)) {

                # if not, we mark the hit group as required
                required[focalHitgroup] <- T
                cat(".") # progress indicator (though it is fast)

                # and we jump to the next hit group
                next
            }
        }

        # if the focal hit group is not required to explain any other, we validate its removal from
        # explanarory hit groups
        explainedSpecies <- explainedSpecies_without

        # and we consider the hit group as explained
        explained[focalHitgroup] <- T
    }
}

# WE ARE NOW FINISHED WITH THE COUNT OF HTT EVENTS -------------------------------------------------------------------



# we only retain the columns we need in the htt hit tablee
httHits <- httHits[, .(query, subject, q, s, sp1, sp2, f1, f2, superF, mrca, pID, length, qStart, qEnd, sStart, sEnd, ka, ks, length.aa, com, hitGroup)]

# we add two logical columns:
# - "keep" is TRUE for hit groups we retained (based on the contamination filter , sufficient divergent time and low Ks)
# - "independent" is TRUE is a hit group is not explained by other (i.e. a htt event that may be seen as "independent" or "direct")
httHits[, c("keep", "independent") := .(hitGroupStats$keep[hitGroup], !explained[hitGroup])]
writeT(httHits, "HTThitsAssessed.txt") # we save this whole table, including hitgroups that are not retained


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

# we remove hits and columns we don't keep
retainedHits <- httHits[keep == T, -c("keep", "mrca", "superF")]

# we replace column names with more user-friendly ones
setnames(
    x = retainedHits,
    old = c("query", "subject", "f1", "f2", "com", "superFName"),
    new = c("copy1", "copy2", "consensus1", "consensus2", "community", "superfamily")
)


writeT(retainedHits, "supplementary-data4-retained_hits.txt")
