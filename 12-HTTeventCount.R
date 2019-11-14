#counting HTT events required to explain the hit groups (see methods for details, as this would be quite long to summarize in a comment)
source("HTvFunctions.R")
httHits = fread("oc200HitGroup.txt")				


##### we need to find copies that are similar between hit groups to assess whether they can be "exlpained" by others. We do it per superfamily 
hitList = split(httHits[,.(q, s, pID = as.integer(pID * 1000), hitGroup)], gsub("/",".", httHits$superF, fixed = T))	#so we slit the hit table per super family, retaining only the columns we need (copy ids) and converting to integer as appropriate
saveRDS(hitList, file = "TEs/clustering/forHitGroupPairs.RDS")

system("Rscript hitGroupPairHomologousCopies.R 5")				#this script finds homologous TE copies between any 2 hit groups
# we can explain an apparent htt (hit group) by others that may have brought TEs in both clades involved in hit group.
hitGroupPairs=fread("TEs/clustering/networks/all.hitGroupPairSharedCopies.txt")			#retreiving results from hitGroupPairHomologousCopies.R. "hitGroup1" is a hit group with a TE "copy" having some degree of similarity ("pID") with a (unnammed) copy of "hitGroup2". hitGroup2 may therefore "explain" hitGroup1. But for this, hitGroup2 has to fulfil some conditions

#First, we assess whether hitGroup2 involves species related to those of hitGroup1. We assess if a clade of hitGroup1 is the same as, nested in, or encompassing, a clade of hitGroup2. Note that this doesn't require that species are shared. "clade 1" is always the clade on the left of table, (hence containing the "query" copies) in the hit group hits
tree = read.tree("timetree.nwk")
clades = httHits[, .(node1 = MRCA(tree, unique(sp1)), node2 = MRCA(tree, unique(sp2))), by = hitGroup]		#the pairs of clades involved each hit group (transfer), encoded as node numbers in the species tree
mrcas = httHits[match(1:max(hitGroup), hitGroup), mrca] #this vector will make the correspondance between a hit group (its index) and an mrca (its value) 
node1 = clades[match(1:max(hitGroup), hitGroup), node1] ; node2 = clades[match(1:max(hitGroup), hitGroup), node2]		#same principle as above for the 2 clades composing the hit groups. Doing this speeds up the command below as this avoid performing a match on the huge hitGroupPairs table

hitGroupPairs[,c("n1_1","n2_1", "n1_2","n2_2", "mrca1") := .(node1[hitGroup1], node2[hitGroup1] ,node1[hitGroup2], node2[hitGroup2], mrcas[hitGroup1])]		#we retreive these elements for each hitGroupe pair
mat = nestedClades(tree)			#matrix that is TRUE if a clade (column or row index) is nested in, or includes, another (row or column index). The diagonal is also TRUE
mat2 = nestedClades(tree, F)
hitGroupPairs[,c("nested1","nested2","nestedMRCA") := .(mat[cbind(n1_1, n1_2)] | mat[cbind(n1_1, n2_2)], mat[cbind(n2_1, n1_2)] | mat[cbind(n2_1, n2_2)], mat2[cbind(n1_2, mrca1)] | mat2[cbind(n2_2, mrca1)])] #these new colums will be TRUE if the clade1 of hitGroup1 is nested, or indludes, a clade of hitGroup2. nested2 does the same for clade2 of hitGroup1, and if the mrca of hitGroup1 is nested within a clade of hitGroup2 (in which case we will consider that hitGroup2 cannot exlpain hitGroup1)
hitGroupPairs[,c("n1_1","n1_2", "n2_1","n2_2","mrca1") := NULL]


#We also impose that the best identity that a copy of hitGroup1 has with any copy of hitGroup2 is not lower that the lowest identity the copy has within hitGroup1
minCopyIDs = rbind(httHits[,min(pID), by = .(hitGroup, copy = q)], httHits[,min(pID), by = .(hitGroup, copy = s)])		#min pID each copy has within each hit group
minCopyIDs[,copyHitGroup := copy*10000 + hitGroup]			#integer to identify the copy-hit group pair
hitGroupPairs[,copyHitGroup1 := copy*10000 + hitGroup1]
hitGroupPairs[, minID := minCopyIDs[match(copyHitGroup1, copyHitGroup), as.integer(V1*1000)]]		# this allows placing this pID in the hitGroupPair table

temp = httHits[,.(cop = c(q, s), sp = chmatch(c(sp1, sp2), tree$tip.label))]		#we also retreive the species of each copy of hit group 1, coded as integer
hitGroupPairs[,sp := temp[match(copy, cop), sp]] ; rm(temp)

hitGroupPairs[,c("inClade1", "inClade2") := .(copyHitGroup1 %in% httHits[,q*10000 + hitGroup], copyHitGroup1 %in% httHits[,s*10000 + hitGroup])]		#we also determine if a copy belongs to a species from the left or right clade for each hit group

# we also determine if hitGroup1 and hitGroup2 were not considered independant based on the criterion 1 of hit clustering that separated hit groups (not counting the absence of homology at protein regions). If that is the case, we will decide that hitGroup2 cannot explain hitGroup1. 
mStats = fread("TEs/clustering/comPairStats.Ks05.OCC200.txt")		#file generated in step 10-hitClusteringRound2.R
crit1 = mStats[,sum(links)/sum(tot), by = .(hitGroup1, hitGroup2)]
crit1[, pair := ifelse(hitGroup1 < hitGroup2, hitGroup1*10000 + hitGroup2, hitGroup2*10000 + hitGroup1)]
hitGroupPairs[, hitGroupPair := ifelse(hitGroup1 < hitGroup2, hitGroup1*10000 + hitGroup2, hitGroup2*10000 + hitGroup1)]
hitGroupPairs[, indep := hitGroupPair %in% crit1[V1 < 0.05, pair]] ; rm(mStats, crit1)

#to determine if a hit group can be explained by others, we will inspect the less "reliable" hit groups first
hitGroupStats = fread("TEs/clustering/hitGroupStats.txt")		#imports file generated in the previous step
setorder(httHits, hitGroup, -pID)   							#puts best hits on top, for each hit group 
pIDQ = httHits[!duplicated(data.table(hitGroup, q)), sum(pID), by = hitGroup]		#the reliability score of a hit group is based on the pID of best hits of copies involved (see Methods)
pIDS = httHits[!duplicated(data.table(hitGroup, s)), sum(pID), by = hitGroup]
hitGroupStats[,score := pmin(pIDQ$V1, pIDS$V1)]
orderedHitGroups = hitGroupStats[keep == T, hitGroup[order(score)]]		#the retained hitGroups, ordered by reliability score

#we now select rows passing filters. A copy of a hit group can be brougth by another HTT (hitGroup2) only if the hitGroup2 copy belongs to a species that is the same (or related) to the one carrying the copy. We also consider that hitGroup2 cannot explain hitGroup1 if this pair of hit groups did not pass criterion 1 or if the mrca of hitGroup1 is nested within a clade of hitGroup2
#We also discard hit groups that were considered unreliable in the previous step. This could have been done earlier, but we prefer applying at the last moment.
sel = hitGroupPairs[pID >= minID & indep == F & ((inClade1==T & nested1==T) | (inClade2==T & nested2==T)) & nestedMRCA == F & hitGroup2 %in% orderedHitGroups & hitGroup1 %in% orderedHitGroups, .(hitGroup1, hitGroup2, sp, c1 = nested1 & inClade1, c2 = nested2 & inClade2)]	

#as we will evaluate whether hitGroup2 could have brought TE copies in species composing a clade of hitGroup1, we do the following:
#for each retained hit group, we list the species (for each clade) that have copies similar to those of other hit groups harboring related species. This is done in the function below
sharedSpeciesOfClade = function(clade) {		#clade is either 1 (left clade) or 2
	if(clade == 1) {
		temp = unique(sel[c1 == T, .(hitGroup1, hitGroup2, sp)]) #a table of species of left clade that have copies similar to other hit groups (having species related to clade 1)
	} else {
		temp = unique(sel[c2 == T, .(hitGroup1, hitGroup2, sp)])			#same as above for right clade
	}
	spList = split(temp$sp, temp[, hitGroup1*10000 + hitGroup2])			#splits species by pair of hit groups (pair encoded as an integer for speed)
	pair = as.integer(names(spList))			#we will re-split this list for each hitGroup1, so we retreive hit group 1 from the integer encoding the hit group pair
	hitGroup1 = as.integer(pair/10000)
	names(spList) = pair - hitGroup1*10000			#same for hit group 2, which we use to name elements of the list
	spList = reList(split(spList, hitGroup1), max(orderedHitGroups))		#so spList[[n]] gives the "explanatory" hit groups (hitGroup2) that have copies similar to those of left-clade species of hit group #n, and within each of these hit groups (2nd level of the list), the species of hit group #n (left clade) to which the explanatory hit group have similar copies
	temp = unique(temp[,.(hitGroup1, hitGroup2)]) 			#We do more or less the reciprocal of the above list, expect it does not list species
	contributedClade = reList(split(temp$hitGroup1, temp$hitGroup2), max(orderedHitGroups))	#so contributedClade[[n]] list the hit groups whose left-clade copies are "explained" or "contributed" by hit group #n (i.e., are similar to those of hit group #n). 
	list(spList, contributedClade)		#we return both lists
}

res = sharedSpeciesOfClade(1) ; clade1Sp = res[[1]] ; contributedClade1 = res[[2]]
res = sharedSpeciesOfClade(2) ; clade2Sp = res[[1]] ; contributedClade2 = res[[2]]

#to evaluate whether a hit group can be adequately explained by others, we also need to list all the species pairs (constituting the hits) it involves. Species of the left clade are on the left in pairs. Pairs are encoded as integers, to speed up subsequent commands
temp = httHits[,unique(chmatch(sp1, tree$tip.label)*1000L + chmatch(sp2, tree$tip.label)), by = hitGroup]		
spPairsForHitGroup = reList(split(temp$V1, temp$hitGroup)) ; rm(temp)		#for each hit group (index of this list), gives the species pair
NspPairsForHitGroup = sapply(spPairsForHitGroup, length)					#we will also need to know the number of species pair per hit group

#two functions we use to evaluate hit groups : 
isExplained = function(hitGroup) {									#tells if a focal hit group (integer number) can be explained, returns T or F
	exp1 = clade1Sp[[hitGroup]]										#the species of its left clade whose copies are similar to other hit groups
	exp2 = clade2Sp[[hitGroup]]										#same for its right clade
	sp1 = unique(exp1, use.names = F))								#we unlist these species (there were splitted by hit groups)		
	sp2 = unique(unlist(exp2, use.names = F))						
	if(!any(sp1) | !any(sp2) | length(union(names(exp1), names(exp2))) < 2L) return(F)		#we can already return F if at last one clade of the focal hit group has no such species, or if these species are not present in at least 2 other hit groups (a hit group must be explained by at least two others)
	spPairs = rep(sp1, each = length(sp2))*1000L + rep(sp2, length(sp1))		#based on these species, creates all possible pairs of species (left-right)
	found = sum(spPairs %in% spPairsForHitGroup[[hitGroup]])					#the number of these pairs that are actually represented in the hits of the hit group
	NspPairsForHitGroup[hitGroup] <= found | found == NspPairsForHitGroup[hitGroup]			#applies our criteria for "explanation". Basically, enough pair of species must be explained by other hit groups
}

removeHitGroupFrom = function(hitGroup, ls) {						#removes a hit group (integer) from a list (argument ls) of potential explanatory groups
	charHitGroup = as.character(hitGroup)
	lapply(ls, function(x) x[!names(x) %chin% charHitGroup])
}

#we are now ready to go
explained = required = logical(nrow(hitGroupStats)) 			#logical vectors that will be TRUE for hit groups that are "explained" hence removed (as indices of these vectors) or "required" to explain other hit groups (see below)
for (i in 1:length(orderedHitGroups)) {			#takes hit groups by order or reliability
	hitGroup = orderedHitGroups[i]
	if(isExplained(hitGroup)) {				#if a hit group can be explained by others
		expl1 = intersect(which(explained | required), contributedClade1[[hitGroup]])		#retreives the already explained or required hit groups that are partly explained by this hit group, for the left clade
		expl2 = intersect(which(explained | required), contributedClade2[[hitGroup]])		#same for the right clade
		expl = union(expl1, expl2) 
		if(any(expl)) {											#if there are any
			clade1Spb = clade1Sp ; clade2Spb = clade2Sp			#makes backups of the lists we use (we could do differently, but this allows reducing the number of arguments passed to isExplaiend())
			if(any(expl1)) clade1Sp[expl1] = removeHitGroupFrom(hitGroup,  clade1Sp[expl1])		#and removes the hit group from the explanatory ones (left clade)
			if(any(expl2)) clade2Sp[expl2] = removeHitGroupFrom(hitGroup,  clade2Sp[expl2])		#same for the right clade
			if(!all(sapply(expl, isExplained))) {		#and check if all the hit groups it contributes to can still be explained. If not, we will not remove the hit group, since it is required to explain others
				clade1Sp = clade1Spb ; clade2Sp = clade2Spb		#so we restore the lists
				required[hitGroup] = T 		#mark the hit group as required
				cat(".") ; next			#and jump to the next hit group
			}
		}		#else, we remove the hit group from explanatory ones
		toChange = intersect(orderedHitGroups[-(1:i)], contributedClade1[[hitGroup]])	#the other hit groups that the hit group explains (left clade). We only list those that are not inspected yet (as the others have alreayd been dealt with above)
		if(any(toChange)) clade1Sp[toChange] = removeHitGroupFrom(hitGroup, clade1Sp[toChange])
		toChange = intersect(orderedHitGroups[-(1:i)], contributedClade2[[hitGroup]])
		if(any(toChange)) clade2Sp[toChange] = removeHitGroupFrom(hitGroup, clade2Sp[toChange])
		explained[hitGroup] = T
	}
}

##### now cleaning the httHits table

httHits = httHits[,.(query, subject, q, s, sp1, sp2, f1, f2, superF, mrca, pID, length, qStart, qEnd, sStart, sEnd, ka, ks, length.aa, com, hitGroup)]		#removes columns we no longer need
httHits[,c("keep","independent") := .(hitGroupStats$keep[hitGroup], !explained[hitGroup])]		#keep is TRUE for hit groups we retained (based on the contamination filter , sufficient divergent time and low Ks) and independnat is TRUE is a hit group is not explained by other
writeT(httHits, "HTThitsAssessed.txt")		#we save this whole table, including hitgroups that are not retained

#prepares the supplementary dataset of retained hits

corres = fread("superF.txt", header = F, col.names = c("superFam","subClass","newName"))		#correspondance between  repeat modeler super families and more common super family names, which we will use from now on. This file is provided with the scripts
httHits[,superFName := corres[chmatch(superF, superFam), newName]]			#adds common super family names to the hits
retainedHits = httHits[keep == T, -c("keep","mrca","superF")]				#removes hits and columns we don't keep
setnames(retainedHits,c("query","subject","f1","f2","com","superFName"), c("copy1","copy2","consensus1","consensus2","community","superfamily"))		#replace column names with more user-friendly ones
writeT(retainedHits, "data4-retained_hits.txt")								#the list of HTT hits associated with the paper (supplementary dataset). This one only includes retained hit groups

