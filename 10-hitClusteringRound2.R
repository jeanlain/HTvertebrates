
#### for the 2nd round of clustering, we apply criterion 1 to hits from different communities within the same super family and MRCA (a "group"). The procedure is similar to what we did in the previous step
source("HTvFunctions.R")
httHits = fread("occ200_comm.txt")			#HTT hits with community IDs from the previous step

httHits[,group := stri_c(gsub("/",".", superF, fixed = T), "_", mrca)]  #we will confront any two hits from the same superfamily and mrca. Here, the two hits need bot be between the same pair of young clades.
hitList = split(httHits[,.(hit = 1:.N, query = q, subject = s, pID = as.integer(pID*1000L), group, com)], f = httHits$group)
nCom = sapply(hitList, function(x) length(unique(x$com)))
hitList = hitList[nCom > 1L]									#there is nothing to cluster if there is just one community of hits per group
hitList = hitList[order(sapply(hitList, nrow), decreasing = T)]

dir.create("TEs/clustering/round2")
saveRDS(hitList, file ="TEs/clustering/round2/hitsFor2ndClustering.RDS")

#we now assess criterion 1 for hits belonging to different communities :
system("Rscript iterativeSecondClustering.R 20")		#this job is shorter (no clustering algorithm involved), so launched on the Poitiers server

###### collecting results. We now evaluate the connexions between hit communities that may represent the same HTTs
stats = fread("TEs/clustering/round2/all.txt")	#imports statistics for every relevant pair of hit commmunities (number of hit pairs compared and so on. See iterativeSecondClustering.R)
stats[,links := tot - allDiff]			#the number of hit pairs that pass criterion 1 between two communities

## we evaluate whether 2 communities can represent the same HTT, based on "citerion 2", which relies on the inferred age of the transfer. To reflect the same HTT the transfers (represented by the 2 clusters of hits) should not be more recent that both clades and 
#we infer the maximal times of transfers from the mean Ks between copies at hits of the same community
#we compare it do the Ks of buscos of corresponding clades
#so we determine the clades involved. For 2 clusters of hits between young clades A-B (community 1) and C-D (community 2) both coalescing to the same MRCA (mrca of A and B is the same that of C and D), 2 clades are involved. One is composed of A-C and the order of species B-D, if A and C belong to one of the 2 subclade diverging from the MRCA. Remember that we took care to place species of one subclade in column sp1 and species of the other subclade in column sp2, for each mrca. 

#for each community, we compute the mean Ks of hits, and other some associated metrics (that we won't actually use afterwards)
httHits[ks < 0, ks := 0]		#for some reasons rare ks values can be very sightly < 0, we fix that
ksStats = httHits[, .(nq = length(unique(query)), ns = length(unique(subject)), maxKs = max(ks), meanKs = mean(ks)), by = .(com, mrca)]

#we get the core gene Ks threshold we used for to filter hits between species of a given MRCA
Ks = fread("gunzip -c Ks200AAnoRedundancy.txt.gz")					#file of filtered BUSCO Ks (200 AA alignments and one score per BUSCO gene per pair of clades) generated in step 5-coreGeneKs.R
quant = Ks[,.(min = min(Ks), q05 = quantile(Ks, 0.005), m = mean(Ks), divTime = divTime[1]), by = clade]
ksStats[, c("threshold", "minKs") := quant[match(mrca, clade), .(q05, min)]]

#we retrieve the 2 clades (A-C and B-D in the explanation above) that are involved in every pair of communities
tree = read.tree("timetree.nwk")
mrca = mrca(tree)
species = cbind(httHits[match(stats$com1, com), cbind(sp1, sp2)], httHits[match(stats$com2, com), cbind(sp1, sp2)])		#for each community pair, we retreive the species involved (2 per community). These 2 species are just the ones that appear in the first hit in the httHits table for a given community, but it's ok since sp1 species from a hit community are all from the same young clade (same for sp2 species)  
stats[,c("c1","c2") := .(mrca[species[,c(1,3)]], mrca[species[,c(2,4)]])]
quant = rbind(quant, data.table(clade = 0L, min = 0, q05 = 0, m = 0, divTime = 0))			#we add this last row to ensure that the clades not in the quant table (those that were too young) will have Ks stats of zero. (see how the "nomatch" argument of match() below is set)
stats[,c("Ksc1", "Ksc2") := .(quant[match(c1, clade, .N), q05], quant[match(c2, clade, .N), q05])]				#adds the 0.5% quantile  BUSCO Ks for the clades in the table of community pairs
stats[,c("Ks1",  "Ks2") := data.table(ksStats[match(com1, com), meanKs], ksStats[match(com2, com), meanKs])]	#and we do the same for the mean Ks of hit communities involved

stats[,crit2 := (pmin(Ks1, Ks2 ) >= pmin(Ksc1, Ksc2) & pmax(Ks1, Ks2) >= pmax(Ksc1, Ksc2))]    #applying criterion 2
stats[,crit1 := links/tot > 0.05]															#and criterion 1 needs to be passed by more than 5% of hit pairs to consider the communities as representing the same HTT. But we are not done yet (see next)

########### to check whether different clusters of hits could represent the retention of different parts of a TE transferred once, we compare their protein regions. For this we reuse the blastx results of copies against repeat proteins (done in step 6-filterTEhits.R)

blastx = fread("TEs/blastx/all.copies.successiveBlastx.out", header = T, sep ="\t")				
copyIDs = fread("TEs/clustering/selectedCopiesKs05occ200.IDs.txt")				#we retreive copy integer IDs, to speed up some functions 
blastx[,query := copyIDs[chmatch(query, copy), id]]
blastx = blastx[!is.na(query)]					#we only retain relevant hits (not those involving copies that we did not retain in the pipeline)

#we blast proteins similar to copies against themselves. This will be used to determine whether two copies have some homologies, that is, is they have homology to protein regions that are themselves homologous. This should be quite sentitive.
prot = blastx[!is.na(query),unique(subject)]					#we ignore copies that were not from ourr selection
protSeqs = readAAStringSet("TEs/blastx/RepeatPeps.lib")			#the protein database used by repeat modeler (retreive from repeatmasker)
names(protSeqs) = splitToColumns(names(protSeqs), " ",1)		#discards sequence descriptions
dir.create("TEs/blastp")
writeXStringSet(protSeqs[prot], "TEs/blastp/involvedProtOCC200Ks05.fas")	#we create a fasta off all proteins that are similar to TEs 
system("makeblastdb -in TEs/blastp/involvedProtOCC200Ks05.fas -dbtype prot")
system("blastp -query TEs/blastp/involvedProtOCC200Ks05.fas -db TEs/blastp/involvedProtOCC200Ks05.fas -max_target_seqs 1500 -outfmt 6 -num_threads 5 -evalue 1e-4 -out TEs/blastp/involvedProtOCC200Ks05.self.out")		#we set a larger max_target_seqs than the default to ensure that all possible hits are reported

dt = httHits[, .(copy = unique(c(q, s))), by = .(com, mrca)]							#all copies involved in each community. We also need the mrca column 
protHits = dt[, blastx[query %in% copy], by = .(com, mrca)]								#retreiving blastx hits for copies in each community. This takes a bit of time
comPerMRCA = protHits[,length(unique(com)), by = mrca]									#counts the number of communities that have copies matching proteins for each MRCA (within which hits were clustered). 
protHits = protHits[mrca %in% comPerMRCA[V1 > 1L, mrca]]								#if there is just one such community, there is nothing to do so we can ignore these

#imports results of the self blastp launched above
#since we did a self blastp, we may remove reciprocal hits. We retain the ones with best score (reciprocal hits should have the same score, but we never know)
blastp = fread("TEs/blastp/involvedProtOCC200Ks05.self.out", header = F, col.names = c("query","subject","pID", "length","mismatches","gapOpen","qStart","qEnd","sStart","sEnd","evalue","score"))
#we make every hit reciprocal (forward and reverse, this will be useful later). In theory, they should be all be there already since we did a self-blast. But some reciprocal hits are not found for some reasons, and some may have different score than "forward" hits
blastp = removeReciprocal(blastp, removeSelf = F)				#this selects the best HSP in the same subject-query pair (including reciprocal hits)
blastpR = copy(blastp[query != subject])						#we now make a blast result data table with query and subject reversed
blastpR[, c("query","subject","qStart","qEnd","sStart","sEnd"):=.(subject, query, sStart, sEnd, qStart, qEnd)]
blastp = rbind(blastp, blastpR); rm(blastpR)
blastp[, pair := paste(query, subject)]

#for a community, copies may have several protein-coding regions that are overlapping or adjacent on a given protein. We can combine these regions, so we get larger, non redundant, regions that should more likely represent longer ancestral elements
protRegions = combineRegions(protHits[,.(stri_c(com, subject, sep=" "), sStart, sEnd)], distance = 10L)		#combining protein regions that are distant by 10 aas or less on a protein in a given community
protRegions = combineHomologous(protRegions, blastp, protSeqs)										#see function for details
writeT(protRegions, "TEs/clustering/protRegionsInCommunities.txt")

### now comparing communities at protein regions 
setorder(protRegions, commID)
nr =nrow(protRegions) ; rows = 1:(nr-1L)				#we will generate every possible pair of rows from this table in order to compare protein regions with each others and assess how they overlap
pairs = data.table(V1 = rep(rows, nr-rows), V2 = unlist(lapply(2:nr, function(x) x:nr)))		#this takes a lot of ram as there are lots of possible pairs
pairs[,c("com1","com2") := data.table(protRegions[V1, commID], protRegions[V2, commID])]
pairs[,pair := com1*10^5 + com2]			#creates unique identifier for community pairs
stats[,pair := com1*10^5 + com2]			#and simular ifentifiers for the other talbe containing statistics on community pairs. Note that com1 < com2 is both cases
pairs = pairs[pair %in% stats$pair] ; gc()

regionComp = pairs[, data.table(pair, protRegions[V1, .(prot, start, end)], protRegions[V2, .(prot2=prot, start2=start, end2=end)])]			#gets coordinate of protein regions we want to determine the homology
regionComp = proteinOverlap(regionComp, blastp, protSeqs)		#see function in HTvFunctions.R for details

compStats = regionComp[,.(max=max(interAll), sum=sum(interAll), maxL=max(length, na.rm=T), maxID=max(pID, na.rm=T)), by = pair]				#for a pair of communities, gets the longest inferred alignment and total length of all
mStats = merge(stats, compStats, by = "pair", all = T)															#merges these with stats
mStats[com1 != com2 & is.na(max), c("max", "sum", "maxL","maxID") := 0L]										#if a community was not even present in the protRegions, the values are NAs. We put zeros (we leave NAs for t1 == t2 as there was nothing to compare)

#obtain useful information per community: number of different proteins (some homologous may be not merged), length of longest protein region and length or longest protein
protRegions[,lp := nchar(protSeqs[prot])]																		#length of protein
comStats = protRegions[,.(nProt=length(unique(prot)), lr=max(end-start+1), lp=max(lp)), by = commID]				#number of different proteins, max protein length, etc.
mStats[,c("nProt1","lr1","lp1", "nProt2","lr2","lp2") := data.table(comStats[match(com1, commID), .(nProt, lr, lp)], comStats[match(com2, commID), .(nProt, lr, lp)])]
mStats[is.na(nProt1), c("nProt1","lr1","lp1") := 0L]
mStats[is.na(nProt2), c("nProt2","lr2","lp2") := 0L]

mStats[max < 100 & valid == 0, crit1 := T]			#two communities pass criterion 1 if there is no nucleotide identity of copies within clade (no self blastn hit) and the longest aligned protein part ("max") between copies of the 2 transfers is less than 100 aa. We leave the possiblity that copies that apparently have little to nothing in common could be different parts of the same TEs, which were retained by different clades after a single HTT. 


#we link communities that pass criteria 1 and 2 as being possibly the result of the same HTT. An HTT can only include communities that are all linked with each other (= complete-linkage clustering). Because linkage doesn't mean 2 communities MUST reflect the same HTT, it could just mean a lack of information to distinguish HTTs (while two communities not passing the criteria cannot possibly be included in the same HTT). So we have to take decisions as to which group to from. Within each clade pair and super family, we want to group first: communities that cover very similar (or the same) proteins (high maxID)
mStats = mStats[order(-maxID)]						#we place these pairs of transfers on top, since we will regroup from top to bottom	

hitGroups = mStats[,cLinkFromPairs(com1, com2, linked = crit2 & crit1)]						#does the clustering
mStats[,c("hitGroup1", "hitGroup2") := .(hitGroups[as.character(com1)], hitGroups[as.character(com2)])]	
mStats[hitGroup1 == hitGroup2, all(crit2 & crit1)]													#checks that all pairs of communities belonging to the same group passed both criteria

writeT(mStats, "TEs/clustering/comPairStats.Ks05.OCC200.txt")

httHits[,hitGroup := hitGroups[as.character(com)]]							#assigns hits to hit groups (new column added)
httHits[is.na(hitGroup), hitGroup := 1:.N + max(hitGroups)]					#any hit that is not in a hit group constitutes its own hit group (but these small hit groups will be removed by our filter in the next step)
writeT(httHits, "oc200HitGroup.txt")
