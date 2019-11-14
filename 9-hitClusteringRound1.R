
#Cluster TE-TE hits iteratively to delinate HTT events. This is teh first step where clusters hits from the same pair of "young" caldes: clades of less than 40 My, in which we did not even look for HTT and did not measure any BUSCO Ks
#then we will agregate the clusters that may reflect the same HTT. 
#clustering and agregation will be based on criterion 1: pID of TEs between clade is lower than pID within clades

source("HTvFunctions.R")
httHits = fread("occ200Ks05.txt")			#hits to cluster, filtered in step 7-TEKsAndHTTfilter.R

#we will only compare hits that involve the same TE super family and species that coalesce to the same MRCA, or else they cannot result from the same HTT, and there is not even the notion of a shared "clade" between hits to use criterion 1. 
#to compare TE copies wihtin clade to avaluate criterion 1, if is easier to swap sp1 and sp2 so that sp1 always corresponds to the same clade (and sp2 to the other) among the two sister caldes that diverged from an MRCA. So we will compare "queries" (copies from the left-hand clade) with each others and "subjects" (copies from the right-hand clade) with each others (not queries to subjects, that is already done with the pID in the HTT hits)
# so for each hit between species coalescing to an MRCA, we find its clade (the one among the two that diverged from the mrca)

tree = read.tree("timetree.nwk")
mrcas = httHits[, unique(mrca)]			#the different MRCAs of the two species involved in each HTT hits
edges = data.table(tree$edge)			#this 2-column table gives the children (subclades, in column V2) of each node (MRCA, column V1)
subClades = edges[V1 %in% mrcas, V2]	#these nodes are the subclades of each mrca in the httHits
subClades = tipsForNodes(tree, subClades, names = T)	#retreives the species of these subclades
subClades[,MRCA := edges[match(node, V2), V1]]		#and the parent node (mrca) of each subclade
subClades[,comb := stri_c(tip, MRCA, sep = " ")]	#we combine tip (species) name and MRCA number, to match these combinations in the httHits

httHits[,c("sc1","sc2") := .(subClades[chmatch(stri_c(sp1, mrca, sep = " "), comb), node], subClades[chmatch(stri_c(sp2, mrca, sep = " "), comb), node])]		#assigns subclade numbers (node numbers of the tree) to species in each hit, and below we ensure that the subclade with lower number is always on the left.
#This requires swapping many columns of the table:
httHits[sc1 > sc2, c("query","subject","sp1","sp2","f1","f2","qStart","qEnd","sStart","sEnd","sc1","sc2") := .(subject, query, sp2, sp1, f2, f1, sStart, sEnd, qStart, qEnd, sc2, sc1)]

writeT(httHits, "occ200_subClades.txt")

#to get within-clades pID of TE hits, we use the self-blastn outputs, which are huge (hundreds of GBs total). So we filter them to select only the hits we need: those between copies of the same (sub)clade, to evaluate criterion 1 (within-clade pID lower than between-clade pID). We already have between-clade pIDs, since these are in the httHits table
#for this, we need copy integer identifiers as these were used in the self blast to save space
copyIDs = fread("TEs/clustering/selectedCopiesKs05occ200.IDs.txt")				#we retreive copy integer IDs
copyIDs[,name := copyName(copy)]
httHits[,c("q","s") := .(copyIDs[chmatch(query, name), id], copyIDs[chmatch(subject, name), id])]

#since we will compare only copies within superFamilies and MRCA (the combination of which we call a "group"), we extract them
queries = httHits[, unique(q), by = .(superF, mrca)]
queries = split(queries$V1, queries[,paste(superF, mrca)])  #the "query" copies for each superfamily and mrca, that we we compare together as these are from the same (left-hand) subclade
ql = sapply(queries, length)								#we will sort groups by decreasing number of copies in queries
queries = queries[order(ql, decreasing = T)]
subjects = httHits[,unique(s), by = .(superF, mrca)]		#extracting subjects TE copies
subjects = split(subjects$V1, subjects[,paste(superF, mrca)])
subjects = subjects[order(ql, decreasing = T)]

#we extract self-blast hits per super family, so we sub-split the TE copies (subjects and queries) lists by super Families
superF = gsub("/",".", splitToColumns(names(queries), " ",1), fixed = T)		#extract super family names and replaces slashes with periods as they will be used to generate file names
queries = split(queries, superF) ; subjects = split(subjects, superF)

nq = sapply(queries, function(x) length(unique(unlist(x))))			#and we sort these list by decrasing numbers of copies (not required, but I prefer starting with the large ones)
queries = queries[order(nq, decreasing = T)] ; subjects = subjects[order(nq, decreasing = T)] 

blastFiles = list.files("TEs/clustering/blastn/done", pattern = "", full.names = T)		#path of self blast result files generated in set 8-prepareClustering.R
blastFiles = blastFiles[order(file.size(blastFiles), decreasing = T)]					#again, ordered by decreasing size
blastFiles = split(blastFiles, splitToColumns(basename(blastFiles), "_",1))				#splits blast files per super families (large superfamilies were processed in parallele self-blast searches)

filteredHitFolder = "TEs/clustering/selectedHits/"	#were the selected hits will go
dir.create(filteredHitFolder)		

getHits = function(superF) {			#extracts the relevant TE-TE hits for a given super family
	files = blastFiles[[superF]]
	import = function(file) {
		hits = fread(file, sep ="\t", header = F, nThread = 1, select = c(1:3, 9))
		hits = Map(function(q, s)  hits[(V1 %in% q & V2 %in% q) | (V1 %in% s & V2 %in% s),], queries[[superF]], subjects[[superF]])		#select hits between copies of the same subclade
		hits = rbindlist(hits)
		hits = hits[!duplicated(data.table(V1, V2)),]		#some rare hits may have been selected severale times, we remove duplicates
		cat(".")											#some progress indicator as files can be quite big
		hits
	}
	hits = mclapply(files, import, mc.cores = min(20, length(files)), mc.preschedule = F)		#applies the function for each blast file of the super family, in parallel
	writeT(rbindlist(hits), stri_c(filteredHitFolder, superF,".out"), col.names = F)			#writes a single file of selected hits
	print(paste("done", superF))
	NULL	
}
	
m = lapply(names(queries), getHits)


#we now cluster hits between species of the same pair of "young" clades that diverged within the last 40 My. Wihtin these clades, we did not even search for HTT (no blast) and we didn't measure Ks at buscos.
# These clusters of hits will be the initial ones (the "communities") that we will agregate in a second round
#so we assign species involved hits to the young clades
clades = cladesOfAge(tree, 40, withTips = T)				
httHits[,c("C1","C2") := .(clades[chmatch(sp1, tip), node], clades[chmatch(sp2, tip), node])]		#assigning species to clades

# to cluster hits, we prepare a table with only the columns we need.
httHits[,c("hit", "group") := .(1:.N, stri_c(stri_replace(superF, ".", fixed = "/"), C1, C2, sep = "_"))]		#clustering of hits will be done within each "group" (hits involving a given pair of young clades in a super family). Hits are indentified by their row indices. We again replace slashes in super family names with periods
conv= httHits[,data.table(hit, pID = as.integer(pID*1000L), group, query = q, subject = s)]		# We get the column we need to the clusering convert blast pIDs to integers for efficiency in the clustering 
hitList = split(conv, conv$group)				#split hits in groups 
hitList = hitList[order(-sapply(hitList, nrow))]
hitList = hitList[sapply(hitList, nrow) > 1L]			#if there is just one hit, no clustering needs to be done for this group

dir.create("TEs/clustering/round1/")
saveRDS(hitList, file ="TEs/clustering/round1/hitsForFirstClusteringKs05occ200.RDS")			#this list will be used in the clustering jobs. Below we plan the jobs:

groupList = data.table(group = names(hitList), nHits = sapply(hitList,nrow), superF = splitToColumns(names(hitList), "_", 1))		#number of hits to compare 2by2 per "group"
perSuperF = groupList[,.(nHits = sum(nHits), .N), by = superF]																#number of hits and groups per super family. We'll do the clustering by super families, so that a blast output (which concerns a super family) need only be imported once for all the groups within it
perSuperF[,batch := 0L]; perSuperF[nHits > 30000, batch := 1:.N]						#attributes batches (job numbers) to super families. Note that several super families with low number of hits may be processed in the same job
writeT(perSuperF, "TEs/clustering/round1/batchesKs05occ200.txt")							#this will be used to manage jobs in iterativeFirstClustering.R

# example of job number 6 using 5 cpus:
system('sbatch --mail-type=BEGIN,END,FAIL --cpus-per-task=5 --mem=100G --wrap="Rscript iterativeFirstClustering.R 6 5"')

###### collecting results
files = list.files("TEs/clustering/round1", pattern = "allGroups.txt", full.names = T)
groups = rbindlist(lapply(files, fread, header = T, sep = "\t"))
groups[, ucomm :=  toInteger(stri_c(comm, group, sep = " "))]			#generates unique integer community ids (accross all groups)

httHits[groups$hit, com := groups$ucomm]					#adds these identifiers to the HTT hit table 
httHits[is.na(com), com := 1:.N + max(groups$ucomm)]		#com is NA for hits were not even clustered because there are just one per super family and clade pair, so we give them unique community numbers.

writeT(httHits, "occ200_comm.txt")
