

##%######################################################%##
#                                                          #
####           This scripts applies cliterion           ####
####              1 to hits from different              ####
####       communities, to potentially merge them       ####
####            into hit group (htt events).            ####
#                                                          #
##%######################################################%##

#This shares a lot of code with iterativeFirstClustering.R
#many of the commands replicate what the other script does, but there are subtleties.

source("HTvFunctions.R")
args = commandArgs(trailingOnly = TRUE)
nCPUs = as.integer(args[1])

hitList = readRDS("TEs/clustering/round2/hitsFor2ndClustering.RDS")

superF = splitToColumns(names(hitList), "_", 1)
nHits = sapply(hitList, nrow)
nHits = sort(tapply(nHits, superF, sum), T)		#to process the super families with more hits first

blastFiles = list.files("TEs/clustering/selectedHits", full.names = T)
names(blastFiles) =  gsub(".out", "", basename(blastFiles), fixed = T)		#names blast files with superfamily names (extracted from their names)

communityPairStats = function(superFam) {
  blastFile = blastFiles[superFam]
  blast = fread(blastFile,
                sep = "\t",
                header = F,
                drop = 4)
  blast[, V3 := as.integer(V3 * 1000L)]										#converts pIDs to integers as we did for the htt hits (to save memory)
  blast[, copyPair := 10 ^ 6 * V1 + V2]				#here we don't create a matrix of pIDs where subject and queries are rows and columns numbers, as it would be too big for this round. We create a unique subject-query identifier based on the fact that these are integers, and both lower than 10^6.
  groups = hitList[superF == superFam]			#the group of hits (one per mrca) for this superfamily
  
  processHitsOfGroup = function(group) {
    hits = copy(group)
    setorder(hits, com)													#sorting hits per community will ensures that the communities pairs are formed in the same "direction" (com1 vs com2 and never the opposite)
    uCopies = hits[, unique(c(query, subject))]							#the copies involved in the hits to cluster
    sub = blast[V1 %in% uCopies &
                  V2 %in% uCopies, .(copyPair, V3)]		#we selected the TE hits involving these copies
    sub = rbind(sub, data.table(copyPair = 0, V3 = 0L))					#we add a copyPair that does not actually exist (with 0 pID) at the end of this table. This trick will help us later
    nHits = nrow(hits)
    nBatches = ceiling(nHits ^ 2 / 2 ^ 28)			#we anticipate the number of hits we will have to compare and create batches that will perform 2^28 comparison at once (lower than the limit of list sizes in R)
    hitBatches = list(1:(nHits - 1L))
    if (nBatches > 1)
      hitBatches = splitEqual(hitBatches[[1]], n = nBatches)
    criterion_1 = function(hitBatch) {
      #this "connects" hits according to the criterion, similar to what we did in the first round
      pairs = data.table(V1 = rep(hitBatch, nHits - hitBatch),
                         V2 = unlist(lapply(hitBatch[1]:max(hitBatch) + 1L, function(hit)
                           hit:nHits)))					#makes all possible pairs of hits for this batch
      pairs[, c("q1", "s1", "com1", "inter1", "q2", "s2", "com2", "inter2") := data.table(hits[V1, .(query, subject, com, pID)], hits[V2, .(query, subject, com, pID)])]	#we retreive the ids of copies involved in the 2 hits (2 per clade)
      pairs[q1 > q2, c("q1", "q2") := .(q2, q1)]			#to create query-subject pair idenfiers, we first make sure that the left query number is always lower than the rigth one (as it is the case for the blast files)
      pairs[s1 > s2, c("s1", "s2") := .(s2, s1)]			#same for subjects
      pairs[, c("qPair", "sPair") := .(10 ^ 6 * q1 + q2, 10 ^ 6 * s1 + s2)]		#creates the same pair identifiers we created for the blast results
      pairs[, c("intra1", "intra2") := .(sub[match(qPair, copyPair, nomatch = .N), V3], sub[match(sPair, copyPair, nomatch =
                                                                                                    .N), V3])]			#so we can retreive their identities. The nomatch argument makes it so that the pID of the last row in sub (wich has pID 0) is return. This seeds things up a little
      pairs[q1 == q2, intra1 := 100000L]
      pairs[s1 == s2, intra2 := 100000L]
      #if copies are the same, they are not in the blast hits so we set their identity to 100%
      pairs[, maxIntra := pmax(intra1, intra2)]
      pairs[, connected := inter1 < maxIntra |
              inter2 < maxIntra]			#the 2 hits will be connected in the best intra-clade identity is higher than one inter-clade identity (that of hits)
      cat("*")
      pairs[, .(
        allDiff = sum(!connected),
        diff = sum(maxIntra > 0L &
                     !connected),
        valid = sum(maxIntra > 0L),
        tot = .N
      ), by = .(com1, com2)]		#returns statistics for each community pair (note that a pair may also involve the same community twice, which helps to assess how much hits within a community are connected)
      #allDiff is the number of pairs of hits not passing the criterion, diff is
      #the same, but for hits whose copies have some degree of intra-clade
      #indentity, valid is the number of pairs of hits whose copies have non-zero
      #intra-clade indentity, tot is the number of pairs of hits evaluated
    }
    statsPerCom = rbindlist(mclapply(
      hitBatches,
      criterion_1,
      mc.cores = min(nCPUs, length(hitBatches)),
      mc.preschedule = F
    ))
    statsPerCom[, .(
      allDiff = sum(allDiff),
      diff = sum(diff),
      valid = sum(valid),
      tot = sum(tot)
    ), by = .(com1, com2)]
  }
  stats = rbindlist(lapply(groups, processHitsOfGroup))
  writeT(stats,
         stri_c("TEs/clustering/round2/", superFam, ".all.txt"))			#writes to disk for safety, this is also returned by he function
  print(paste("done", superFam))
  stats
}

res = lapply(names(nHits), communityPairStats)
writeT(rbindlist(res), "TEs/clustering/round2/all.txt")
