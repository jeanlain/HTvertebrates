##%######################################################%##
#                                                          #
####             This scripts performs the              ####
####           first clistering of hits into            ####
####        "communities" according to criteron         ####
####              1 of the study's method               ####
#                                                          #
##%######################################################%##


source("HTvFunctions.R")
library(igraph)

args = commandArgs(trailingOnly = TRUE)
job = as.integer(unlist(strsplit(args[1], ",")))
nCPUs = as.integer(args[2])

print(paste("job:", job))

hitList = readRDS("TEs/clustering/round1/hitsForFirstClusteringKs05occ200.RDS")				#the list of tables of hits (with conversion to integers)

batches = fread("TEs/clustering/round1/batchesKs05occ200.txt")

superF = batches[batch %in% job, superF]								#the superfamily(ies) to process in this batch
superFam = splitToColumns(names(hitList), "_", 1)
groups = names(hitList)[superFam %in% superF]							#selects the groups of hits to cluster for the corresponding super family(ies)

blastFiles = stri_c("TEs/clustering/selectedHits/", superF, ".out")			#imports corresponding blast file(s) of copies against themselves
blastFiles = blastFiles[file.exists(blastFiles)]
blastFiles = blastFiles[file.size(blastFiles) > 0]
blast = rbindlist(lapply(blastFiles, function(x)
  fread(
    x,
    header = F,
    drop = 4,
    sep = "\t"
  )))
blast[, V3 := as.integer(V3 * 1000L)]										#converts pID to integer as we did for the htt hits (to save memory)

hitCommunities = function(group) {
  #we do the clustering per 'group' (hits of pair of clades in a super family)
  hits = copy(hitList[[group]])									#retreives the hits to cluster. We copy the data table to avoid some side effect related to how data.table functions work
  ucopies = hits[, sort(unique(c(query, subject)))]				#all the different copies (ids) in these hits
  hits[, c("qid", "sid") := .(match(query, ucopies), match(subject, ucopies))]				#we convert these to smaller integers so that we can make matrices of pID between copies, where a copy id will be a row/column index in the matrix
  blast[, c("qid", "sid") := .(match(V1, ucopies), match(V2, ucopies))]				#we thus need to use the same ids for the copies in the blast output
  sub = blast[!is.na(qid) &
                !is.na(sid)]							#extract only the hits we need in these
  pIDmat = matrix(0L, length(ucopies), length(ucopies))			#makes a matrix of pID for each copy pair (with pID zero by default, when there is no hit)
  diag(pIDmat) = 100000L											#the diagonal is set to 100% pID (*1000), to represent perfect pID for a copy with itself
  pIDmat[cbind(sub$qid, sub$sid)] = sub$V3							#we can now fill the matrix, (both semi matrices)
  pIDmat[cbind(sub$sid, sub$qid)] = sub$V3
  
  nHits = nrow(hits)
  nBatches = ceiling(nHits ^ 2 / 2 ^ 28)		#we cannot always cluster all the hits at once due to the max size of vectors in R (and RAM required)
  hitBatches = list(1:(nHits - 1L))			#we may split the hits into several batches
  if (nBatches > 1)
    hitBatches = splitEqual(hitBatches[[1]], n = nBatches)
  
  criterion_1 = function(batch) {
    #this function "connects" hits according to criterion 1
    pairs = data.table(V1 = rep(batch, nHits - batch),
                       V2 = unlist(lapply(batch[1]:max(batch) + 1L, function(hit)
                         hit:nHits)))					#for a batch of hits, we make pairs of hits. The left-hand hit (V1) is from the batch, and the right-hand hit (V2) includes all other hits (including other matches)
    pairs[, c("q1", "s1", "q2", "s2") := data.table(hits[V1, .(qid, sid)], hits[V2, .(qid, sid)])]										#we retreive the ids of copies involved in the 2 hits (2 per clade)
    pairs[, c("inter1", "inter2", "intra1", "intra2") := data.table(hits[V1, pID], hits[V2, pID], pIDmat[cbind(q1, q2)], pIDmat[cbind(s1, s2)])]			#and retreive their pIDs within each clade (intra) of of the 2 hits (inter)
    pairs[, maxIntra := pmax(intra1, intra2)]
    cat("*")
    pairs[inter1 < maxIntra |
            inter2 < maxIntra, maxIntra, .(V1, V2)]			#returns pairs of hits where best intra-clade identity is higher than one inter-clade identity (that of hits)
  }
  
  pairs = rbindlist(lapply(hitBatches, criterion_1))						#does the above for batches of hit pairs and concatenates the results
  if (nrow(pairs) > 0) {
    #if there is not connected hits, there is not clustering to make
    cls = graph_from_data_frame(pairs, directed = F)
    cls = data.tableFromCommunities(cluster_fast_greedy(cls))			#we clustering with clauset et al. algorithm
    cat("-")															#to monitor progress
    com = integer(nrow(hits))
    com[cls$member] = cls$community			#we attribute hits to "communities"
    com[com == 0L] = 1:(sum(com == 0L)) + max(com)						#any hit that is not in a community now forms its own community
    hits[, comm :=  com]													#adds community column to the original table of htt hits
  } else
    hits[, comm := 1:.N]
  writeT(hits[, .(hit, comm, group)],
         stri_c("TEs/clustering/round1/", group, ".groups.txt"))		#saves results to disk for safety
  hits[, .(hit, comm, group)]
}



res = mclapply(groups,
               hitCommunities,
               mc.cores = nCPUs,
               mc.preschedule = F)
res = rbindlist(res)
if (length(superF) > 1)
  superF = "others"
writeT(res,
       stri_c("TEs/clustering/round1/", superF, ".allGroups.txt"))
print("finished")
