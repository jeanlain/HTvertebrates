
source("HTvFunctions.R")

#we shuffle species to generate random species pairs involved in transfers, for each super family. However, we may not generate "illegal pairs" for which no HTT can be infered, due to our requirement that species must be divergent by 2*120 My.
tree = read.tree("timetree.nwk")			#so we need to import the timetree
retainedHits = fread("gunzip -c data4-retained_hits.txt.gz")

args = commandArgs(trailingOnly=TRUE)	
nCPUs = as.integer(args[1])										#number of CPUs to use
node = ifelse(is.na(args[2]), which.max(node.depth(tree)), as.integer(args[2]))		#to permutate only species within a certain clade, denoted by the tree node number. Default to the basal node to permutate all species
n = ifelse(is.na(args[3]), 1000L, as.integer(args[3]))			#number of replicates we want, defaults to 1000
maxn = ceiling(n/nCPUs)		#the number of legal permutation we need per job launched in parallel

#as several species are involved in a hit group even though this would just represent one transfer, we only retain the two species for the best hit (best pID) per hit group. We also replace their names with tip numbers of the tree, as integer numbers speed things up and save memory considerably
setorder(retainedHits,-pID)				#puts best hits on top
HTThits = retainedHits[independent == T, .(sp1 = chmatch(sp1[1], tree$tip.label), sp2 = chmatch(sp2[1], tree$tip.label), repl = 0L), by = .(hitgroup, superfamily)]		#repl is a replicate number (= 0 for original, non-permuted species)

nTr = HTThits[,.N, by = superfamily]		#number of independent transfers per super family, as we combine super families (within classes) that are involved in less than 20 transfers
nTr[,class := ifelse(grepl("CMC|hAT|Mariner|Maverick|Merlin|PIF|PiggyBac", superfamily),"DNA","RNA")]
nTr[, combined := ifelse(N < 20, paste("other", class, sep = " "), superfamily)]

HTThits[,superfamily := nTr[match(HTThits$superfamily, superfamily), combined]]		#replaces small super families with the combined names
HTThits[,occ := occurences(superfamily)]				#as there may be too many hits per superfamily to obtain only "legal" permutations, we will split certain super families in batches of hits (when they encompass more than 120 transfers)

nTr = nTr[,.(N =sum(N)), by=combined]
nTr[,maxi := N / ceiling(N/121)]					#the max number of transfers we allow per superfamily batch
HTThits[,maxi := nTr[match(HTThits$superfamily, combined), maxi]]	#which we transfer to the hits table

HTThits[superfamily == "Tc1/Mariner", maxi := 61L]		#for mariners, we need even smaller batches of ≤ 61 transfers (for some reason, it is harder to obtain legal permutations in these TEs)
HTThits[,batch := ceiling(occ / maxi)]		
hitList = split(HTThits[,.(sp1, sp2, repl)], f= list(HTThits$superfamily, HTThits$batch), drop = T)		#splits the hits by superfamily and batch

l = length(tree$tip.label)						#the total number of species
divTimeMat = cophenetic(tree)					#matrix of divergence times
tooClose = divTimeMat < 240						#true for species pairs that are "too close" to be involved in HTT (a species here is a row/column index of the logical matrix)
toShuffle = tipsForNode(tree, node)				#integer numbers of species we will permutate (tip numbers of the tree)
newsp = matrix(rep(1:l, 10^5), ncol = 10^5)		#this matrix will contain the permuted "replacement species" during the work, one set of species (permutation) per column. Note that we anticipate 10^5 permutations per batch, although we keep fewer. This is because many permutations may be "illegal". The row indices of this matrix = species (integer ids)


shuffle = function(hits, superfamily) {		#permutates species for htts (hits) of a superfamily
	print(superfamily)						#to print progress, which is the only use of the superfamily argument
	nr = nrow(hits)
	#to speeds things up, we generate 10^5 permutations in a row as the vast majority lead to invalid transfers. This allows taking advantage of vectorisaton after that. 
	sp1 = rep(0:(10^5-1), each = nr)*l + hits$sp1		#to do this, we replicate the hits 10^5 times, but incrementing species ids at each replicate (by the total number of species = l)
	sp2 = rep(0:(10^5-1), each = nr)*l + hits$sp2
	shuffleWork = function(job)	{					#performs the permutations in parallel. job is a simple integer identifier
		legal = NULL								#will be a matrix of "legal permutations" (same format as the newsp matrix, where row numbers = species ids)
		g = 0										#indicator to tell when to stop
		repeat {									#until we obtain maxn legal permutations
			newsp[toShuffle,] = replicate(10^5, sample(toShuffle))			#the workhose intruction that performs the permutation (10^5 at a time). This creates a matrix of permuted species, with 10^5 columns (each is a vector of shuffled species)
			newPairs = cbind(newsp[sp1], newsp[sp2])						#replacing species by the sampled ones in the transfers
			illegal = colSums(matrix(tooClose[newPairs], nrow = nr))		#number of illegal transfers per permutation (we use a matrix to quickly count them via colSums)
			if(any(illegal == 0L)) {
				legal = cbind(legal,newsp[,illegal == 0])					#we extract the columns corresponding to permutations with no illegal transfer and add them to the retained permutations
				g = ncol(legal)
				cat(".")													#progress indicator
			}
			if(g >= maxn) break
		}
		newSp = as.vector(legal[,1: maxn])			#"unrolls" the matrix of legal permutations (useful to replace species with the sampled ones, as the position in this vector will be the original species identifier, and its value is the replacement species). We do not retain more than maxn legal permutations (there may actually be up to 10^5 if there are few hits).
		repl = rep(1:maxn, each = nr)				#generates a permutation id number
		pairs = data.table(sp1 = rep(hits$sp1, maxn) + (repl-1L)*l, sp2 = rep(hits$sp2, maxn) + (repl-1)*l, repl)	#we replicate the transfers n times but incrementing species numbers, so that we can easly perform the replacements with random species. This is similar to what we did to create the newPairs matrix, except this time we use a data table
		pairs[,c("sp1","sp2","repl") := .(newSp[sp1], newSp[sp2], repl+job*maxn)]		#we use the job id to generate final replicate identifiers, which will have to differ between jobs
		pairs
	}
	pairs = mclapply(1:nCPUs-1L, shuffleWork, mc.cores = nCPUs, mc.preschedule = F)	#does the permutations in parallel 
	rbindlist(pairs)
}

randomHTTs = Map(shuffle, hitList, names(hitList))		#applies permutations to superfamilies successively

for (i in 1:length(randomHTTs)) {			#We generate unique repliate identifiers within eahch super family
	nrep = nrow(randomHTTs[[i]]) / nrow(hitList[[i]])		#the number of replicates is the number of simulated HTTs in the batch divided by the number of observed HTTs
	randomHTTs[[i]][,repl := rep(1:nrep, each = nrow(hitList[[i]]))]
	randomHTTs[[i]] = rbind(hitList[[i]], randomHTTs[[i]])			#we also concatenates real and random HTTs (we can differentiate the two because repl == 0 for real transfers)
}

saveRDS(randomHTTs, file = stri_c("allPermutations_Node.", node, ".RDS"))