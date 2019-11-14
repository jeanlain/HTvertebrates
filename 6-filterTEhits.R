
####### Filters the TE-TE hits
####### As there are too many hits to process in a single table, we first filter-out all hits whose pID (not Ks) is lower than the 1 - the 0.5% quantile of Ks of the corresponding clade pair. Hits not passing this flter would have been very unlikey to pass downsteam filters based on TE Ks

source("HTvFunctions.R")
Ks = fread("Ks200AAnoRedundancy.txt")						#file of core gene Ks generated in step 5-coreGeneKs.R
quant = Ks[, .(q = quantile(Ks, 0.005)), by = clade]		#computes the 0.5% Ks quantile for every clade
tree = read.tree("timetree.nwk")
mrcaMat = mrca(tree)										#matrix of mrcas for each pair of species

#we add this quantile to for every species pair that was blasted
pairs = fread("pairsToBlastn.txt", select = c(1,2, 8))		#pairs of species whose TEs were compared by blast, generated in 4-blastTEs.R
pairs[, q := 100 - 100 * quant[match(mrcaMat[cbind(sp1, sp2)], clade), q]]						#adds a column indicating the minimum pID to retain blast hits
writeT(pairs, "pairsToFilter.txt")

#the script below also remove any hits involving a TE that belongs to a family whose consensus may involve a non-TE gene (see step 3-findDubiousTEs.R)
system('sbatch --cpus-per-task=10 --mem=30G --wrap="Rscript filterTEblastnHits.R 10"')	

#importing filtered hits that have been collapsed in one tabular file (which is quite big)
TEhits =fread("all.quantile005score200.blastn.out", header = F, sep ="\t", col.names = c("query","subject","sp1","sp2","f1","superF1", "f2","superF2","pID","length","qStart","qEnd","sStart","sEnd","score"))		
#exporting TE copies involved in hits to blast against rep proteins
copies = TEhits[,data.table(sp = c(sp1, sp2), seq = c(query, subject), fam = c(f1, f2), superF = c(superF1, superF2))]		#copies are identified by species, sequence name, family and superfamily names, for both copies and subjects of the hits
copies = copies[!duplicated(seq)]		#copies are typically involved in several hits, so we remove duplicates
pasted = copies[,stri_c(sp, ":", seq, ":", fam, "#", superF)]		#we paste the field and separate them with colons, except the super family name
spl = split(pasted, copies$sp)

dir.create("TEs/selectedCopies")
m = Map(function(x, y) write(x, file = y), spl, stri_c("TEs/selectedCopies/", names(spl), ".txt"))		#copies names are saved in separate files for each species, as the blastx will be run in parallel

#### blasting copies against the diamond repeatPep database, as we use this to compute Ks between TEs

system('sbatch --cpus-per-task=10 --mem=50G --wrap="Rscript successiveBlastx.R 10"')

######## using the blastx information to select hits involving protein regions

blastx = fread("TEs/blastx/all.copies.successiveBlastx.out", header = T, sep ="\t", col.names = c("query","subject","pID","length","mismatch", "gapopen","qStart","qEnd","sStart","sEnd","evalue","score","ex"))				#the combined blastx results generated above

blastx[,copy := copyName(query)]									#generating copy names from query names in the blastx file (i.e., excluding species and family info), to match the TE-TE hits

perCopy = blastx[ex == 2L,.(start = min(c(qStart, qEnd)), end = max(c(qStart, qEnd))), by = copy]			#computing the first and last position of all aligments on "expected" proteins (ex==2) for each copy, which we infer as the protein region of each copy (even if it may encompasse several ORFs)

#below we select TE-TE hits involving TE protein regions  >= 300 bp, to avoid computing unnecessary Ks (since we require kS computed on ≥300 bp)

TEhits[,c("qSt", "qEn", "sSt","sEn") := data.table(fixRanges(cbind(qStart, qEnd)), fixRanges(cbind(sStart, sEnd)))]			# = blastn coordinates where starts always < ends
TEhits[,c("qS","qE") := perCopy[chmatch(query, copy),.(start, end)]]									#we add portein ranges for query and subject in blastn hits
TEhits[,c("sS","sE") := perCopy[chmatch(subject, copy),.(start, end)]]
TEhits[is.na(qS), c("qS","qE") := 0L]																	#prevents bugs that are caused by NAs in interesection()
TEhits[is.na(sS), c("sS","sE") := 0L]	
TEhits[,c("qS", "qE") := data.table(intersection(qSt, qEn, qS, qE, T))]									#the region of the query that aligns on the subjects and also aligns on proteins
TEhits[,c("sS", "sE") := data.table(intersection(sSt, sEn, sS, sE, T))]									#same for the subject
TEhits[,c("qsS","qsE") := .(qSt + sS - sSt, qEn + sE - sEn)]											#converts the above into query coordinates
TEhits[, inter := pmin(qE, qsE) - pmax(qS, qsS) + 1L]													#so we can compute the lenght of the TE alignment that also aligns on proteins
TEhits[is.na(inter) | inter < 0L, inter := 0L]	
TEhits[,c("qSt","qEn","sSt","sEn","qS","qE","sS","sE","qsS","qsE") := NULL]								#removes columns we no longer need

TEhitsCoveringProteins = TEhits[superF1 == superF2 & inter >=300L]										#selects hits with sufficient protein region (300 bp) and between TEs assigned to the same super families
setnames(TEhitsCoveringProteins, "superF1", "superF") ; TEhitsCoveringProteins[,superF2 := NULL]		#so now we only need one superfamily name per hit

writeT(TEhitsCoveringProteins, "all300aaSameSuperF.txt")

#there are still too many hits and copies to compute Ks. Selecting among redundant hits that may represent the same HTT. This is based on the fact that many hits share a copy (see Method section).
setorder(TEhitsCoveringProteins, -inter, -score)										#puts bet hits on top (longest protein region then higher score)
TEhitsCoveringProteins[, cl := clusterFromPairs(query, subject)]						#single-linkage clusters of copies based on hits
TEhitsCoveringProteins[,occ := occurences(data.table(cl, sp1, sp2))]					#for each hit, the number of time we have seen a cluster in a species pair. Within a species pair, hits of the same cluster would represent the same transfer and are thus redundant. 

sel = TEhitsCoveringProteins[occ <= 2000L,]		#we select at most 2000 hits between copies from a given cluster for any species pair
sel[,mrca := mrcaMat[cbind(sp1, sp2)]]			#adds the mrca of each species pair, for later use
writeT(sel, "allOCC2000.txt")
