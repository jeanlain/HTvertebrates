
# We evaluate the reliability of hit groups. A reliable hit group (putative HTT) must involve enough TE copies per clade (to reduce the risk that it actually reflects contamination). But since we filtered a lot of hits (hence copies) in previous steps, we evaluate whether copies involved in HTT have close relatives in their respective genomes in order to retreive similar copies. So we blast copies against the TEs of their host genomes.

source("HTvFunctions.R")
httHits = fread("oc200HitGroup.txt")				
fas = readDNAStringSet("TEs/clustering/selectedCopiesKs05occ200.fas")			#copies involved in the putative HTTs (hits)
sp = extractSpeciesNames(names(fas))											#we will save fasta files for seperate species, to make one blast per sp. So we split the sequences
spl = split(fas, sp)		
dir.create("TEs/clustering/testConta")													
fasFiles = stri_c("TEs/clustering/testConta/", names(spl), ".selectedCopies.fas")
l = Map(writeXStringSet, spl, fasFiles)

system('sbatch --mail-type=BEGIN,END,FAIL --cpus-per-task=15 --mem=100G --wrap="Rscript blastAgainstAllTEs 15"')

#collecting results
blast = fread("TEs/clustering/testConta/blastn/all.out")		#concatenated results of the intraspecific blast performed above, where the query is a copy involved the candidate HTT hits and the subject is a TEs of the same species (and may not be in the copies aming HTT hits, which is the whole point of the blast)
blast = blast[V3 > 75 & V4 >= 200L]						#removes hits of low quality since none of the retained HTT hits had pID < 75%
blast[,sp := splitToColumns(V1, ":",1)]					#extracts species names. 
blast[,c("V1","V2") := .(copyName(V1), copyName(V2))]
copyIDs = fread("TEs/clustering/selectedCopiesKs05occ200.IDs.txt")				#we retreive copy integer IDs
copyIDs[,copy := copyName(copy)]

spForCopy = blast[,unique(c(V1,V2)), by = sp]			#makes correspondance between copy and species
blast = blast[V1 != V2]									#removes self hits
blast[,q := copyIDs[chmatch(V1, copy), id]]

spl = reList(split(1:nrow(blast), blast$q))						#this list all hits (by their row number) involving a given copy
copiesInHitGroups = httHits[,unique(c(q, s)), by = hitGroup]										#removes redundancies wihtin hit groups	
copiesInHitGroups[,name := copyIDs[match(V1, id), copy]]											#retreives full copy names (useful later)
hitsInHitGroups = copiesInHitGroups[,.(hit = unlist(spl[V1])), by = hitGroup]								#so we get all intraspecific hits involving each copy of each hit group
hitsInHitGroups[,c("query","subject", "pID", "score") := blast[hit, .(V1, V2, V3, V12)]]			#add useful hit information in new columns
setorder(hitsInHitGroups, hitGroup, -pID)
#hitsInHitGroups = hitsInHitGroups[!duplicated(data.table(hitGroup, subject)),]								#we now decide if a subject hitted by copies can be assigned to the hit group. 
bestIDs = httHits[,.(mean(pID), mean(score)), by = .(hitGroup, copy=query)]						#we will compare the similarity of each hit between a subject and copy to that of hits involving the copy within the hit group (those hits are supposed to reflect HTT, while the former hit should reflect within-genome transposition). Se we retreive the mean pID and scores per copy per hit group
bestIDs = rbind(bestIDs ,httHits[,.(mean(pID), mean(score)), by = .(hitGroup, copy=subject)])		#we need to do that for copies that are in "query" and "subject" columns
bestIDs[,pair := stri_c(hitGroup, "_", copy)]														#hit group-copy pair indentifier
hitsInHitGroups[,c("meanID", "meanScore") := bestIDs[chmatch(stri_c(hitsInHitGroups$hitGroup, "_", query),pair),.(V1, V2)]]		#so we retreive the mean pID and score for hits involving a query, for each hit group

perSubject = hitsInHitGroups[,mean(pID > meanID), by = .(hitGroup, subject)]					#computes the poportion of copies of a hit group that are more similar to a subject than they are to TEs in the other clade (supposed HTT)
perSubject = perSubject[V1 > 0.5]													#and we retain subjects that are on average more similar to same-genome copies than these are to TEs of the other clade

copiesInHitGroups = rbind(copiesInHitGroups[,.(hitGroup, copy = name, new = F)], perSubject[,.(hitGroup, copy = subject, new = T)])	#we can now concatenate these "new" copies with the "old" ones
copiesInHitGroups = copiesInHitGroups[!duplicated(data.table(hitGroup, copy))]				#removes "new" copies that were in fact already present in hit groups

copiesInHitGroups[,sp := spForCopy[match(copy, V1), sp]]			#we also retreive the species to which copies belong, as we will need to know their clade (our contamination filter imposes a certain number of copies per clade per hit group)
sp1perHitGroup = httHits[,unique(sp1), by = hitGroup]				#to determine the clade to which the copy belongs, we make a list of species of the "left" clade in each hit group
sp1perHitGroup = split(sp1perHitGroup$V1, sp1perHitGroup$hitGroup)
setorder(copiesInHitGroups, hitGroup)
copiesInHitGroups$inClade1 = copiesInHitGroups[,sp %chin% sp1perHitGroup[[hitGroup]], by = hitGroup]$V1		#this new column should be TRUE for copies belonging to the "left" clade (sc1) of a hit group

############################# excluding hit groups that may result from contamination or VT based on additionnal filters
#to exlclude VT, we evaluate the distribution of Ks of TEs in a hit group, in lights ot the thresholds we imposed (the 0.5% quantile of the BUSCO Ks distribution, or 0.5)
Ks = fread("gunzip -c Ks200AAnoRedundancy.txt.gz")					#file of filtered BUSCO Ks (200 AA alignments and one score per BUSCO gene per pair of clades)
quant = Ks[Ks < 10,.(min = min(Ks), q05 = quantile(Ks, 0.005), m = mean(Ks), divTime = divTime[1]), by = clade] #collects useful information, including the 0.5% Ks quantile that we used as a threshold

tree = read.tree("timetree.nwk")
setorder(httHits, hitGroup, -pID)   #so that hitGroupStats below is sorted by hit group
# generate statistiques per hit groups. Mainly, the Ks mode (see function ksMode in HTvFunctions.R) and number of copies per clade
hitGroupStats = httHits[, data.frame(maxkS = max(ks), ksMode(ks + 2*vks), N=.N, nq = length(unique(q)), ns = length(unique(s))), by = .(hitGroup, mrca, superF)]  #nq and ns are number of copies per hit group (and each of the two clades) actually involved in the hits
hitGroupStats = merge(hitGroupStats, copiesInHitGroups[,.(nQ = sum(inClade1), nS = sum(!inClade1)), by = hitGroup], by = "hitGroup")	#nQ and nS = the above + copies retrieved with our blast above
hitGroupStats[,c("qKs", "mKs",  "age") := data.table(quant[match(mrca, clade), .(q05, m)], nodeDepth(tree, mrca))]		#adds 0.5% quantile of busco Ks of the corresponding clade and other useful information

#checks if ka/ks of TEs correlates with the divergence of clades involved. It should not (a priori), especially for DNA transposons, unless there is some spurious VT for which ka/ks should be closer to 1, while ka/ks for TEs transferred by HTT should be < 1 due to selection during HTT
ageBreaks = hitGroupStats[,seq(min(age),max(age), length.out = 20)]				#we will make class of divergence times for better visualisation
httHits[,age :=  nodeDepth(tree, mrca)]											#the divergence time of the species involved in each hit 
httHits[,ageClass := .bincode(age, ageBreaks, include.lowest = T)]				#assigns this time to the classes of divergence times
perAge = httHits[ka/ks <3, .(kaks = weighted.mean(ka/ks, length.aa), age = mean(age)), by = .(DNA = grepl("DNA", superF), ageClass)]		#obtains mean ka/ks per divergence time calss and TE class (DNA or RNA)
p = perAge[DNA ==T, plot(age, kaks, ylim = c(0, max(kaks)), xlab = "Time since divergence (My)", ylab = "Ka/Ks")]		#makes the plot for DNA transposons (there was no obvious trends for class I TEs) = figure S1 of the paper
abline(v = 120, col ="red")

hitGroupStats[,keep := nQ >=5L & nS >=5L & nq >= 2L & ns >=2L & (nMode > nLast +20L | maxkS < qKs - 0.2) & age > 120]	#based on the above, we retain hit groups involving at least 5 copies in each clade (to account for possible contamination), including 2 copies per clade in the retained hits, and whose kS is low enough and involving lineages thare are quite divergent
writeT(hitGroupStats, "TEs/clustering/hitGroupStats.txt")
