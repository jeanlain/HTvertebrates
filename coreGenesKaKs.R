#computes Ks between species at BUSCO CDS

source("HTvFunctions.R")

dir.create("Ks/KaKs")
allAA = readAAStringSet(list.files("CDS", pattern = ".aa.fas", full.names = T))                 #we import all prot and CDS sequences at once
allCDS = readDNAStringSet(list.files("CDS", pattern = ".CDS.fas", full.names = T))
blastpFiles = list.files("Ks/blastp", pattern = ".out", recursive = T, full.names = T)             #results files of diamond blasp searches (tabular)
sp = unique(extractSpeciesNames(names(allAA)))                                                
pairs = allPairs(sp, reciprocal = F)                                            #all pairs of species, no reciprocity
allForward = stri_c("Ks/blastp/", pairs$V1,"/", pairs$V2,".on.",pairs$V1,".out")             #all possible forward searches
allReverse =  stri_c("Ks/blastp/", pairs$V2,"/",pairs$V1,".on.",pairs$V2,".out")            #all reverse searches
pairs = data.table(forward = intersect(allForward, blastpFiles), reverse = intersect(allReverse, blastpFiles))          #all pairs of forward/reverse searches that were undertaken
pairs[,batch :=rep(1:(.N/50), length.out =.N)]                                                #a batch will process approximately 50 species pairs. Big batches may cause segfault during alignments.
pairs = split(pairs, pairs$batch)
m = mclapply(pairs, compute_Ks, mc.cores = 30, mc.preschedule = F)
