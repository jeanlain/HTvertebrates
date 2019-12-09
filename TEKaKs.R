#performs Ka/Ks computations between homologous TEs

source("HTvFunctions.R")
library(seqinr)
args = commandArgs(trailingOnly=TRUE)	
TEhitFile = args[1]							#the TE-TE hits, for which we will compute Ks
blastxFile = args[2]						#the HSPs between these TEs and protein regions, to guide our Ks computations
seqFasta = args[3]							#fasta of copy sequences
outputFolder = args[4]
nCPUs = as.integer(args[5])

TEhits = fread(TEhitFile)					
blastx = fread(blastxFile)			

covPerCopy = blastx[,sum(end-start+1L), by = copy]															#selects HSPs for which there is enough protein regions per copy
TEhits[,hit := 1:.N]																						#to keep track of hits, we assign them an integer identifier in a column
hits = TEhits[, which(query %chin% covPerCopy[V1>30L, copy] & subject %chin% covPerCopy[V1>30L, copy])]		#hits are simply row indices of TEhits table
nJobs = round(length(hits) / 2500)								#nb of jobs to run. The size of a batch must not be too big as this may cause segfaut issues with the pairwise alignments).
hits = splitEqual(sample(hits), n = nJobs)																	#splits the hits in batches for parallel computations. Note that we randomize hits (rows) to reduce differences in job durations (the longer hits are on top of the TEhits table)

seqs = readDNAStringSet(seqFasta)												
names(seqs) = copyName(names(seqs))																		#to match copy names of the TEhits table
qSeqs = TEhits[,subSeq(seqs[query],qStart,qEnd)]
sSeqs = TEhits[,subSeq(seqs[subject],sStart,sEnd)]

qSeqs = lapply(hits, function(x) qSeqs[x])
sSeqs = lapply(hits, function(x) sSeqs[x])
TEhits = lapply(hits, function(x) TEhits[x, .(query, subject, qStart, qEnd, sStart, sEnd)])
rm(seqs, covPerCopy) ; gc()

TEKs = function(job) {																						#this computes Ks on a batch of hits. job is the job number
	aln = alignWithEndGaps(qSeqs[[job]], sSeqs[[job]])														#pairwise aligns copies
	TEhit = TEhits[[job]]
	nuc = splitAlignment(aln, TEhit[,.(query, subject, qStart, qEnd, sStart, sEnd)])						#splits alignments into individual nucleotides and positions
	nuc[, hit1 := assignToRegion(blastx[,.(copy, start, end)], .(seq1, pos1))]								#determines the blastx HSP (with coverage 1) covering each aligned position. Note that positions outside these HSPs would have NA values
	nuc[, hit2 := assignToRegion(blastx[,.(copy, start, end)], .(seq2, pos2))]								#we do it for both copies
	nuc = nuc[!is.na(hit1) & !is.na(hit2)]																	#we can already remove positions outside protein HSPs
	nuc[, c("rev1", "rev2") := .(blastx[hit1, rev], blastx[hit2, rev])]												#we also retreive the direction of the blast HSPs
	nuc[, frame1 := ifelse(!rev1, (pos1 - blastx[hit1, protStart]) %% 3L + 1L, -((blastx[hit1, protStart] - pos1) %% 3L + 1L))]		#so we can determine the position in codon ="frame" (1 to 3 or -1 to -3) of each position
	nuc[, frame2 := ifelse(!rev2, (pos2 - blastx[hit2, protStart]) %% 3L + 1L, -((blastx[hit2, protStart] - pos2) %% 3L + 1L))]
	nuc[TEhit[aln, qStart > qEnd], frame1 := -frame1]							#reverses the frame if the copy is itself reversed in the megablast hit
	nuc[TEhit[aln, sStart > sEnd], frame2 := -frame2] 							
	
	nuc = nuc[frame1 == frame2]																					#removes aligned position where the frames don't match
	nuc[,c("hit1","hit2","rev1","rev2", "frame2") := NULL]														#as well as columns that are no longer necessary
	
	#we now compute Ka and Ks values with seqinr. for this we first discard incomplete codons from alignments, as we will concatenate nucleotides back into sequences. Codons will be identified by integers
	#we determine whenever a new codon starts (jumps in position within codon, new alignment, indel)
	n = 2:nrow(nuc)
	shift = nuc[, c(T, (abs(frame1[n]-frame1[n-1]) != 1L | aln[n] != aln[n-1L] | abs(pos1[n]-pos1[n-1L]) != 1L | abs(pos2[n]-pos2[n-1L]) != 1L))]
	codon = cumsum(shift)																					#generates an integer for each codon
	nuc = nuc[codon %in% which(tabulate(codon) == 3L)]														#we retain only complete codons (found 3 times)
	nuc = nuc[aln %in% which(tabulate(aln)> 30L)]															#we discard alignments that are too short

	flat = nuc[,.(query = stri_flatten(base1), subject = stri_flatten(base2), forward = frame1[1] > 0L), by = aln]	#we now concatenate the positions backs into strings. We also note if the copies are in reverse in respect to the protein (frame is negative)
	flat[forward == F, c("query","subject") := .(revCom(query), revCom(subject))]										#we reverse complement the sequence of these copies
	spl = split(flat[,cbind(query, subject)], 1:nrow(flat))																#splits the data table into individual pairs of strings
	alns = lapply(spl, seqinrAlignment)																					#so we can produce seqinr alignment objects
	KaKs = lapply(alns, kaks)
	dt = rbindlist(lapply(KaKs, as.data.table))															#concatenate the result list into a data.table
	dt[,c("hit","length") := .(hits[[job]][flat$aln], nchar(flat$query))]
	writeT(dt, stri_c(outputFolder,"/KaKs.", job,".txt"))							#write results for the finished job for safety (we write the combined results at the end)
	dt
}	

res = mclapply(1:length(TEhits), TEKs, mc.cores = 10, mc.preschedule  = F)
writeT(rbindlist(res), stri_c(outputFolder,"allKaKs.txt"))
