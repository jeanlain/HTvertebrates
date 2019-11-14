
#performs blastn search between TEs in many pairs of species of a given batch

library(data.table)
library(parallel)
args = commandArgs(trailingOnly=TRUE)	
job = as.integer(unlist(strsplit(args[1], ",")))		#the batch number
print(paste("job:", args[1]))
nCPUs = as.integer(args[2])								#number of searches to launch in parallel
sel = fread("pairsToBlastn.txt")
done= unique(list.files("TEs/blastn/done"))		#the output file names of jobs that are already done (to avoid redoing the search, in case the job is relaunched due to a bug, memory issue, and so on)

sel = sel[batch %in% job]			#selects the species pairs of the batch to process
sel[, out:=paste("TEs/blastn/out/", sp1,"-on-", sp2,".blastn.out", sep = "")]		#generates output file names
sel = sel[! basename(out) %chin% done,]												#discard species pairs if output files are present 
sel[, fas:=paste("TEs/copies/", sp1,".TEs.fasta.gz", sep ="")]						#query and db file names
sel[, db:=paste("TEs/blastn/db/", sp2,".TEs.fasta", sep = "")]
sel = sel[! out %chin% done,.(fas,db,out,size)]
setorder(sel, -size)														#we will launch the larger blast first
sel[,job := rep(1:(.N/nCPUs), nCPUs, length.out=.N)]
setorder(sel, job)
nt = max(1, as.integer(nCPUs/nrow(sel)))

blastn = function(fas, db, out) {
	system(paste("gunzip -c", fas, "| blastn -query - -db", db, "-max_target_seqs 1 -outfmt 6 -num_threads", nt, "| awk '{if ($3>=75 && $4>=300 && $12>=200) print $0}' >", out))		#we do not retain HSP with pID < 75 and score < 200 and alignement shorter than 300 pb. Filtering these hits from the begining was required, else we might have run out of disk space on our account
	file.rename(out, gsub("out/","done/", out, fixed = T))
	system(paste("gzip", out))			#compresses the output as they take a lot of space
	return(NULL)
}

res = sel[,mcMap(function(x,y,z) blastn(x, y,z), fas, db, out, mc.cores = nCPUs, mc.preschedule = F)]
print("finished")
