## self blasts TE copies of super family against themselves (launched in step 8)
library(data.table)
library(parallel)
args = commandArgs(trailingOnly=TRUE)	
job = as.integer(unlist(strsplit(args[1], ",")))
nCPUs = as.integer(args[2])

print(paste("job:", job))
searches = fread("copiesToSelfBlast.txt")

done= list.files("TEs/clustering/blastn/done")

f = searches[,batch %in% job & !basename(out) %in% basename(done)]		#in case jobs need to be relaunched, we avoid redoing already completed searches
nt = max(1, as.integer(nCPUs/sum(f)))				#number of threads per blast search (usually 1, but if only a very few searches are done, it may be higher)

blastn = function(query, db, out) {
	done = gsub("/out/","/done/", out, fixed =T)
	system(paste("blastn -query", query, "-db", db, "-max_target_seqs 100000 -max_hsps 1 -outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore' -num_threads", nt,"| awk '{if ($4>=100 && $1<$2) print $0}' >", out))		#because this is a self blast, we select results where the query name (an integer) is < subject name, to reduce output file size by two (we only retreive the hit copy1 vs copy2, not the reverse). 
	#note also the max_target_seqs argument, which is larger than the number of sequences in the db and ensures that all hits are returned
	file.rename(out, done)
	return(NULL)
}

res = searches[f,mcMap(blastn, query, db, out, mc.cores = min(nCPUs, sum(f)), mc.preschedule = F)]
print("finished")
