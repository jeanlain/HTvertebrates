#lauches diamond blast search  between pairs of species, the list of which is read from a file ("pairsToBlastp.txt")

library(data.table)
library(parallel)
args = commandArgs(trailingOnly=TRUE)	
job = as.integer(unlist(strsplit(args[1], ",")))
print(paste("job:", job))
nCPUs = as.integer(args[2])
pairs = fread("pairsToBlastp.txt", header = T)
pairs = pairs[batch %in% job]

blastp = function(db, aa, out, daa) {
	system(paste("diamond blastp --quiet --sensitive -p 1 -k 1 --max-hsps 1 -d", db, "-q", aa, "-a", daa))
	system(paste("diamond view -a", daa, "-o", out))		#converting output file to tabular format
	file.remove(daa)
	print(paste(out))
}

res = pairs[!file.exists(out), mcMap(blastp, db, aa, out, daa, mc.cores = nCPUs, mc.preschedule = F)]
