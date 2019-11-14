### prepares and launches the pair-wise reciprocal blastn searches of TEs (one genome against another)

source("HTvFunctions.R")

# we generate pairs of species for which we will blast TEs, since we will not search for HTT between species that are too closely related
tree = read.tree("timetree.nwk")
sp = tree$tip.label			
distMat = cophenetic(tree)			#the divergence time between all species (as a matrix). We use it to avoid comparing species that are too closely related, and between which we won't search for HTT
pairs = data.table(allPairs(sp, reciprocal =T, same = F))		#a fonction that generates pairs (faster than combn())
setnames(pairs, c("sp1","sp2"))
pairs[,divTime := distMat[cbind(sp1, sp2)]]
selectedSpeciesPairs = pairs[divTime > 80]			#only retained species that have diverged in the last 40 My

#we will launch the longer blast searches first (those involving the bigger fasta files), to better use the CPUs
selectedSpeciesPairs[,size := file.size(stri_c("TEs/copies/", sp1,".TEs.fasta.gz")) + file.size(stri_c("TEs/copies/", sp2,".TEs.fasta.gz"))]
setorder(selectedSpeciesPairs, -size)				
f = selectedSpeciesPairs[,sp1=="Ambystoma_mexicanum" | sp2 =="Ambystoma_mexicanum"]		#these pairs will constitute a dedicated batch (Ambystoma_mexicanum has much more TEs than the others)
selectedSpeciesPairs[f, batch := 0L]
selectedSpeciesPairs[!f,batch := rep(1:(.N/4000), length.out=.N)]		#a single batch contains approximately 4000 searches. This was determined given the number of CPUs per node on the server (30 cores) and the maximal duration of jobs (4 days). Using an array of jobs via sbatch --array would have been better, but I didn't know how to do it by then)
writeT(selectedSpeciesPairs, "pairsToBlastn.txt")

#makes blastn database for blast between TE copies of two species
library(parallel)
files = list.files("TEs/copies", pattern = ".gz", full.names =T)		#list compressed fasta files of TE copies generated in step 2-TEextractionAndComposition.R
dir.create("TEs/blastn/db", recursive = T) ; dir.create("TEs/blastn/done") ; dir.create("TEs/blastn/out")

db = gsub("copies", "blastn/db", gsub(".gz","",files))
f = !file.exists(paste(db, ".nin", sep=""))
makedb = function(fas, db)  {
	system(paste("gunzip -c", fas, "| makeblastdb -in - -dbtype nucl -title", db, "-out", db), ignore.stderr = T)
	return(NULL)
}

res = mcMap(function(x,y) makedb(x, y), files[f], db[f], mc.cores = 10, mc.preschedule = F)


###### we launch the blastn on the selected species pairs (several jobs sent for different batches)
system('sbatch --mail-type=BEGIN,END,FAIL --cpus-per-task=20 --mem=60G --wrap="Rscript pairwiseSpeciesBlastn.R 1 20"')		#this is the command for batch 1 with 20 CPUs

