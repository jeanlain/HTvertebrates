

##%######################################################%##
#                                                          #
####       This stage prepares the clustering of        ####
####            TE hits to count HTT events.            ####
#                                                          #
##%######################################################%##


source("HTvFunctions.R")
dir.create("TEs/clustering")
httHits = fread("occ2000ds05.txt")		#selected TE-TE hits from the previous step

###### there are too many hits for this clutering, we only retain 200 of hits
###### per single-linkage cluster, as hits from the same cluster are very likely
###### to represent the same HTT

setorder(httHits,-length.aa)				#we favor hits that involve the longest protein coding regions
httHits = httHits[occ <= 200L, ]
writeT(httHits, "occ200ds05.txt")			#and write them to file

# extracting TE copies sequences in the retained hits, as we will compare copies
# with super families with blastn all vs all (referred to as "self blast"
# afterwards). This will be used to determine of two hits represent the same HTT
# (identity at the hits is lower than identity between copies in each clade,
# refered to as "criterion 1" in the paper), and also for the anlysis of TE
# evolution wihtin genomes
copies = httHits[, union(
  stri_c(sp1, ":", query, ":", f1, "#", superF),
  stri_c(sp2, ":", subject, ":", f2, "#", superF)
)]  #we add species and classification information to copy names, as specified in the fasta file of TE copies, so we can extract them
write(copies, "temp.bed")		#temporary file used by seqtk
seqtk(
  "TEKs/selectedCopiesAA300cl2000.fas",
  "temp",
  "TEs/clustering/selectedCopiesKs05occ200.fas"
)
file.remove("temp.bed")

############### we prepare self blastn within each super family of the copies extracted above
## we first rename copies with integers, to reduce the size of blastn files and help with various things

seqs = readDNAStringSet("TEs/clustering/selectedCopiesKs05occ200.fas")
copies = names(seqs)

superF = gsub("/", ".", stri_extract_last(copies, regex = "[^#]+"), fixed = T)	#extracts super family names and replaces slashes with periods as file names will containes super family names
copyID = data.table(copy = copies, id = 1:length(copies))						#attributes integer number to copies, which we will used instead of long copy names
writeT(copyID, "TEs/clustering/selectedCopiesKs05occ200.IDs.txt")				#writes this correspondence to disk as it will be reused many times afterwards
names(seqs) = copyID[chmatch(names(seqs), copy), id] 							#replaces copy names with the integers in the sequecnes

dir.create("TEs/clustering/blastn/db", recursive = T)
all.db = split(seqs, superF)													#splits copy sequences by super families, as blastn searches will be done within super families
fasNames = stri_c("TEs/clustering/blastn/db/", names(all.db), ".fas")									#prefixes of blast db files
m = Map(writeXStringSet, all.db, fasNames)
m = lapply(fasNames, function(x)
  system(paste("makeblastdb -dbtype nucl -in", x)))		#makes blast databases

#splitting query files for better blastn parallelisation (num_threads of blastn isn't very efficient so we leave it to 1) and better management of jobs + outputs (which are very big)
dir.create("TEs/clustering/blastn/split")
dir.create("TEs/clustering/blastn/out")
dir.create("TEs/clustering/blastn/done")

nD = as.numeric(sapply(all.db, length))													#number of sequences per super family, which we used to decide how much we split the queries
nComp = data.table(query = names(all.db), nc = nD * as.numeric(sapply(all.db, length)))	#the criterion is the size of query * size of db
nComp[, nt := round(nc / (max(nc) / 200))]													#this ensures that a query is not split in more than 200 parts
nComp[nt == 0, nt := 1L]

indices = unlist(Map(function(x, y)
  rep(1:x, length.out = y), nComp$nt, nD))			#creates an integer index corresponding to the sub-query (max 200), for each sequences
queryParts = stri_c("TEs/clustering/blastn/split/",
                    rep(nComp$query, nD),
                    "_",
                    indices,
                    ".fas") 				#which is used to attribute each sequence to a sub-query (part)
spl = split(unlist(all.db, use.names = F), queryParts)									#so we split the copy sequences in these sub queries
l = Map(writeXStringSet, spl, names(spl))												#saves them to fastas

dbs = gsub("split",
           "db",
           stri_extract_last(names(spl), regex = "[A-Z]+.[^_]+"))		#the databases corresponding to the split queries
searches = data.table(query = names(spl), db = stri_c(dbs, ".fas"))	#creates a table listing the blast searches to do
searches[, out := gsub("split", "out", query)]
searches[, out := gsub(".fas", ".out", out)]
searches[, batch := rep(1:(.N / 30), length.out = .N)]				#blast searches from a given "batch" will be launch from a given R instance, on a node (using sarray would probably have been better)
writeT(searches, "copiesToSelfBlast.txt")

######## script for the self-blastn searches. This is done by batches. These use a lot of ram and take quite a bit of time, as we don't limite the number of hits
system(
  'sbatch --mail-type=BEGIN,END,FAIL --cpus-per-task=4 --mem=120G --wrap="Rscript TEselfBlastn.R 1 4"'
)  #example for batch 1
