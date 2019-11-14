######  finds TE consensuses (families) that may cover non-TE genes
######	some steps (diamond searches) are not in this script. I (Jean Peccoud) did not perform them

source('HTvFunctions.R')

workFolder = ("TEs/findDubious/")
dir.create(workFolder)

# the first step is to blast consensus of repeat elements against the ncbi non-redundant database of proteins
# we thus collect all the species TE consensuses, and put them in a single fasta file (we prefer doing a single blast and we need to add species names to the sequences)

consensusFiles = list.files(pattern = "families.fa$", full.names = T, recursive = T)	
species = extractSpeciesNames(consensusFiles)		#species names must be present in file paths
import = function(file, sp) {
	seqs = readDNAStringSet(file)
	setNames(seqs, stri_c(sp, "-", names(seqs)))		#adds the species name to each family name
}

seqs = mcMap(import, consensusFiles, species, mc.cores = 10, mc.preschedule = F)
seqs = do.call(c, seqs)		#combines the sequences into a single DNAStringSet (unlist() wouldn't work)

setwd(workFolder)
writeXStringSet(seqs, "allTEconsensuses.fas")

## obtains the nr database (not recommended if you already have it somewhere)
download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz", "nr.gz")
## makes a diamond database
system("diamond makedb --in nr.gz -d nr -p 10")

out = "consensusesOn_nr.out"
fields = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle")
system(paste("diamond blastx -p 1 --more-sensitive -q allTEconsensuses.fas -d nr.dmnd --quiet -k 20 -f 6", paste(fields, collapse = " "), "-o", out))

blast_nr = fread(out, header = F, sep ="\t", col.names = fields)		#imports results
blast_nr[,c("sseqid", "stitle") := .(stitle, NULL)]					#replace subject id by subject title (description)
blast_nr = blast_nr[evalue <= 0.001]								#removes hits of low evalue
setorder(blast_nr, qseqid, qstart,-qend)								
blast_nr[,set := regionSets(data.frame(qseqid, qstart, qend))]		#removes hits nested within others (in the query), as these hits are of lower quality (shorter)
blast_nr = blast_nr[!duplicated(set)]

combined = combineHits(blast_nr, maxDist = 300, maxOverlap = 10000, maxDiff = 300, blastX = T)		#combining relatively adjacent hits (with a lot of margin). This function also renames columns
setorder(combined, query, -score)					#puts best hits on top, per query
combined[,hit := occurences(query)]					#attributes a number to every hit related to the same query (starting at 1)

#we do the same for hits of family consensuses against repeat_pep proteins
#we use RepeatPeps.lib = the repeat modeler classification fasta file (part of the RepeatModeler package)
#this file should be put in current working directory
system("diamond makedb --in RepeatPeps.lib -d repeatPeps")
out = "consensusesOnRepeatPeps.out"

system(paste("diamond blastx -p 1 --more-sensitive -q allTEconsensuses.fas -d repeatPeps.dmnd --quiet -k 20 -f 6", paste(fields, collapse = " "), "-o", out))

blast_repBase = fread(out, header = F, sep ="\t", col.names = fields)

blast_repBase = blast_repBase[evalue <= 0.001]
setorder(blast_repBase, qseqid, qstart,-qend)
blast_repBase[,set := regionSets(data.frame(qseqid, qstart, qend))]
blast_repBase = blast_repBase[!duplicated(set)]
blast_repBase[,sseqid := "repProtein"]			#replaces subject name by this generic term
combinedRep = combineHits(blast_repBase, maxDist = 100000, maxOverlap = 100000, maxDiff = 100000, blastX = T)		#here we combined all hsps of each consensus, regardless of the protein (not very subtle, but it simplifies a lot)

#to find consensuses that may not be TEs. We compare their hits on repbase proteins to hits on nr. But since there may be several hit on NR per consensus, we do it on a per-hit (not per-consensus) basis. The first hit per consensus (query) is processed first, then the second, etc.
dubious = function(occ) {
	merged = merge(combined[hit == occ], combinedRep, by = "query", suffixes = c("",".r"), all = T)		#confronting the hit on repbase to the hit on nr for each consensus
	merged[,leading := pmin(qStart.r, qEnd.r) - pmin(qStart, qEnd)]							#this is the length of the part that aligns on the nr protein before it aligns on the rep protein
	merged[,trailing := pmax(qStart, qEnd) - pmax(qStart.r, qEnd.r)]							#the length the part that aligns on the nr protein after it aligns on the rep protein
	m = merged[is.na(subject.r) | leading > 90 | trailing > 90]						#we retain pairs of hits where there is no alignment on repbase, or where there are parts aligning on nr protein not covered by alignment on rep proteins (of at least 90 bp)
	m
}

dubiousConsensusHits = lapply(unique(combined$hit), dubious)
dubiousConsensusHits = rbindlist(dubiousConsensusHits)
dubiousConsensusHits = dubiousConsensusHits[!grepl("transpo|retro|reverse|integrase|gag|pol-|pol |rna-dependent|polyprot|mobile|jockey|jerky|setmar|copia|recombinase|crypton|mariner|tcmar|tc1|gypsy|helitron|harbi|piggy|maveri|polinton|academ|ltr|cmc|envelop", subject, ignore.case = T) & abs(qStart-qEnd +1) > 100,]		#retains hits to nr proteins that are not annotated as repeat proteins and with alignment length >100
dubiousConsensusHits[,acc := splitToColumns(subject, " ",1)]		#puts accession numbers of nr proteins in this new column

#we now blast proteins similar to TE consensuses against repeat proteis
writeLines(unique(blast_nr$sseqid), "hitted_nrProteins.bed")	#we put their names in a file for seqtk (we use all proteins, not just those in the dubiousConsensusHits table. We can afford it)
seqtk(nr.gz, bed, "hitted_nrProteins.fas")

out = "nrProtsOnRepeatPeps.out"
system(paste("diamond blastp -p 1 -q hitted_nrProteins.fas --more-sensitive -d repeatPeps.dmnd --quiet -k 20 -f 6", paste(fields, collapse = " "), "-o", out))

blastp = fread(out, header = F, sep ="\t", col.names = fields)
blastp = blastp[evalue < 0.001]
combinedP = combineHits(blastp, maxDist = 1000, maxOverlap = 100000, maxDiff = 1000)
repProtAccessions = combinedP[abs(qEnd-qStart)+1 >= 100 & identity >= 35, unique(splitToColumns(query, " ",1))]		#gets accession number of nr proteins that have an acceptable hit with a rebpase proteins. We will consider them as legit TE proteins
dubiousFamilies = dubiousConsensusHits[!acc %in% repProtAccessions, unique(query)]	#so we can extract consensuses that hit to other proteins than these 
write(dubiousFamilies, "familiesToIgnore.txt")			#writes these family names to disk. We still blast these TEs to find HTT. Removal of these TEs is done afterwards