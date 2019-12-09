##%######################################################%##
#                                                          #
####           This stage computes Ks between           ####
####        species for core genes extracted by         ####
####             BUSCO, in order to compare             ####
####           Ks between TEs and core genes            ####
#                                                          #
##%######################################################%##

######### The BUSCO pipeline was not automated in scripts
# example of BUSCO run (launched from bash):
# BUSCO_v3/scripts/run_BUSCO.py -i genomes/Numida_meleagris/Numida_meleagris.fa -o Numida_meleagris.fa -l database/busco/aves_odb9/ -m genome -c 1 -t 10.busco/output/tmp/


######### we first translate BUSCO CDS into proteins in order to run diamond blastp between species, so as to find orthologs between which we compute Ks, and to align proteins afterwards (we use protein alignments as guides since Ka Ks are based on codons)
source("HTvFunctions.R")

folders = list.dirs("BUSCO_results")		#folder assumed to contain BUSCO results
folders = folders[grep("single_copy_busco_sequences", folders, fixed = T)]		#folder path must contain the name of a species corresponding exactly to a tip name of timetree.nwk

dir.create("CDS") 	#we will put BUSCO CDS fasta files here (one file per species), this is done by the function below
m = mclapply(folders,
             getBUSCOs,
             mc.cores = 10,
             mc.preschedule = F)

CDSfiles = list.files("CDS", pattern = ".CDS.fas", full.names = T)		#listing files of BUSCO CDS
aaFiles = gsub(".CDS.fas", ".aa.fas", CDSfiles)
f = !file.exists(aaFiles)		#to avoid redoing translations that were already done if the job needs to be relaunched

translationSummary = mcMap(translateCDS,
                           CDSfiles[f],
                           aaFiles[f],
                           mc.cores = 10,
                           mc.preschedule = F)	#translation done in parallel
translationSummary = rbindlist(translationSummary)
translationSummary[, sp := extractSpeciesNames(cds)]

writeT(translationSummary, "CDS/CDS_translation_summary.txt")

################################### making diamond databases, one per species
dir.create("Ks/blastp", recursive = T)
folders = stri_c("Ks/blastp/", extractSpeciesNames(aaFiles))		#we create one subfolder per species
dbs = stri_c(folders, "/", extractSpeciesNames(aaFiles), ".dmnd")	#paths to diamon database to be generated
lapply(folders, dir.create)
f = !file.exists(dbs)

makedb = function(folder, db, fas) {
  #makes a diamond db
  system(paste("diamond makedb --in", fas, "-d", db))
}
m = mcMap(makedb,
          folders[f],
          dbs[f],
          AAfiles[f],
          mc.cores = 10,
          mc.preschedule = F)  #databases generated in parallel

################################### defining pairwise diamond blastp searches between species
#we use the timetree to reduce the number of searches. HTT is not inferred between closely related species, and there is no need to measure Ks in all possible species pairs between two large clades (only a subset of species is sufficient)
tree = read.tree("timetree.nwk")
distMat = cophenetic(tree)
mrca = mrca(tree)							#the matrix listing the MRCA of every tip (as a row and column index)
#we turn these matrix into a data table for every species pair
pairs = data.table(
  divTime = as.vector(distMat),
  sp1 = rep(colnames(distMat), nrow(distMat)),
  sp2 = rep(rownames(distMat), each = ncol(distMat))
)
pairs[, mrca := mrca[cbind(sp1, sp2)]]
aaFiles = list.files("CDS", pattern = "aa.fas", full.names = T)
aaSpecies = extractSpeciesNames(basename(aaFiles))							#species for which we have translated BUSCO (at least one did not have annotated genes)
pairs = pairs[sp1 %chin% aaSpecies &
                sp2 %chin% aaSpecies]					#only retains pairs of species that have BUSCO genes

pairs[, db := stri_c("Ks/blastp/", sp1, "/", sp1, ".dmnd")]					#generating file names for every blastp. First, the database
pairs[, out := stri_c("Ks/blastp/", sp1, "/", sp2, ".on.", sp1, ".out")]		#then the output
pairs[, aa := stri_c("CDS/", sp2, ".aa.fas")]								#then the inputs
pairs[, daa := stri_c("Ks/blastp/", sp1, "/", sp2, ".on.", sp1, ".daa")]		#and the daa files

#for clades older than 500 My, which are very large, we do not compute Ks between all species as it is overkill. We will select one species per smaller "young" clades of <=30 My, based on the number of BUSCO detected in its genome

youngClades = cladesOfAge(tree, 30, withTips = T)		#all these young clades with their species (see function in HTvFunctions.R)
youngClades[, nCDS := translationSummary[match(tip, sp), nseq - stops]]	#number of translated BUSCO cds for these species
setorder(youngClades,-nCDS)						#we put species that have more BUSCO genes on top
toIgnore = youngClades[duplicated(node), tip]		#so we may ignore all the other species for each clade

# if the divergence of clades is >=500 My, we use only one species per young subclade, else we use all species, except when the divergence is <80 My (which we discard)
pairs = pairs[(divTime >= 500 &
                 !sp1 %chin% toIgnore &
                 !sp2 %chin% toIgnore)  | (divTime > 80 & divTime  < 500)]
pairs[!file.exists(out), batch := rep(1:10, length.out = .N)]		#we create 10 batches of blastp jobs (again, using an sarray of jobs would probably have been better)
writeT(pairs, "pairsToBlastp.txt")

system(
  'sbatch --mail-type=BEGIN,END,FAIL --cpus-per-task=30 --mem=50G --wrap="Rscript coreGenesblastp.R 1 30"'
)  # example of the 1st job with 30 CPUs

########################## computing Ks from results
system("Rscript coreGenesKaKs.R")

############ establishing Ks distribution between sister clades from Ks obtained
Ks = mclapply(list.files("Ks/KaKs", full.names = T), function(x)
  fread(x, header = F, sep = "\t"), mc.cores = 3)	#importing all results off KaKs computations
Ks = rbindlist(Ks)
setnames(Ks,
         c(
           "query",
           "subject",
           "qStart",
           "qEnd",
           "sStart",
           "sEnd",
           "Ks",
           "Ka",
           "alnLength"
         ))
Ks[, sp1 := extractSpeciesNames(query)]			#extracts species names from sequence names
Ks[, sp2 := extractSpeciesNames(subject)]

Ks[, divTime := distMat[cbind(sp1, sp2)]]						#adds the divergence time of each pair of species, for later use
#determining the MRCA (node/clade number) for all pairs of species, which we will use to establish the ditribution of Ks values for each clade
Ks[, clade := mrca[cbind(sp1, sp2)]]							#so we add the MRCA (clade) of each species pair as a new column. These two species de facto belong to the two sister clades that diverged from this MRCA. So the "clade" column also designates the two sister clades composing it

#For a pair of sister clades comprising many species, the same BUSCO (of a given species) yielded many Ks values. To reduce the pseudo-replication, we retain one Ks value per busco gene per species, using its longest alignment.
setorder(Ks,-alnLength)												#puts longest alignments on top
Ks = Ks[alnLength > 200 &
          !duplicated(data.table(clade, subject)) &
          !duplicated(data.table(clade, query))]		#ignores alignments shorter than 200 aa and retains no more than one alignment per BUSCO gene per clade

writeT(Ks, "Ks200AAnoRedundancy.txt")

# making figure S2 to illustrate the rate of synonymous molecular evolution
timeRanges = seq(0, sqrt(max(Ks$divTime, na.rm = T) + 1), length.out = 20) ^
  2					#generates divergence time classes within which we compute average Ks (all species pairs considered). We use sqrt() to have shorter intervals for low divergence times (using quantile() to generate breaks yields too many points at low divergence times, so the plot looks biased)
Ks[, range := .bincode(divTime, timeRanges)]			#assigns divergence times to the classes
KsDistrib = Ks[Ks < 9 &
                 !is.na(divTime), .(time = mean(divTime), Ks = weighted.mean(Ks, w = alnLength)), by = range]		#ignoring Ks values ≥ 9 (oversaturated)

pdf("figure S1.pdf")
p = KsDistrib[, plot(
  time,
  Ks,
  ylim = c(0, 3),
  xlab = "Divergence time (My)",
  ylab = "Ks",
  pch = 16,
  col = "darkgrey"
)]
fit = nls(Ks ~ max * a * time / (max + a * time),
          data = KsDistrib,
          start = list(a = 0.01, max = 3))	#fits a curve assuming initial linear correlation between Ks and divergence time (coeff a) and saturation at Ks = max (as used for insects)
x = KsDistrib[, seq(0, max(time), length.out = 150)]											#to add a smooth curve corresponding to the model, we generate 150 divergence time values (x axis)
lines(x, predict(fit, list(time = x)))
dev.off()
