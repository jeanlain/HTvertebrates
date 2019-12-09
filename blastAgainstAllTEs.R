##%######################################################%##
#                                                          #
####            this script blasts TE copies            ####
####           in retained hit groups againts           ####
####                all TEs of the host                 ####
####        species, to check the reliability of        ####
####            transfers (by the number of             ####
####          copies they may involve), so as           ####
####            to test whether apparent HTT            ####
####       may in fact result from contamination        ####
#                                                          #
##%######################################################%##


source("HTvFunctions.R")
args = commandArgs(trailingOnly = TRUE)
nCPUs = as.integer(args[1])

fasFiles = list.files("TEs/clustering/testConta",
                      pattern = ".fas",
                      full.names = T)		#query file names
sp = extractSpeciesNames(basename(fasFiles))

dbs = list.files("TEs/blastn/db", pattern = ".nin", full.names = T)			#blast databases of TEs for all species, generated in step 4-blastTEs.R
dbs = gsub(".nin", "", dbs)
sp2 = TEs / clustering / testConta / blastn / (basename(dbs))
dbs = dbs[match(sp, sp2)]													#we order db files to match query files
dir.create("TEs/clustering/testConta/blastn/done", recursive = T)
outFiles = stri_c("TEs/clustering/testConta/blastn/",
                  sp,
                  ".selectedCopies.blastn.out")		#blastn output files
doneFiles = stri_c("TEs/clustering/testConta/blastn/done/",
                   sp,
                   ".selectedCopies.blastn.out")	#same, but after a blast search is done

f = !file.exists(doneFiles)				#so we can avoid redoing completed blast searches

blastAgainstTEs = function(fasFile, db, out, done) {
  system(
    paste(
      "blastn -query",
      fasFile,
      "-db",
      db,
      "-outfmt 6 -max_target_seqs 20 -max_hsps 1 -out",
      out
    )
  )
  file.rename(out, done)				#moves result to other folder once finished
  NULL
}

m = mcMap(
  blastAgainstTEs,
  fasFiles[f],
  dbs[f],
  outFiles[f],
  doneFiles[f],
  mc.cores = nCPUs,
  mc.preschedule = F
)
system(
  "cat TEs/clustering/testConta/blastn/done/*.out > TEs/clustering/testConta/blastn/all.out"
)

print("finished")
