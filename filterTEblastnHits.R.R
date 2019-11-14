# this scripts processes "raw" blast output to select hits with sufficient pID
library(data.table)
library(parallel)
library(stringi)
args = commandArgs(trailingOnly=TRUE)   
nCPUs = as.integer(args[1])

dir.create("TEs/blastn/filtered/")
pairs = fread("pairsToFilter.txt")
pairs[, out := gsub("/out/", "/done/", stri_c(out, ".gz"), fixed = T)]				#as output files have been moved and compressed in pairwiseSpeciesBlastn.R
forw = pairs[sp1 < sp2] ; rev = pairs[sp1 > sp2]			#we will process forward and reverse blast searches (as each species was used as query and db) simultaneously, so as to select the best hit among the two reciprocal ones
setnames(rev, 1:3, c("sp2","sp1","rev"))
pairs = merge(forw, rev)
pairs[,filtered := gsub(".gz", "", stri_c("TEs/blastn/filtered/", basename(out)))]		#paths of filtered output files to be generated

toRemove = readLines("TEs/findDubious/familiesToIgnore.txt")			#names of dubious TE families that may not represent true TEs generated in step 3-findDubiousTEs.R



filterHits = function(out, rev, filtered, q) {
        blast = data.table()				#creates an empty data table to avoid a bug in fread() if a blast output file is empty
        nr = 0								#to report the number of hits					
        if(file.size(out) > 0) {			#imports "foward" blast search between two species
                blast = fread(paste("gunzip -c", out), sep = "\t", header = F, drop = c(5, 6, 11))
                nr = nr + nrow(blast)
                blast = blast[V3 > q,]      #removes hits of insufficient pID
        }
        reverse = data.table()				#we do the same for the "reverse" search
        if(file.size(rev) > 0) {
                reverse = fread(paste("gunzip -c", rev), sep = "\t", header = F, drop = c(5, 6, 11))
                nr = nr+nrow(reverse)
                reverse = reverse[V3 > q,]
                reverse[,c("V1","V2","V7","V8","V9","V10") := .(V2, V1, V9, V10, V7, V8)]        #reversing query and subject so we can concatenate forward and reverse hits, and remove reciprocal hits (amongst others)
        }
        if(nrow(blast) == 0 & nrow(reverse)==0) { 
        	file.create(filtered)			#creates an empty output files if there is no retained hit (not sure if this was needed)
        	return(NULL) 
        }
        blast = rbind(blast, reverse)
        rm(reverse)                                             
        setorder(blast, -V12)                                                   #puts best hits on top (column 12 is the score)
        blast = blast[!duplicated(data.table(V1, V2))]                          #retaining the best hsp per pair of copies, which also removes reciprocal hits
        
        #now we remove dubious TE families and create columns that will be used through our pipeline
        mat1 = stri_split(blast$V1, fixed = ":", simplify = T)					#splits copy names into matrices of several columns (fields are separated by colons)
        mat2 = stri_split(blast$V2, fixed = ":", simplify = T)
        f1 = stri_c(mat1[,1], "-", mat1[,6]) ; f2 = stri_c(mat2[,1], "-", mat2[,6])		#creates TE family names from species names (first field) and family names
        f = !f1 %chin% toRemove & !f2 %chin% toRemove							#so we can filter-out hits involving dubious TEs families
        if(sum(f)==0) { file.create(filtered); return(NULL) }					#again creates an empty output file if there is no retained hit
        mat1 = rbind(mat1[f,]) ; mat2 = rbind(mat2[f,])							#rbind() makes sure that these will be matrices and not vectors, even if there is only one hit (to avoid a dimension problem)
        subf1 =  stri_split(mat1[,6], fixed = "#", simplify = T)				#separates TE family name form superfamily name
        subf2 =  stri_split(mat2[,6], fixed = "#", simplify = T)		
        query = stri_c(mat1[,2], mat1[,3], mat1[,4], mat1[,5], sep = ":")		#creates new copy names (not containing species names, as these will be in a separate column) 
        subject = stri_c(mat2[,2], mat2[,3], mat2[,4], mat2[,5], sep = ":")	
        nr2 = nrow(blast)
        blast = data.table(query, subject, mat1[,1], mat2[,1], subf1, subf2, blast[f,-c("V1", "V2")])     
        fwrite(blast, filtered, col.names = F, row.names = F, quote = F, sep ="\t")		#write file of filtered hits (not returned by the function, to save memory and because Map cannot apparently return tables)
        c(basename(filtered), nr, nr2, nrow(blast))			#instead we return the initial number of hits and number of retained hits
}


m = pairs[!file.exists(filtered), mcMap(filterHits, out, rev, filtered, q, mc.cores = nCPUs, mc.preschedule = F)]
m= data.table(do.call(rbind, m))
fwrite(m, "filteredStats.txt", col.names = F, row.names = F, quote = F, sep ="\t")
system("cat TEs/blastn/filtered/*.out > all.quantile005score200.blastn.out")		#concatenates output files

