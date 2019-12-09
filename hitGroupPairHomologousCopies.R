##%######################################################%##
#                                                          #
####   This script finds homologies between TE copies   ####
####      belong to different hit groups, which we      ####
####        will use to infer the minimal number        ####
####       of HTT events to explain these groups        ####
#                                                          #
##%######################################################%##


source("HTvFunctions.R")

args = commandArgs(trailingOnly = TRUE)
nCPUs = as.integer(args[1])
dir.create("TEs/clustering/networks")
hitList = readRDS("TEs/clustering/forHitGroupPairs.RDS")              #a list of hits (data tables with query id, subject id, pID, hit group id) in each super family
blastFiles = list.files("TEs/clustering/blastn/done",
                        pattern = ".out",
                        full.names = T)            #blast files of copies against themselves
superF = testSplit(basename(blastFiles), "_", 1)						#we split these by superfamily (large superfamilies have several blast files)
blastFiles = split(blastFiles, superF)

hitGroupPairStats = function(superF) {
        hits = hitList[[superF]]                        #the HTT hits of the super family being processed
        copies = hits[, unique(c(q, s))]         		#all the copies involved in these hits
        import = function(file) {
                #imports a blast file of copies against themseves
                dt = fread(
                        file,
                        sep = "\t",
                        header = F,
                        select = c(1:3, 9),
                        nThread = 2
                )
                dt = dt[V1 %in% copies & V2 %in% copies]
                dt[, V3 := as.integer(V3 * 1000)]          #converts pID to integer for speed and memory reasons
                dt
        }
        dt = rbindlist(mclapply(
                blastFiles[[superF]],
                import,
                mc.cores = 5,
                mc.preschedule = F
        ))
        dt = rbind(dt[, .(V1, V2, V3)], dt[, .(V1 = V2, V2 = V1, V3)], data.table(
                V1 = copies,
                V2 = copies,
                V3 = 100000L
        ))            #duplicate the hits by creating reciprocal ones (used below) and adds self hits with 100% pIDs
        subjectsForQuery = reList(split(dt$V1, dt$V2))  #so for each query, which is the index of this list, we have the subjects (as integers)
        nSubjects = sapply(subjectsForQuery, length)    #we record the number of subjects per query, used later
        idForQuery = reList(split(dt$V3, dt$V2))                        #we also record the pID scores
        rm(dt)  #reclaims RAM, as there are lots of hits
        copyPerHitGroup = hits[, .(copy = unique(c(q, s))), by = hitGroup]               		#unique copy ids for each hit group
        hitGroupsForCopy = reList(split(copyPerHitGroup$hitGroup, copyPerHitGroup$copy))        #for each copy, which represents an index in this list, we have the hit group(s) where it is found
        copyPerHitGroup = split(copyPerHitGroup$copy, copyPerHitGroup$hitGroup)
        nHitGroups = sapply(hitGroupsForCopy, length)                                         #and the number of hit groups where it is found
        
        hitGroupLinks = function(copies, hitGroup1) {
                #for copies in each hit group (named hitGroup1), this retreives the hit groups with homologous copies
                subject = unlist(subjectsForQuery[copies])              #homologous copies to those in the hit group (subjects of the blast)
                pID = unlist(idForQuery[copies])                         #with associated pIDs
                pairs = data.table(copy = rep(copies, nSubjects[copies]), subject , pID)               #we place them in a data table (which requires replicating copies that have several subjects)
                hitGroup2 = unlist(hitGroupsForCopy[pairs$subject])           #and we get the hit groups to which these subjects belongs
                len = nHitGroups[pairs$subject]
                pairs = pairs[, data.table(copy = rep(copy, len),
                                           pID = rep(pID, len),
                                           hitGroup2)]          #places all this in a data.table
                pairs = pairs[hitGroup2 != hitGroup1]
                setorder(pairs,-pID)
                pairs = pairs[!duplicated(data.table(copy, hitGroup2))] 	#for each hitGroup2, we retain only the most similar subject per copy of hitGroup1, as we don't need more (and saves some memory)
                pairs[, hitGroup1 := hitGroup1]			#adds the hitGroup1 id to the table. Doing it earlier would have increased the memory footprint.
                cat(".")
                pairs
        }
        hitGroupPairs = Map(hitGroupLinks, copyPerHitGroup, as.integer(names(copyPerHitGroup)))         #applies the above for all hit groups of the super family
        writeT(
                rbindlist(hitGroupPairs),
                stri_c(
                        "TEs/clustering/networks/",
                        superF,
                        ".hitGroupPairCopies.txt"
                )
        )
        rbindlist(hitGroupPairs)
}

res = mclapply(names(hitList),
               hitGroupPairStats,
               mc.cores = nCPUs,
               mc.preschedule = F)
writeT(rbindlist(res),
       stri_c("TEs/clustering/networks/all.hitGroupPairCopies.txt"))
