##%######################################################%##
#                                                          #
####      This script computes pairwise Ka/Ks and       ####
####  overal molecular distance between homologous TEs  ####
#                                                          #
##%######################################################%##

# this script is run at stages 7 and 13 of the pipeline

source("HTvFunctions.R")
require(seqinr)

args <- commandArgs(trailingOnly = TRUE)

# the arguments are
# - the path to the tabular file of HSPs between TE sequences, with columns query, subject, qStart, qEnd, sStart, sEnd
TEhitFile <- args[1]

# - the path to a tabular file of HPSs between these sequences and proteins.
# These  HSPs mustn ot overlap on a given TE (see step 7)
# This file as fields: TE copy name, start and end coordinates of HSP on the TE (start < end),
# the start coordinate of the HSP on the protein, and a logical telling whether the HSP is
# in reverse orientation in respect to the TE
# All values in query and subject columns in the first file but be found in the copy column of this file
blastxFile <- args[2]

# - the path to the fasta file containing the TE sequences. 
fastaFile <- args[3]

# -folder where results file will go. 
outputFolder <- args[4]

# - the number of CPUs to use
nCPUs <- as.integer(args[5]) 

dir.create(outputFolder, showWarnings = F)



# STEP ONE, preparation of the TE-TE copy alignments ---------------------------------------------------------------

# We import HSP files
TEhits <- fread(TEhitFile)
blastx <- fread(blastxFile)

# we only retain blastx HSPs that involve sequences in the TE-TE hits
blastx <- blastx[copy %chin% TEhits[,c(query, subject)]]

# to keep track of hits, we give them an integer identifier in a column
TEhits[, hit := 1:.N]

# We select HSPs for which there is enough protein regions per copy
# we first compute the total length of HSPs per copy
covPerCopy <- blastx[, sum(end - start + 1L), by = copy]

# the query and subject TEs must have an alignment on proteins longer than 30 bp
# "hits" below are simply row indices of TEhits table
hits <- TEhits[, which(
  query %chin% covPerCopy[V1 > 30L, copy] & subject %chin% covPerCopy[V1 > 30L, copy]
  )]

# we split the work into several batches (jobs) for parallelisation
# The size of a batch of HSPs (4000) must not be too big due to memory constrains.
nJobs <- max(nCPUs, ceiling(length(hits) / 4000))

# we split the HSPs into the batches
hits <- splitEqual(sample(hits), n = nJobs)
# Note that we randomise HPSs (rows) to reduce differences in job durations 

# we import the TE sequences
seqs <- readDNAStringSet(fastaFile) 

# we modify names to match copy names of the TEhits table (the fasta
# file has longer sequences names that also comprise the host species)
names(seqs) <- copyName(names(seqs))

# we extract the parts of sequences involved in TE hits
qSeqs <- TEhits[, subSeq(seqs[query], qStart, qEnd)]
sSeqs <- TEhits[, subSeq(seqs[subject], sStart, sEnd)]

# we splits these subsequences into batches corresponding to the hit batches to be processed in parallel
qSeqs <- lapply(hits, function(batch) qSeqs[batch])
sSeqs <- lapply(hits, function(batch) sSeqs[batch])

# and we do the same for the hits themselve (we only retain coordinates of hits)
TEhits <- lapply(hits, function(batch) {
      TEhits[batch, .(query, subject, qStart, qEnd, sStart, sEnd)]
  })

rm(seqs, covPerCopy) # we reclaim some RAM




# STEP TWO, TE-TE pairwise sequence alignment and Ka/Ks computation ----------------------------------------------------------------

# below is the core function that pairwise aligns TEs and computes Ka/Ks on a batch of hits (jobs).
KaKsForJob <- function(job) { # the only argument is the job number

    # we align pairs of copies (this uses the Biostrings package)
    aln <- alignWithEndGaps(
        seq1 = qSeqs[[job]],
        seq2 = sSeqs[[job]]
    )
    
    # we retreive the table of coordinates of hits corresponding to this job
    TEhitBatch <- TEhits[[job]]

    # we determine position of bases within codons --------------------------------

    # We split alignments into a table of individual nucleotides and positions
    # see function definition in HTvFunctions.R
    nuc <- splitAlignment(
        aln = aln,
        coords = TEhitBatch[, .(query, subject, qStart, qEnd, sStart, sEnd)]
    )

 
    # in this script, we also measure some overall molecular distance between copies,
    # which we will use later.
    # For this, we need to make DNAbins (used by ape)
    sequencePairs <- Map(
      f = list, 
      pattern = split(nuc$base1, f = nuc$aln), 
      subject = split(nuc$base2, f = nuc$aln)
      )
    
    sequencePairs <- lapply(sequencePairs, as.DNAbin)
    
    rawDistance <- sapply(sequencePairs, dist.dna, model = "raw")
    
    # Kimura 80 is the model on which the Ka Ks estimate is based on
    K80distance <- sapply(sequencePairs, dist.dna, model = "K80")
    
    # to determine the position of bases within codons, we determine 
    # the blastx HSP that covers each aligned position, for each copy.
    # for this we add a new integer column simply denoting the row index of the HSP in the blastx table
    # alignment positions outside these HSPs will get NA values
    # we begin with the "query" TE
    nuc[, hsp1 := assignToRegion(
        bed = blastx[, .(copy, start, end)],
        pos = .(seq1, pos1)
    )]

    # then we process the subject TE copy
    nuc[, hsp2 := assignToRegion(
        bed = blastx[, .(copy, start, end)],
        pos = .(seq2, pos2)
    )]

    
    # we can already discard positions outside protein HSPs (no evidence that they are in ORFs)
    nuc <- nuc[!is.na(hsp1) & !is.na(hsp2)]

    # we now determine the position in codon ="frame" (1 to 3 or -1 to -3) of each position
    # for this, we need to know the the orientation of the protein in respect to the copy
    # we thus add two logical column that tell whether the blastx HSPs are in reversed (for query and subject)
    nuc[, c("rev1", "rev2") := .(
        blastx[hsp1, rev],
        blastx[hsp2, rev]
    )]

    # we determine the frame of the "query" base (base1)
    nuc[, frame1 := ifelse(
      test = !rev1, 
      yes = (pos1 - blastx[hsp1, protStart]) %% 3L + 1L,
      no = -((blastx[hsp1, protStart] - pos1) %% 3L + 1L)
    )]

    # then for the subject base (base2)
    nuc[, frame2 := ifelse(
      test = !rev2,
      yes = (pos2 - blastx[hsp2, protStart]) %% 3L + 1L,
      no = -((blastx[hsp2, protStart] - pos2) %% 3L + 1L)
    )]

    # we reverse the frame of the query bases if it is itself reversed in the TE-TE hit
    nuc[TEhitBatch[aln, qStart > qEnd], frame1 := -frame1]

    # same for the subject
    nuc[TEhitBatch[aln, sStart > sEnd], frame2 := -frame2]

    # we remove aligned position where the frames (position in codons) don't match
    nuc <- nuc[frame1 == frame2]
    
    # as well as columns that are no longer necessary (to free some ram)
    nuc[, c("hsp1", "hsp2", "rev1", "rev2", "frame2") := NULL] 



    # we discard incomplete codons from alignments  -----------------------
    # Codons will be identified by integers (first one, second one, etc.)

    # we first determine that a new codon starts when...
    shift <- nuc[, c(
        T, # it is the very first base of the table or
        # the absolute difference in "frame" between successive bases is not 1, or
        abs(diff(frame1)) != 1L | 
        # the absolute difference in coordinates in the query copy is not 1, or
        abs(diff(pos1)) != 1L | 
        # same for the subject copy
        abs(diff(pos2)) != 1L | 
        # or if we shift to a different alignment
        diff(aln) != 0
    )] 

    # each aligned position gets ab identifier that differs between codons
    codon <- cumsum(shift)

    # we retain only complete codons (found 3 times)
    nuc <- nuc[codon %in% which(tabulate(codon) == 3L)]



    # we prepare, then do, the Ka Ks computations -----------------------------
    # we discard alignments that are too short
    nuc <- nuc[aln %in% which(tabulate(aln) > 30L)]

    
    # we count the number of substitutions per alignment
    nSubstitutions <- nuc[,.(nMut = sum(base1 != base2)), by = aln]
    
    # we concatenate the individual bases back into sequences, in a new table
    flattenedSequences <- nuc[, .(
        query = stri_flatten(base1),
        subject = stri_flatten(base2),
        frame = frame1[1] # we need this column below
    ), 
    by = aln
    ]

    # we reverse-complement the sequence that are not from the coding strand	(frame is negative)
    flattenedSequences[frame < 0L, c("query", "subject") :=
        .(revCom(query), revCom(subject))]

    # we produce seqinr alignment objects that are required for kaks()
    alns <- Map(
      f = function(seq1, seq2) {
        seqinrAlignment(c(seq1, seq2))
      }, 
      flattenedSequences$query,
      flattenedSequences$subject
    )

    # and compute Ka and Ks with seqinr
    KaKs <- lapply(alns, kaks)

    # we stack the results into a data.table
    KaKs <- rbindlist(lapply(
        X = KaKs,
        FUN = as.data.table
    ))

    
    # we add columns for HSP identifiers and alignment lengths
    alnID <- flattenedSequences$aln
    
    KaKs = data.table(
        hit = hits[[job]][alnID],
        KaKs,
        length = nchar(flattenedSequences$query),
        nMut = nSubstitutions[match(alnID, aln), nMut],
        K80distance = K80distance[alnID],
        rawDistance = rawDistance[alnID]
    )

    # we write results for safety (we write the combined results at the end)
    writeT(KaKs, stri_c(outputFolder, "/KaKs.", job, ".txt"))

    KaKs
}

# we apply the above function for the different jobs in parallel
res <- mclapply(
    X = 1:length(TEhits),
    FUN = KaKsForJob,
    mc.cores = nCPUs,
    mc.preschedule = F
)


writeT(data = rbindlist(res), stri_c(outputFolder, "/allKaKs.txt"))

# if we have reached this step, we may remove intermediate files
unlink(stri_c(outputFolder, "/KaKs.*.txt"))
