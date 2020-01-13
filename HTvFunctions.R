## %######################################################%##
#                                                          #
####               functions used by the                ####
####        different scripts used the analysis         ####
#                                                          #
## %######################################################%##

# some of the functions are used only once in the analysis, but they were not written specifically for it


# We also use this script to check and install the required packages
packages <- c(
  "parallel", "stringi", "data.table", "matrixStats", 
  "Biostrings", "ape", "igraph", "seqinr", "RColorBrewer"
  )

missing <- setdiff(packages, rownames(installed.packages()))

repository = "https://cran.us.r-project.org"

if(length(setdiff(missing,"Biostrings")) >0) {
  install.packages(
  pkgs = setdiff(missing,"Biostrings"), 
  repos = repository
  )  
}

if ("Biostrings" %in% missing) {
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = repository)
  }
  
  BiocManager::install("Biostrings")
}

stillMissing <- setdiff(packages, rownames(installed.packages()))

if(length(stillMissing) > 0) {
  stop(paste(
    "Package(s)", 
    paste(stillMissing), 
    "could not be installed. Consider installing it/them manually.")
    )
}

# we load the packages we need for the functions below
l = lapply(setdiff(packages, c("igraph","seqinr","RColorBrewer")), require, character.only = TRUE)

options("stringsAsFactors" = F)
options("scipen" = 20)




rescale <- function(x, range) {
    # rescales a numeric vector so that its range (min and max) 
    # corresponds to the range argument (a 2-number vector)
    
    ampl <- max(range) - min(range)
    amplx <- max(x, na.rm = T) - min(x, na.rm = T)
    xb <- x * ampl / amplx
    xb - min(xb, na.rm = T) + min(range)
}

mids <- function(v) {
    # returns the midpoints between any two successive elements of v (a numeric vector)
    
    i <- 2:length(v)
    (v[i - 1L] + v[i]) / 2L
}

Split <- function(x, f, drop = FALSE, sep = ".", recursive = F,...) {
  # splits a vector/table recursively by each factor of "f" if f is a list
  # and if "recursive" is TRUE
  
  ls <- split(x, f, drop = drop, sep = sep, ...)
  
  if(recursive & is.list(f)) {
    for(i in 2:length(f)) {
      fields = splitToColumns(names(ls), sep)
      names(ls) <- fields[,ncol(fields)]
      ls = split(
        x = ls, 
        f= data.frame(fields[,1:(ncol(fields)-1L)]), 
        drop = drop,
        sep = sep,
        ...)
    }
  }
  
  ls
}

reList <- function(list, len = NA) {
    # rearranges a list whose names are integers so that the names of 
    # elements of the list are placed at equivalent indices of a new list.
    # this can be used to optimise the speed of access to the list elements
  
    indices <- as.integer(names(list))
    
    if (any(is.na(indices))) {
        return(list)
        stop("some names are not convertible to integers")
    }
    
    m <- max(indices)
    
    if (!is.na(len) & len >= m) {
        m <- as.integer(len)
    }
    
    temp <- vector(mode = "list", m)
    temp[indices] <- list
    temp
}

toInteger <- function(string,
                      sorted = F,
                      unique = T) {
    # turns a string vector into integers
    
    if (is.numeric(string)) {
        return(string)
    } else {
        if (unique) {
            u <- unique(string)
        } else {
            u <- string
        }
        if (sorted) {
            u <- sort(u)
        }
        return(chmatch(string, u))
    }
}



splitEqual <- function(x, n, ...) {
    # splits an object (x, a vector or a table) into n classes of equal sizes, returns a list
  
    l <- length(x)
    if (is.data.frame(x)) {
        l <- nrow(x)
    }
    if (n == 1L) {
      f <- rep(1L, l)
    } else {
      f <- cut(1:l, n)
    }
    split(x, f, ...)
}




occurrences <- function(table) {
    # returns the occurence of each element of a vector or table (e.g. 1 if it's
    # the first time it appears, 2 if the 2nd time, etc). If vect is a table, it
    # gets the occurence of each row 
    
    dt <- data.table(table)
    dt[, pos_ := 1:.N]
    counts <- dt[, .(occ = 1:.N, pos_), by = setdiff(colnames(dt), "pos_")]
    res <- integer(nrow(dt))
    res[counts$pos_] <- counts$occ
    res
}


regionSets <- function(regions) {
    # gives a number ("set" identifier) to all regions that are included in the
    # same genomic region. "regions" is a 3-column data frame with sequence id,
    # start and end positions, ordered by sequence name, start coordinate, and -end
    
    
    pb <- txtProgressBar(
        char = "*",
        width = 100,
        style = 3
    )
    nrows <- nrow(regions)
    rows <- 2:nrows
    newSeq <- c(FALSE, regions[rows, 1] != regions[rows - 1, 1], T)
    ends <- c(regions[, 3], 0)
    set <- rep(NA, nrows)
    end <- ends[1]
    n <- 1
    for (i in 2:(nrows + 1)) {
        if (ends[i] > end | newSeq[i]) {
            set[n:(i - 1)] <- n
            n <- i
            end <- ends[i]
            setTxtProgressBar(pb, i / nrows)
        }
    }
    set
}


assignToRegion <- function(bed, pos, olv = T) {
    # assigns positions (like SNPs) to genome regions defined in a data frame in
    # bed-style format (sequenceID, start, end) and pos is a data frame specifying
    # sequenceID and position. If a position is in the area of overlap between 2
    # regions, it will be assigned to the second one
    
    bed <- as.data.table(bed)
    pos <- as.data.table(pos)
    setnames(bed, 1:3, c("sequenceID", "start", "end"))
    setnames(pos, 1:2, c("sequenceID", "pos"))

    # we will identify regions by their rows in the original bed
    bed[, id := 1:.N]
    bed <- bed[sequenceID %in% pos$sequenceID, ]
    
    # we convert coordinates into global coordinates to call .bincode once
    # This requires estimating the length of sequences, which is the max coordinate per sequenceID 
    seqLengths <- rbind(bed[, .(sequenceID, pos = end)], pos)[, max(pos) + 1L,
        by =
            sequenceID
    ]
    
    # we conver the coordinates for the bed and pos tables combined. 
    starts <- globalCoords(
        bed$sequenceID,
        bed$start,
        seqLengths$V1,
        seqLengths$sequenceID
    )
    
    pos <- globalCoords(
        pos$sequenceID,
        pos$pos,
        seqLengths$V1,
        seqLengths$sequenceID
    )
    
    # creates a data table with start and ends in the same "boundary" column.
    # "+1" gives a bit or headroom (probably not needed).
    # These will be the intervals for binning. We use "id" to keep track of 
    # the rows (regions id) of the original bed
    
    dt <- data.table(
        boundary = c(starts, starts + bed[, end - start + 1L]),
        id=rep(bed$id,2)
    )
    
    setorder(dt, boundary)

    # "region" will be the region (id) covering an interval (only for the start of this interval).
    dt[, c("row", "region") := .(1:.N, as.integer(NA))]
    
    # The idea is that when a region ends, the interval should get the region that main contain the one that ended or be NA (if between 2 regions)
    rows <- dt$row[-1]

    # first rows regions that do not overlap with others (whose ids appear successively in the table)
    short <- dt[, which(id[rows - 1L] == id[rows])]

    # all the rows of these regions (first and last)
    shorts <- c(short, short + 1L)
    
    if (length(shorts) < nrow(dt)) {
        # if there are overlapping regions
        
        # for each of these overlapping regions, we create a vector that represent all the rows of dt that are covered by that region
        if (length(shorts) == 0) {
            rows <- dt[, .(row = row[1L]:(row[.N] - 1L), size = row[.N] - row[1L]), by = id]
        } else {
            rows <- dt[-shorts, .(row = row[1L]:(row[.N] - 1L), size = row[.N] - row[1L]), by = id]
        }

        # puts the longest regions on top
        setorder(rows, -size)

        # so we can add a region id for each row in dt.
        dt[rows$row, region := rows$id]
        # If there are duplicates in rows$row, ids of shorter regions will overwrite longer ones. This takes advantage from the fact that r operates in vector order
    }

    # we now copy the ids for short regions.
    dt[short, region := id]
    dt[.bincode(pos, boundary, F, T), region]
}


globalCoords <- function(sequences, pos, seqLengths, seqNames) {
    # converts sequence/contig positions into unique genome positions assuming contiguous sequences whose lengths and name are provided in order

    
    names(seqLengths) <- NULL
    seqStart <- cumsum(as.numeric(c(1, seqLengths[-length(seqLengths)])))
    
    if (!is.numeric(seqNames)) {
        return(pos + seqStart[chmatch(sequences, seqNames)] - 1L)
    } else {
        return(pos + seqStart[match(sequences, seqNames)] - 1L)
    }
}

clusterize <- function(vect, limit, group = "n") {
    # takes a numeric vector and assigns elements to groups based on the maximum distance between them. Returns the groups as integers
    
    dt <- data.table(g = group, v = vect)
    dt$n <- 1:nrow(dt)
    dt <- dt[order(g, v), ]
    vect <- dt$v
    group <- dt$g
    pos <- 2:length(vect)
    
    w <- c(
        1,
        which(vect[pos] - vect[pos - 1] > limit |
            group[pos] != group[pos - 1]) + 1,
        length(vect) + 1
    )
    
    pos <- 2:length(w)
    groupSize <- w[pos] - w[pos - 1]
    group <- rep(1:length(groupSize), groupSize)
    group[order(dt$n)]
}



intersection <- function(start,
                         end,
                         start2,
                         end2,
                         range = F,
                         negative = T,
                         olv = 1L) {
    # returns the length (or range, if range is TRUE) of intersections for pairs of
    # ranges defined by their start and end coordinates (numeric vectors, typically
    # integers). "negative" allows for negative intersection lengths (when ranges
    # do no intersect). If FALSE, intersection lengths cannot be < 0. If the end of
    # a range is the same of the start of the other range, the overlap is considered
    # to be olv (a number)
    
    minEnd <- pmin(pmax(start, end), pmax(start2, end2))
    maxStart <- pmax(pmin(start, end), pmin(start2, end2))
    if (range) {
        r <- cbind(maxStart, minEnd)

        # when ranges do not intersect, the range of the intersection is NA
        r[maxStart > minEnd, ] <- NA
        return(r)
    }
    inter <- minEnd - maxStart + olv
    if (!negative) {
        inter[inter < 0] <- 0L
    }
    inter
}


genomeCov <- function(bed,
                      minCov = 1L,
                      seqLength = NULL,
                      seqNames = NULL,
                      combine = F,
                      successive = T) {
    # computes sequencing depth of a genome by aligned reads, given alignment
    # starts and ends of reads in a data table in bed-like format. Outputs
    # successive interval for a given depth. Ignores intervals covered by less than
    # minCov reads. seqLength specify the sequence length in case one wants
    # to record intervals with 0 depth (required for last intervals of contigs that
    # may have 0 detph). "combine" combines successive intervals covered by at
    # least 1 read (i.e. regions that were sequenced)

    # remembers column names of the bed
    nams <- names(bed)[1:3]
    setnames(bed, 1:3, c("sequenceID", "start", "end"))
    
    if (minCov == 0L &
        !is.null(seqLength) &
        !is.null(seqNames)) {
       
        # if one wants to record intervals with zero depth, we apply a trick where we add one artificial
        # read totally covering each contig. This adds 1X coverage to all positions, but afterwards we substract this.
        bed <- rbind(bed, data.table(
            sequenceID = seqNames,
            start = 1L,
            end = seqLength
        ))
    }
    
    # "melts" the bed to obtain 2 rows per read: alignment start and alignment end + an index (n) of 1 for starts and -1 for ends. 
    #We also converts sequenceID names to integers to save memory. We then add 1 to end positions, as this will facilitate merging intervals
    dt <- bed[, data.table(
        seqID = rep(toInteger(sequenceID, unique = F), 2),
        pos = c(start, end + 1L),
        n = rep(c(1L, -1L), each = .N)
    )]

    # orders the table by sequenceID and position (regardless of start and end)
    setorder(dt, seqID, pos)

    # this is the workhorse function to compute coverage at each position since reads starts add 1 to the coverage and ends remove 1
    dt[, cov := cumsum(n)]

    # reads starting or ending at the same position creates duplicate positions.
    # We only retain the last row among each set of duplicates as this is the one with accurate depth
    dt <- dt[!duplicated(data.table(seqID, pos), fromLast = T)]
   
     # to combine regions with at least 1X depth, we convert all depth > 1 into 1
    if (combine) {
        dt[cov > 1L, cov := 1L]
    }
 
    # we combine successive rows in the same sequenceID ifthey have identical depth
    # (may happen if a read ends just before another read starts). 
    # This effectively combines region with ≥ 1X if combine == T
    if (successive) {
        dt <- dt[c(T, diff(seqID) != 0L |  diff(cov) != 0L)]
    }
    
    # we convert the table back to interval format
    nr <- 2:nrow(dt) - 1L
    res <- dt[, data.table(
        seqID = seqID[nr],
        start = pos[nr],
        end = pos[nr + 1L] - 1L,
        cov = cov[nr]
    )]
    
    if (minCov == 0L) {
        # if we want intervals with 0 depth
        
        if (!is.null(seqLength) &
            !is.null(seqNames)) {
            res[, cov := cov - 1L]
        } else {
            nr <- 2:nrow(res)
            res <- res[c(seqID[nr - 1L] == seqID[nr] | cov[nr] == 0L, T)]
        }
    }
    
    res <- res[cov >= minCov]
    res[, seqID := bed[seqID, sequenceID]]
    setnames(res, 1:3, nams)
    setnames(bed, 1:3, nams)
    
    res
}

combineRegions <- function(bed, distance = 0L) {
    # combines regions (a 3-column data table with sequence name, start and end
    # position) that are distant by a certain 'distance'. By default it will
    # aggregate contiguous or overlapping regions, but distance could be > 0
    # (regions separated by distance pb or less) or negative (requires a certain
    # amount of overlap)
    
    bed <- as.data.table(bed)
    setnames(bed, c("sequenceID", "start", "end"))
    bed[end < start, c("start", "end") := .(end, start)]
    bed[, end := end + distance]
    res <- genomeCov(bed, combine = T)
    res[, end := end - distance]
    
    res[, -4, with = F]
}


compl <- function(bases) {
    # complements DNA/RNA bases in a character string
    
    chartr(
        "acgtmrwsykvhdbACGTMRWSYKVHDB",
        "tgcakywsrmbdhvTGCAKYWSRMBDHV",
        bases
    )
}

revCom <- function(seqs) {
    # reverse complements DNA sequences. Ambiguities are not treated
    
    com <- compl(seqs)
    seqs[] <- stri_reverse(com)
    seqs
}

writeT <- function(data, path, ...) {
    # convenience function to write a data.table to disk with commonly used settings
    
    fwrite(
        data,
        path,
        quote = F,
        row.names = F,
        na = "NA",
        sep = "\t",
        ...
    )
}


splitToColumns <- function(vect,
                           split,
                           columns,
                           empty = "",
                           mode = "character") {
    # split a character vector and returns a matrix. Characters are split upon each
    # occurence of the split argument (a character of length one). empty defines
    # what the empty cell should be, and mode is the mode of the returned matrix.
    # This basically does what stri_split_fixed() does, but may sometime be more
    # convenient or faster
    
    vect <- as.character(vect)
    nc <- stri_length(split)
    out <- NULL
    
    if (!missing(columns)) {
        cols <- as.integer(columns)
        cols <- cols[cols >= 1]
        out <- matrix(empty, length(vect), length(cols))
    } else {
        cols <- 1:1000
    }
    
    col <- 1
    
    repeat {
        temp <- vect
        pos <- as.vector(regexpr(split, vect, fixed = T)) - 1
        f <- pos >= 0
        temp[f] <- stri_sub(vect[f], 0, pos[f])
        
        if (!missing(columns)) {
            if (col %in% cols) {
                out[, match(col, cols)] <- temp
            }
        } else {
            out <- cbind(out, temp)
        }
        
        col <- col + 1
        
        if (col > max(cols) | !any(f)) {
            break
        }
        
        vect[f] <- stri_sub(vect[f], pos[f] + nc + 1, stri_length(vect[f]))
        vect[!f] <- empty
    }
    
    if (ncol(out) == 1) {
        out <- as.vector(out)
    }
    
    storage.mode(out) <- mode
    
    out
}


patternLeadingGap <- function(x) {
        # determines the length of pattern leading gaps of pairwise alignments generated by biostrings
    
        stopifnot(is(x, "PairwiseAlignments"), type(x) == "global")
        start(subject(x)) - 1L
    }

subjectLeadingGap <- function(x) {
    # determines the length of subject leading gaps of pairwise alignments generated by biostrings
    
    stopifnot(is(x, "PairwiseAlignments"), type(x) == "global")
    start(pattern(x)) - 1L
}


alignWithEndGaps <- function(seq1, seq2, ...) {
    # pairwise aligns sequences (DNAStringSets) with Biostring while preserving the
    # end gaps, adding "-" when necessary (Biostrings does not preserve end gaps by
    # default). Returns a data.table with two columns (character vectors) : aligned
    # pattern and aligned subject
    
    aln <- pairwiseAlignment(seq1, seq2, ...)
    
    aln <- data.table(
        pattern = as.character(pattern(aln)),
        subject = as.character(subject(aln)),
        pgap = patternLeadingGap(aln),
        sgap = subjectLeadingGap(aln)
    )
    
    f <- aln$sgap > 0
    missing <- stri_sub(seq1[f], 1, aln$sgap[f])
    aln[f, subject := stri_join(gsub("[A-Z,a-z]", "-", missing), subject)]
    aln[f, pattern := stri_join(missing, pattern)]
    
    f <- aln$pgap > 0
    missing <- stri_sub(seq2[f], 1, aln$pgap[f])
    aln[f, pattern := stri_join(gsub("[A-Z,a-z]", "-", missing), pattern)]
    aln[f, subject := stri_join(missing, subject)]

    subjectTrailingGap <- nchar(seq1) - nchar(gsub("-", "", aln$pattern))
    patternTrailingGap <- nchar(seq2) - nchar(gsub("-", "", aln$subject))

    f <- subjectTrailingGap > 0
    
    missing <- stri_sub(seq1[f], nchar(seq1[f]) - subjectTrailingGap[f] +
        1, nchar(seq1[f]))
    
    aln[f, subject := stri_join(subject, gsub("[A-Z,a-z]", "-", missing))]
    aln[f, pattern := stri_join(pattern, missing)]
    
    f <- patternTrailingGap > 0
    
    missing <- stri_sub(seq2[f], nchar(seq2[f]) - patternTrailingGap[f] +
        1, nchar(seq2[f]))
    
    aln[f, pattern := stri_join(pattern, gsub("[A-Z,a-z]", "-", missing))]
    aln[f, subject := stri_join(subject, missing)]
    
    return(aln[, .(pattern, subject)])
}



splitAlignment <- function(aln, coords, fillGaps = F) {
    # splits a pairwise sequence alignment data table returned by
    # alignWithEndGaps() into nucleotides, and generates coordinates of each
    # nucleotide according to "coords", which is a data.table with sequence names,
    # starts and ends of aligned regions for each sequence. deleted nucleotides are
    # not given coordinates, except if fillGaps is TRUE
    
    setnames(coords, c("seq1", "seq2", "start", "end", "start2", "end2"))
    nc <- nchar(aln$pattern)
    
    nuc <- data.table(
        base1 = unlist(strsplit(stri_flatten(aln$pattern), "")),
        base2 = unlist(strsplit(stri_flatten(aln$subject), "")),
        seq1 = rep(coords$seq1, nc),
        seq2 = rep(coords$seq2, nc)
    )
    
    nuc$pos1 <- 0L
    nuc$pos2 <- 0L
    nuc[base1 != "-", pos1 := unlist(Map(":", coords$start, coords$end))]
    nuc[base2 != "-", pos2 := unlist(Map(":", coords$start2, coords$end2))]
    nuc$aln <- rep(1:nrow(coords), nc)
    
    if (fillGaps) {
        n <- 2:nrow(nuc)
        newAln <- nuc[, c(T, aln[n] != aln[n - 1L], T)]
        
        fillG <- function(pos) {
            inGap <- pos == 0L
            startGap <- which(inGap[n] & (!inGap[n - 1L] | newAln[n])) + 1L
            endGap <- which(inGap[n - 1L] & (!inGap[n] | newAln[n]))
            
            if (tail(inGap, 1)) {
                endGap <- c(endGap, max(n))
            }
            
            previousPos <- pos[startGap - 1L]
            
            if (startGap[1] == 1L) {
                previouPos <- c(0L, previousPos)
            }
            
            nextPos <- pos[endGap + 1L]
            
            if (is.na(tail(nextPos, 1))) {
                nextPos[length(nextPos)] <- 0L
            }
            
            toTake <- ifelse((previousPos < nextPos &
                !newAln[startGap]) | newAln[endGap + 1L],
            previousPos,
            nextPos
            )
            
            pos[inGap] <- rep(toTake, endGap - startGap + 1L)
            
            pos
        }
        
        nuc[, c("pos1", "pos2") := .(fillG(pos1), fillG(pos2))]
    }
    
    nuc
}



splitStringInParts <- function(string, size) {
    # splits a character string into a character vector whose word length is defined by "size" (an integer)
    
    pat <- paste0(".{1,", size, "}")
    res <- stri_extract_all_regex(string, pat)
    names(res) <- names(string)
    res
}


aaToCDS <- function(aas, cds, collapse = T) {
    # generates a codon alignment based on a protein alignment and the
    # corresponding cds, both being character vectors. Returns a character vector
    # is collapse is TRUE. If collapse is FALSE, keeps codon separate and returns
    # a list of vector of codons. This function does NOT check whether the proteins
    # actually correspond to the cds according to a genetic code
    
    if (any(stri_length(gsub("-", "", aas, fixed = T)) != stri_length(cds) /
        3)) {
        stop("CDS length must be exactly 3 times protein length")
    }
    
    codons <- splitStringInParts(cds, 3)
    aas <- strsplit(aas, "", fixed = T)
    
    fillCodons <- function(aas, codons) {
        res <- character(length(aas))
        deleted <- aas == "-"
        res[!deleted] <- codons
        res[deleted] <- "---"
        res
    }
    
    res <- Map(fillCodons, aas, codons)
    
    if (collapse) {
        res <- sapply(res, stri_flatten)
    }
    
    res
}

subSeq <- function(string, start, end, reverse = T) {
    # builds upon Biostrings' subseq() to automatically reverse complement the
    # sequences (by default) if end coordinates are lower than starts. string is an
    # XStringSet, and start and end are integer vectors. Works like substring()
    
    f <- start > end
    temp <- end
    end[f] <- start[f]
    start[f] <- temp[f]
    subs <- subseq(string, start, end)
    
    if (reverse & any(f)) {
        subs[f] <- reverseComplement(subs[f])
    }
    
    subs
}


seqinrAlignment <- function(seqs, names = NULL) {
    # generates a seqinr alignment object from a character vector
    
    if (is.null(names)) {
        if (!is.null(names(seqs))) {
            names <- names(seqs)
        } else {
            names <- paste(rep("seq", length(seqs)), 1:length(seqs), sep = "_")
        }
    }
    
    if (var(nchar(seqs)) > 0) {
        stop("sequence lengths differ")
    }
    
    aln <- list(
        nb = length(seqs),
        nam = names,
        seq = as.character(seqs),
        com = NA
    )
    
    class(aln) <- "alignment"
    
    aln
}

clusterFromPairs <- function(V1,
                             V2,
                             reciprocal = F,
                             int = F) {
    # assigns elements to clusters via single-linkage clustering, based on pairs
    # (e.g., query - subject). Pairs are elements of V1 and V2 (vectors of equal
    # length) at the same position If an element is found in several pairs (and
    # clusters, iteratively), the elements of these pairs (or clusters) will be
    # assigned to the same cluster
    
    if (length(V1) != length(V2)) {
        stop("provide vectors of equal length")
    }
    
    m <- 0
    if (!int) {

        # we assign unique cluster ids (integers) to elements (faster than union())
        uniq <- unique(c(V1, V2))
        
        if (is.character(V1)) {
            V1 <- chmatch(V1, uniq)
        } else {
            V1 <- match(V1, uniq)
        }
        
        if (is.character(V2)) {
            V2 <- chmatch(V2, uniq)
        } else {
            V2 <- match(V2, uniq)
        }

        # will give the correspondance between current cluster ids (the indices of m) and new 
        # clusters ids (values of m), for replacement. Initially, they are the same
        m <- 1:length(uniq)
    } else {
        m <- 1:max(V1, V2)
    }
    
    if (!reciprocal) {
        # if there is no reciprocity in pairs

        # we now have reciprocal pairs of clusters
        pairs <- data.table(q = c(V1, V2), s = c(V2, V1))
        
    } else {
        pairs <- data.table(q = V1, s = V2)
    }
    
    f <- pairs[, q > s]
    
    while (any(f)) {
        # clustering proceeds until the 2 elements of every pair are assigned to the same cluster

        # for a given query cluster we need lowest cluster id in all pairs where it is found, including itself.
        # So there's no need to use rows for which the cluster is lower than the query
        mins <- pairs[f, min(s), by = q]

        # and we replace each cluster id by this lowest id in both query and subject clusters 
        # (these 2 lines are much faster than a match())
        m[mins$q] <- mins$V1
        pairs[, c("q", "s") := .(m[q], m[s])]
        f <- pairs[, q > s]
        cat("*")
    }

    # we get cluster ids corresponding to V1 elements since they are the same for V2
    res <- pairs$q[1:length(V1)]

    # this removes possible gaps in cluster numbering
    match(res, unique(res))
}

pastePair <- function(string1, string2, sep = "_") {
    # paste elements of string1 and string2 together (character vectors) so that
    # the lower of the two (in alphabetical order) is at the left of the pasted
    # result
    
    f <- string1 < string2
    
    if (!is.character(string1)) {
        string1 <- as.character(string1)
        string2 <- as.character(string2)
    }
    
    ifelse(f,
        stri_join(string1, string2, sep = sep),
        stri_join(string2, string1, sep = sep)
    )
}

allPairs <- function(v,
                     sort = T,
                     noDups = T,
                     reciprocal = F,
                     same = F,
                     v2 = NA) {
    # returns all possible pairs for elements of a vector, or two vectors (2dn
    # vector specified in v2). Much faster than combn(v, 2). reciprocal = T means
    # reciprocal pairs are returned. sort = T means pairs are sorted alphabetically
    # or numerically. same tells if pairs comprising the same element twice should
    # be returned (only valid if v2 = NA)
    
    if (noDups) {
        v <- unique(v)
    }
    
    v2 <- unique(v2)
    if (sort) {
        v <- sort(v)
    }
    
    v2 <- sort(v2)
    
    if (all(is.na(v2)) & !reciprocal) {
        i <- 0L
        
        if (!same) {
            i <- 1L
        }
        
        lv <- length(v)
        V2 <- unlist(lapply((i + 1):lv, function(x) {
            v[x:lv]
        }))
        
        return(data.table(V1 = rep(v, lv:1 - i), V2))
    }
    else {
        self <- F
        if (all(is.na(v2))) {
            v2 <- v
        }
        self <- T
        cb <- data.table(
            V1 = rep(v, each = length(v2)),
            V2 = rep(v2, length(v))
        )
        if (!self & reciprocal) {
            cb <- rbind(cb, cb[, 2:1])
        }
        cb
    }
}


data.tableFromCommunities <- function(comm) {
    # creates a data table from communities generated by iGraph
    
    cNumbers <- unlist(Map(rep, 1:length(comm), sapply(comm[1:length(comm)], length)))
    members <- unlist(comm[1:length(comm)])
    data.table(member = as.integer(members), community = cNumbers)
    
}



combineHits <- function(blast,
                        maxDist = 50,
                        maxOverlap = 30,
                        maxDiff = 50,
                        blastX = F) {
    # combines nearby HSPs in a bast tabular output, if these HSPs involve the same
    # query and subject. Do not combine HSPs that are distant by more than maxDist
    # and overlap by more than maxOverlap. bastX should be TRUE if hits are
    # obtained by blastx
    
    blast <- blast[, 1:12, with = F]
    setnames(
        blast,
        c(
            "query",
            "subject",
            "identity",
            "length",
            "mismatches",
            "indels",
            "qStart",
            "qEnd",
            "sStart",
            "sEnd",
            "evalue",
            "score"
        )
    )
    
    blast[evalue > 1, evalue := 1]
    
    if (blastX) {
        blast[, c("sStart", "sEnd") := .(sStart * 3, sEnd * 3)]
    }
    
    blast[, pair := stri_c(query, subject)]
    blast[sEnd < sStart, c("sStart", "sEnd") := .(-sStart, -sEnd)]
    blast[qEnd < qStart, c("qStart", "qEnd") := .(-qStart, -qEnd)]
    setorder(blast, pair, qStart, sStart)
    rows <- 2:nrow(blast)
    samePair <- c(F, blast[, pair[rows] == pair[rows - 1]])
    distQuery <- c(0, blast[, qStart[rows] - qEnd[rows - 1]])
    distQuery[!samePair] <- NA
    distSubject <- c(0, blast[, sStart[rows] - sEnd[rows - 1]])
    distSubject[!samePair] <- NA
    
    w <- which(
        distSubject <= maxDist &
            distQuery <= maxDist &
            distSubject >= -maxOverlap &
            distQuery >= -maxOverlap & abs(distSubject - distQuery) <= maxDiff
    )
    
    if (length(w) == 0) {
        print("no HSP can be combined")
        return(blast[, -13, with = F])
    }
    
    c <- clusterize(w, 1)
    blast[w, group := c]
    blast[setdiff(w - 1, w), group := (unique(c))]
    blast[, c("qStart", "qEnd", "sStart", "sEnd") := .(abs(qStart), abs(qEnd), abs(sStart), abs(sEnd))]

    queryRegions <- combineRegions(blast[!is.na(group), 
                                         data.frame(group, qStart, qEnd)], distance = maxDist * 2)
   
     subjectRegions <- combineRegions(blast[!is.na(group), 
                                           data.frame(group, sStart, sEnd)], distance = maxDist * 2)
   
      regionStats <- blast[!is.na(group), .(
        query = query[1],
        subject = subject[1],
        identity = mean(identity),
        length = sum(length),
        mismatches = sum(mismatches),
        indels = sum(indels) + .N - 1,
        evalue = prod(evalue),
        score = sum(score) - .N * 10
    ), by = group]

    combined <- data.table(regionStats[, -1, with = F], queryRegions[, -1,
        with = F], subjectRegions[, -1, with = F])
    
    setcolorder(combined, c(1:6, 9:12, 7:8))
    setnames(combined, names(blast)[1:12])
    
    all <- rbind(combined, blast[is.na(group), -(13:14), with = F])
    
    if (blastX) {
        all[, c("sStart", "sEnd") := .(sStart / 3, sEnd / 3)]
    }
    
    all
}

removeReciprocal <- function(blast, removeSelf = T) {
    # removes reciprocal and self hits (if removeSelf = T) from blast results, format 6. Score is expected in column 12.
   
     blastb <- copy(blast)
    names <- colnames(blast)
    setnames(blastb, stri_c("V", 1:ncol(blastb)))
    blastb[, pair := pastePair(V1, V2)]
    setorder(blastb, -V12)
    f <- T
    
    if (removeSelf) {
        f <- blastb[, V1 != V2]
    }
    
    blastb <- blastb[!duplicated(pair) & f]
    blastb[, pair := NULL]
    setnames(blastb, names[1:ncol(blastb)])
    blastb
}

cLinkFromPairs <- function(V1, V2, linked) {
    # makes cliques (complete-linkage clusters) from pairs of element (V1 and V2
    # are vectors of equal lengths). linked is a logical that tells which pairs of
    # elements are linked. The algorithm is described in peccoud et al. 2017 PNAS
    
    u <- unique(c(V1, V2))
    
    if (length(V1) != length(V2)) {
        stop("provide vectors of equal length")
    }
    
    if (is.integer(V1)) {
        i1 <- match(V1, u)
    } else {
        i1 <- chmatch(V1, u)
        # convertion to integer number
    }
    
    if (is.integer(V2)) {
        i2 <- match(V2, u)
    } else {
        i2 <- chmatch(V2, u)
    }

    # for each element number (index of clique), the element(s) belonging to its clique (values of gr). Originally, an element is just grouped with iself, so this integer vector is suited. But it will become an integer list
    clique <- 1:length(u)

    # all pairs of adjacent elements, with reciprocity.
    allReciprocalPairs <- cbind(c(i1[linked], i2[linked]), c(i2[linked], i1[linked]))

    # we add to these pairs all elements paired with themselves
    allReciprocalPairs <- rbind(allReciprocalPairs, cbind(clique, clique))

    # so the splitting yields, for each element, a "set" (=value of the list) of all the elements
    set <- split(allReciprocalPairs[, 1], allReciprocalPairs[, 2])
    # that are linked to it (and are potentially allowed to be grouped with it), including itself.
    # Because all element numbers were used to split that list, these sets are ordered by element number (set[[n]] is the set for element number n)
    # so we do not need to fetch a set by name matching => :-)
    set <- sapply(set, unique)
    
    pb <- txtProgressBar(
        char = "*",
        width = 100,
        style = 3
    )
    
    n <- 1
    w <- which(linked)
    
    for (i in w) {
        # for each pair of adjacent elements,

        # gets left element number
        t1 <- i1[i]

        # right element number
        t2 <- i2[i]
        
        if (!t2 %in% clique[[t1]]) {
            # if the two are not already grouped
            if (all(clique[[t2]] %in% set[[t1]]) &
                all(clique[[t1]] %in% set[[t2]])) {
                # and if elements of their clique (which includes themselves) are in the set of allowed elements of the other	element

                # we can merge the t1 and t2 cliques into a new clique
                newGroup <- union(clique[[t1]], clique[[t2]])

                # and this new clique becomes the clique of all elements composing it (it's like a square)
                clique[newGroup] <- list(newGroup)

                # then each element of the new clique has its set become the intersection of the t1 and t2 sets, i.e,
                set[newGroup] <- list(intersect(set[[t2]], set[[t1]]))
                # it now contains only elements that are adjacent to BOTH t1 and t2.
                # Other elements cannot be admitted in the clique as it would break the rule. Note that all sets for elements of the same clique become the same.
            }
        }
        setTxtProgressBar(pb, n / length(w))
        n <- n + 1
    }
    cliqueElements <- unique(lapply(clique, sort))

    # makes the correspondance between a element number (vector position) and a clique number (vector value)
    cliqueNumber <- rep(1:length(cliqueElements), sapply(cliqueElements, length))
    names(cliqueNumber) <- u[unlist(cliqueElements)]
    cliqueNumber
}



seqtk <- function(fas,
                  bed,
                  out,
                  ex = "seqtk subseq",
                  formated = F) {
    # calls seqtk to return sequence from a fasta fas, based on bedfile bed.
    # Frite to fasta file out (if specified) or return to R (if not). 
    #if formatted is TRUE, return sequences as a DNAStringset (if not, returns a character vector 
    # corresponding to the fasta) ex specifies how seqtk subseq should be executed
    
    if (missing(out)) {
        seqs <- system(paste(ex, fas, bed), intern = T)
        
        if (formated) {
            f <- stri_sub(seqs, 1, 1) == ">"
            dt <- data.table(content = seqs[!f], id = cumsum(f)[!f])
            concat <- DNAStringSet(dt[, stri_flatten(content), by = id]$V1)
            names(concat) <- stri_sub(seqs[f], 2L, nchar(seqs[f]))
            return(concat)
            
        } else {
            return(seqs)
        }
        
    } else {
        system(paste(ex, fas, bed, ">", out))
    }
}



# FUNCTION THAT ARE MORE SPECIFIC TO THE HTT ANALYSIS ----------------------------------------------------------

extractSpeciesNames <- function(x) {
    # extracts species names from file names (character vector x) used at various
    # steps of the analysis. genus and species in x must be separated by
    # underscores
    
    stri_extract_last(x, regex = "[A-Z][a-z]+_[a-z]+[_]*[a-z]*")
}


extractCopies <- function(sp, gff, cons, genome, out, minLen = 300L) {
    # extracts TE copies on length ≥ minLen of a given species (sp) from tis genome (genome) to gzipped fasta file (out), 
    #based on a gzipped repeat masker gff (gff), write also a report of the percentage of consensus sequences in "cons" 
    # (the repeat modeler fasta file) that have masked copies as a filename ending by ".percentMasked.txt".
    # this function returns metrics on the TE compositions as a data.table (see end of function)
    
    gff <- fread(gff,
        sep = "\t",
        header = F,
        skip = 1
    )

    # extracts TE family name, and super family name
    gff[, family := stri_extract_first(V9, regex = "rnd-[^ ]+")]
    gff[, superF := gsub(
        "Class=",
        "",
        stri_extract_first(V9, regex = "Class=[^;]+"),
        fixed = T
    )]

    gff[, family := stri_c(family, "#", superF)]

    # we will compute the proportion of consensuses that have masked sequences (to check if repeat masker worked properly)

    # imports consensus (=family) names from the fasta file generated by repeat modeler (faster than readDNAStringset)
    cons <- system(paste("grep '>'", cons), intern = T)
    cons <- gsub(">", "", splitToColumns(cons, " ", 1))
    prop <- mean(cons %chin% gff$family)
    
     # selects TE copies that are long enough and assigned to a precise super family (which has a slash in its name, like "DNA:Tc1/Mariner")
    selectedCopies <- gff[grepl("/", family, fixed = T) & V5 - V4 + 1L >= minLen]

    # temorary bed file used by seqtk in the extraction process
    bed <- stri_c(sp, ".TEcopies.bed")

    # writes TE coordinates in the bed used by seqtk (not that we correct start coordinates to match the system used by seqtk)
    writeT(selectedCopies[, .(V1, V4 - 1L, V5)], bed, col.names = F)

    # imports sequences directly from seqtk (as a fasta string)
    seqs <- seqtk(genome, bed)

    # locates sequences names
    f <- stri_sub(seqs, 1, 1) == ">"
    if (sum(f) != nrow(selectedCopies)) {
        stop(paste("missing sequences", sp))
    }

    # extract these names
    seqNames <- stri_extract_first(seqs[f], regex = ">[^ ]+")

    # recreates sequences names that are generated by seqtk, to find were the copies are 
    # in this fasta (because the order is not that in the bed, it's in the genome order)
    selectedCopies[, seqName := stri_c(">", V1, ":", V4, "-", V5)]

    # and creates copy names based on species, contig, coordinates, direction and TTEe family (also containing super family)
    newNames <- selectedCopies[chmatch(seqNames, seqName), stri_c(">", stri_c(sp, V1, V4, V5, V7, family, sep = ":"))]

    # replace sequence names in the fasta with these new names, then writes to gzipped fasta file
    seqs[f] <- newNames
    file.remove(bed)
    gz1 <- gzfile(out, "w")
    writeLines(seqs, gz1)
    close(gz1)
    
    # writes a report about the percent of consensuses that masked TEs
    writeT(
        data.table(sp = sp, prop, nCopies = nrow(selectedCopies)),
        stri_c("TEs/TEcomposition/", sp, ".percentMasked.txt")
    )

    # we return metrics on the TE composition of the genome, per superfamily, for selected copies
    # turns family names to integer numbers, to speed up the command below

    selectedCopies[, family := toInteger(stri_c(family, "#", superF))]
    
    res <- selectedCopies[, .(
        nCopies = .N,
        nFam = length(unique(family)),
        bp = sum(V5 - V4 + 1L)
    ), by = superF]

    # adds species name in a new column
    res[, sp := sp]
    res
}

getBUSCOs <- function(folder) {
    # concatenates busco genes of a species (having a specific folder) in 
    # single fasta and renames sequences by adding the species name to each
    
    cds <- list.files(folder, pattern = ".fna", full.names = T)

    # readDNAStringSet takes a vector of file paths, but lapply was 
    # necessary to avoid a bug due to too many open connections
    seqs <- lapply(cds, readDNAStringSet)

    # concatenating sequences (unlist doesn't reliably work on DNAStringSets)
    concat <- do.call("c", lapply(seqs, as.character))

    # because somehow the c() function above would add the object number in the list to all sequence names
    nams <- unlist(lapply(seqs, names))
    sp <- extractSpeciesNames(folder)

    # renaming sequences
    names(concat) <- stri_c(stri_extract_first(nams, regex = "[^:]+:"), sp)
    
    # and writing them to disk
    writeXStringSet(DNAStringSet(concat), stri_c("CDS/", sp, ".CDS.fas"))
    
    print(paste("done", sp))
}


translateCDS <- function(cds, aa) {
    # translates CDS into proteins, cds is the single path to CDS fasta, and aa is the path to output protein fasta
    # and returns some reports regarding number of sequences with internal stop codons
    
    seq <- readDNAStringSet(cds)

    # checks for final stop codons
    f <- lastChars(seq, 3) %chin% c("TAA", "TAG", "TGA")

    # if present, removes them
    seq[f] <- stri_sub(seq[f], 1, nchar(seq[f]) - 3)

    # does the translation
    tr <- suppressWarnings(translate(seq, if.fuzzy.codon = "solve"))

    # find if stop codons are present
    stops <- stri_count(tr, fixed = "*")

    # we only write the sequences without stop codons
    writeXStringSet(tr[stops == 0L], aa)
    
    # we return a table on the presence of stop codons for each CDS
    data.table(
        cds = basename(cds),
        nseq = length(seq),
        stops = sum(stops > 0)
    )
}


seqtkDiamond <- function(fas, bed, db, out) {
    # extract sequences from fasta file fas using bedfile bed and blasts them search against diamond protein database db, to output file out
    
    system(
        paste(
            "seqtk subseq",
            fas,
            bed,
            "| diamond blastx -p 1 --sensitive -d",
            db,
            "--quiet -k 1 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen -o",
            out
        )
    )
    
    if (file.size(out) > 0) {
        # because fread() cannot read empty files
        return(fread(out, header = F, sep = "\t"))
        
    } else {
        # returns an empty table if there is not hit
        return(data.table())
    }
}

copyName <- function(x) {
    # in certain files (fastas, mostly), TE copy names must have species and TE
    # family names attached. This function removes them to match copy names in
    # certain tables where species and family name are in other columns

    # fieles are separated by colons. We split them
    mat <- stri_split(x, fixed = ":", simplify = T)

    # extacts the fields we want and concatenates them back
    stri_c(mat[, 2], mat[, 3], mat[, 4], mat[, 5], sep = ":")
}


expectedProtein <- function(query, subject) {
    # determines whether a TE copies has a hit against a repBase protein of the same super family. 
    # This is a bit tricky since the super family information has no standard format in repbase

    # gets superfamily name from query name

    superFam <- toupper(gsub("?", "", splitToColumns(query, "#", 2), fixed = T))

    # we split superfamily names into elements, since they are classes, subclasses....
    superFam <- stri_split_regex(superFam, "[:punct:]", simplify = T)
    superFam[superFam == ""] <- NA

    # same for the subject protein
    superFamProt <- toupper(gsub("?", "", splitToColumns(subject, "#", 2), fixed = T))
    superFamProt <- stri_split_regex(superFamProt, "[:punct:]", simplify = T)
    test <- rowSums(superFam == superFamProt[, 1:min(ncol(superFam), ncol(superFamProt))],
        na.rm =
            T

        # determines the number of elements that are identical between superFamily names (e.g., DNA + TCMAR + Mariner)
    )

    # > 0 if the number of identical elements is lower than the number of elements
    diff <- rowSums(!is.na(superFam)) - test
    exp <- integer(length(test))

    # set to 2 when all elements are the same
    exp[diff == 0] <- 2
    exp[diff > 0 &

        # set to 1 when 2 elements are the same, but not all (ex DNA-TCMAR-MARINER vs DNA-TCMAR-TC1 or just DNA-TCMAR).
        test >= 2] <- 1
    exp
}


processBlastX <- function(blast) {
    # processes results (tabular file of hits) of a blastx of TE copies against repeat proteins

    # in case copies subparts were blasted, we convert qstart and qend coordinates into initial copy coordinates, 
    # based on where the blasted part starts in the copy (field 8 of query name, as added by seqtk)
    # for this, we need to split the query names into fields
    fields <- as.data.frame(splitToColumns(blast$V1, ":"))
    
    # detects whether copy parts and not full copies (first round) were blasted, based on the number of fields
    if (ncol(fields) == 7) {
       
        # we split two numbers separated by a dash: the start and end of the copy part that was used. We get the start
        ps <- as.integer(splitToColumns(fields[, 7], "-", 1))

        # if NA (no field), it means that we blasted the full copy, so the start position is 1
        partStarts <- ifelse(is.na(ps), 1L, ps)

        # converting start and end query hsp coordinates into original copy coordinates
        blast[, c("V7", "V8") := .(V7 + partStarts - 1L, V8 + partStarts - 1L)]

        # replacing copy part names by copy names, thus ignoring fields 8 and 9
        blast[, V1 := do.call(stri_c, c(fields[, 1:6], sep = ":"))]
    }

    # we now combine successive HSPs related to the same hit, since blastx creates a new HPS where there is a frame shift
    blast <- combineHits(blast, blastX = T)
    # since blastX finds many spurious hits and since valid hits may also have low score, we 
    # determine if the protein belongs to the same super family as the copy (suggesting a valid hit)
    blast[, ex := expectedProtein(query, subject)]
    blast
}


unalignedParts <- function(blast) {
    # this extracts copy parts that did not align to proteins from the current blastx round, to be aligned with the next blastx round.
    
    # we do it for copies that had an acceptable HSP (=expected protein family and at least 30 amino-acids aligned), 
    #since any other hit should be of lower quality (we set max_taget_seqs 1, so we have the best hit possible).
    selectedCopies <- blast[ex > 0 &  length >= 30, unique(query)]
    
    # in row passes our criteria, we return an empty table
    if (length(selectedCopies) == 0) {
        return(data.table())
    }

    # we retrieve the alignment coordinates for these selected copies, in all blastx rounds to avoid blasting 
    # the same parts several time, which is why we needed to concatenate all blastx rounds in the previous function
    coords <- blast[query %in% selectedCopies, .(query, qStart, qEnd)]

    # making sure qStart < qEnd
    coords[qStart > qEnd, c("qStart", "qEnd") := .(qEnd, qStart)]

    # getting the overall start and end of all HSPs (combined) for a copy. 
    # We're not interested in parts between HSPs on the same protein since they should not encompass a different protein
    coords <- coords[, .(start = min(qStart), end = max(qEnd)), by = query]

    # retrieving copy length from start and end positions in contigs = fields 3 and 4 of copy names
    copyCoords <- splitToColumns(coords$query, ":", 3:4, mode = "integer")

    # this will be used as the end coordinate of each trailing unaligned part
    copyEnd <- as.integer(copyCoords[, 2]) - as.integer(copyCoords[, 1]) + 1
    
    # the unaligned leading part of each copy
    leading <- coords[, data.table(
        copy = query,
        start = 0,
        end = start - 1
    )]
    
    # and the unaligned trailing part
    trailing <- coords[, data.table(
        copy = query,
        start = end,
        end = copyEnd
    )]
    
    #we stack them in a single table
    parts <- rbind(leading, trailing)

    # returns parts that are at least 100-bp long
    parts[end - start >= 100]
}


nestedClades <- function(tree, symmetrical = T) {
    # from a tree ("phylo" class), returns a matrix were the value is TRUE if the
    # clade at the column (the "child") is the same as, or nested in, the clade at
    # the row (the "parent"). If symmetrical is TRUE, the matrix is also TRUE if
    # the clade at the row is nested within the clade at the column. Uses the fact
    # that nodes of a phylo object are integer numbers
    
    nodes <- unique(as.vector(tree$edge))
    
    childs <- lapply(nodes, function(x) {
        subNodes(tree, x)
    })
    
    childOf <- cbind(
        parent = rep(nodes, sapply(childs, length)),
        child = unlist(childs)
    )
    
    mat <- matrix(F, max(childOf), max(childOf))
    diag(mat) <- T
    mat[childOf] <- T
    
    if (symmetrical) {
        mat[childOf[, 2:1]] <- T
    }
    
    mat
}

MRCA <- function(tree, tips) {
    # the same as getMRCA function from ape, except it returns the tip if only one tip is specified, instead of NULL
    
    if (length(tips) == 1) {
        return(match(tips, tree$tip.label))
    }
    getMRCA(tree, tips)
}


nodeDepth <- function(tree, nodes = sort(unique(as.vector(tree$edge)))) {
    # returns the age, or depth, of nodes of a timetree (nodes specified as integers)
    
    len <- node.depth.edgelength(tree)
    max(len) - len[nodes]
}



tipsForNode <- function(tree,
                        node,
                        names = F,
                        itself = F) {
    # returns all tips for a node from tree, or NULL if node is a tip (except if itself is TRUE). 
    # Returns tip names if names is TURE, else tip numbers
    
    tips <- NULL
    n <- node
    
    while (length(n) > 0) {
        descendants <- tree$edge[tree$edge[, 1] %in% n, 2]
        tips <- c(tips, descendants[descendants <= length(tree$tip)])
        n <- setdiff(descendants, tips)
    }
    
    if (itself & length(tips) == 0) {
        tips <- node
    }
    
    if (names) {
        return(tree$tip.label[tips])
    }
    
    tips
}


tipsForNodes <- function(tree,
                         nodes,
                         names = F,
                         itself = T) {
    # returns all tips of several tree nodes, based on the function above. 
    # returns a data.table whose firt column is the node number and second the tips
   
     tips <- lapply(nodes, function(x) {
        tipsForNode(tree, x, names = names, itself = itself)
    })
     
    nSP <- sapply(tips, length)
    dt <- data.table(node = rep(nodes, nSP), tip = unlist(tips))
    
    dt[order(-rep(nSP, nSP))]
}

subNodes <- function(tree, node) {
    # returns all subnodes of a tree node, including tips
    
    descendants <- node
    subNodes <- NULL
    
    repeat {
        descendants <- tree$edge[tree$edge[, 1] %in% descendants, 2]
        
        if (length(descendants) == 0) {
            break
        }
        
        subNodes <- c(subNodes, descendants)
    }
    
    subNodes
}


cladesOfAge <- function(tree,
                        age,
                        withTips = F,
                        names = T,
                        singletons = T) {
    # returns the oldest clades of age ≤ age, from a tree. None of these clade should be nested within another
   
    ages <- nodeDepth(tree)

    # all nodes younger than the requested age
    clades <- which(ages <= age)

    # to obtain children of nodes
    edges <- as.data.table(tree$edge)

    # removes clades that are children of nodes younger than the age. This leaves only the oldest clade of age ≤ minAge
    clades <- setdiff(clades, edges[V1 %in% clades, ]$V2)
    
    if (!singletons) {
        clades <- setdiff(clades, 1:length(tree$tip.label))
    }
    
    if (!withTips) {
        return(clades)
    }
   
     tipsForNodes(tree, clades, names = names)
}


combineHomologous <- function(protRegions, blastp, protSeqs) {
    # combines homologous protein regions among copies in a community. Several
    # similar copies may have a hit to different proteins, which are homologous. We
    # would rather have a single protein than several homologous ones.
    n <- 1

    # we do it iteratively, this is the first round
    oneProt <- NULL
    repeat {
        print(paste("round", n, nrow(protRegions), "regions"))
        n <- n + 1

        # retrieving information about community from fields in protein names (the column is called "sequenceID" due to how combineRegions() works)
        protRegions[, c("commID", "prot") := as.data.table(splitToColumns(sequenceID, " ", 1:2))]
        protRegions[, commID := as.integer(commID)]
        # for each community, we find regions belonging to different, homologous,
        # proteins. We generate all possible pairs of regions in a community. Based
        # on the protein homology, we convert start/end region coordinates of the
        # less represented protein in the community (between the two) into
        # coordinates of the other protein. Then the converted region, which also
        # takes the other protein name, may be combined with other regions of this
        # protein with combineRegions()
        nProt <- protRegions[, length(unique(prot)), by = commID]

        # we store protein regions that belong to communities associated to a single protein
        oneProt <- rbind(oneProt, protRegions[commID %in% nProt[V1 == 1L, commID]])
        protRegions <- protRegions[commID %in% nProt[V1 > 1L, commID]]

        # row identifier for each region, used to generate pairs of regions to check homology between proteins, as below
        protRegions[, row := 1:.N]

        # for each tranfer, we generate all pairs of protein regions (referenced by their rows)
        allP <- protRegions[, data.table(allPairs(row, sort = F)), by = commID]
        
        # gets  protein names and region coordinates for these pairs
        allP <- allP[, data.table(
            commID, 
            protRegions[V1, .(prot, start, end)], protRegions[V2, .( prot2 =  prot, start2 = start, end2 = end)]
        )]
        
        # if protein regions belong to the same protein, there is no coordinate to convert
        allP <- allP[prot != prot2]

        # for 2 homologous regions, we will convert the one belonging to the protein
        # that is less represented in the community

        # the non overlapping length of a given protein in copies in a community
        lenPerCom <- protRegions[, sum(end - start + 1), by = sequenceID]
        allP[, l1 := lenPerCom[chmatch(paste(commID, prot), sequenceID), V1]]
        allP[, l2 := lenPerCom[chmatch(paste(commID, prot2), sequenceID), V1]]


        # puts data from the less represented protein (whose region coordinates will be converted) at the left
        allP[l1 > l2, c("prot", "prot2", "start", "end", "l1", "l2") := .(prot2, prot, start2, end2, l2, l1)]

        # creates an identifier for protein pairs and removes now uncessary columns at the same time
        allP[, c("pair", "start2", "end2") := .(paste(prot, prot2), NULL, NULL)]

        # we use this id to discard protein pairs not in the blastp as not homologous
        allP <- allP[pair %chin% blastp[, pair]]

        # we will assess how much a region to be converted is covered by the blastp 
        # HSP with the homologous protein (prot2), so we retrieve HSP coordinates.
        allP[, c("qStart", "qEnd") := blastp[match(allP$pair, pair), .(qStart, qEnd)]]
       
        # cov = the proportion of a region that is covered by the alignment with prot2
         allP[, cov := intersection(start, end, qStart, qEnd, negative = F) / (end - start + 1)]

        # we select region start/end coordinates that are inside (or close to) the HSP (else coordinate conversion would be risky
        sel <- allP[cov >= 0.9, .(commID, prot, start, end, prot2, l2, pair)]

        # a given protein region may have several homologous proteins in a community.
        # We select the one with highest presence (l2) in the community, as below
        sel[, region := paste(commID, prot, start, end)]
        setorder(sel, region, -l2)
        sel <- sel[!duplicated(region)]

        # we align homologous proteins we've selected, to obtain precise coordinate correspondance.
        # we get the blast HSP for protein pairs we need to align
        # (we may align the same pair twice reciprocally, but we can afford it, as it simplifies the code)
        f <- blastp[, pair %chin% sel[, pair]]
       
        # if there is less than 200 protein to convert, we stop
        if (sum(f) < 200) {
            break
        }
       
        # parts of protein sequences involved in HSPs. We don't use subseq() as the alignment has issues with AAStringSets
        qSeqs <- blastp[f, stri_sub(protSeqs[query], qStart, qEnd)]
        sSeqs <- blastp[f, stri_sub(protSeqs[subject], sStart, sEnd)]
        cat("aligning...")
        aln <- alignWithEndGaps(qSeqs, sSeqs)

        cat("done.\n")

        # splits the alignment into nucleotides and assigns coordinates to each aligned position, using HSP coordinates
        corres <- splitAlignment(aln, blastp[f, .(query, subject, qStart, qEnd, sStart, sEnd)], fillGaps = T)

        # selects the columns we need and names them more appropriately
        corres <- corres[, .(prot = seq1, prot2 = seq2, pos1, pos2)]

        # gets start (min) and end (max) coordinates of the left protein in each alignment
        ends <- corres[, .(mini = min(pos1), maxi = max(pos1)), by = .(prot, prot2)]
        ends[, pair := paste(prot, prot2)]

        # these indicate the minimal and maximal coordinate that is present in corres for the protein, since we allowed starts 
        #and ends of regions to be sighly outside HSPs (cov >= 0.9). 
        #We cannot just fetch starts and ends as some are not present in corres
        sel[, c("mini", "maxi") := ends[match(sel$pair, pair), .(mini, maxi)]]

        # if the start and end are within these mini and maxi, we can fetch them directly, 
        # or else we use mini and maxi (which should be close)
        sel[, c("mini", "maxi") := .(pmax(mini, start), pmin(maxi, end))]

        # id for each particular position of a region in each alignment
        corres[, id := stri_c(prot, prot2, pos1, sep = " ")]

        # which we use to fecth start and end positions to convert in each protein region
        sel[, c("startID", "endID") := .(paste(prot, prot2, mini), paste(prot, prot2, maxi))]
        
        sel[, c("newStart", "newEnd") := .(corres[chmatch(startID, id), pos2] - (mini - start), 
                                           corres[chmatch(endID, id), pos2] + end - maxi)]

        # negative starts may occure if the start of the region was lower than mini 
        # and the homologous protein was shorter at the start. We correct that
        sel[newStart < 1, newStart := 1L]

        # for similar reasons, the new end position may be longer than the homologous protein.
        sel[newEnd > nchar(protSeqs[prot2]), newEnd := nchar(protSeqs[prot2])]

        # now update these positions in the original data table or protein regions.
        protRegions[, region := paste(sequenceID, start, end)]
        protRegions[chmatch(sel$region, region), c("prot", "start", "end") := sel[, .(prot2, newStart, newEnd)]]

        # re-combining new regions
        protRegions <- combineRegions(protRegions[, .(paste(commID, prot), start, end)], distance = 10L)
        rm(corres, ends, aln)
    }
    rbind(protRegions[, -"row"], oneProt)
}

proteinOverlap <- function(regionComp, blastp, protSeqs) {
    # determines how much TEs from different communities overlap in respect to the protein-coding regions they cover
    
    regionComp[, protPair := paste(prot, prot2)]

    # retrieves HSP coordinates of protein pairs (including self-alignments)
    regionComp[, c("qStart", "qEnd", "sStart", "sEnd", "pID", "length") := blastp[match(protPair, pair), .(qStart, qEnd, sStart, sEnd, pID, length)]]

    # intersection of left region with the HSP (returned as a range)
    inter1 <- regionComp[, intersection(start, end, qStart, qEnd, range = T)]

    # intersection of right region
    inter2 <- regionComp[, intersection(start2, end2, sStart, sEnd, range = T)]
    # now converting the intersection coordinates of one protein into the coordinate system of the other protein to compute the overlap. 
    # This requires properly aligning the proteins and see which position corresponds to which in the other protein
    
    # we only realign different proteins (not the same with itself) that are involved in blastp
    toConvert <- regionComp[, !is.na(sStart) &
        prot != prot2 &
        !is.na(inter1[, 1]) & !is.na(inter2[, 1])]
    
    realign <- blastp[pair %chin% regionComp[toConvert, pastePair(prot, prot2, sep = " ")]]
    
    aln <- realign[, alignWithEndGaps(
        stri_sub(protSeqs[query], qStart, qEnd),
        stri_sub(protSeqs[subject], sStart, sEnd)
    )]
    
    nuc <- splitAlignment(aln, realign[, .(query, subject, qStart, qEnd, sStart, sEnd)], fillGaps = T)
    
    nuc <- rbind(nuc, nuc[, data.table(
        base1 = base2,
        base2 = base1,
        seq1 = seq2,
        seq2 = seq1,
        pos1 = pos2,
        pos2 = pos1,
        aln
    )])
    
    nuc[, pairPos1 := stri_c(seq1, seq2, pos1, sep = " ")]
    startID <- regionComp[toConvert, stri_c(prot, prot2, inter1[toConvert, 1], sep = " ")]
    endID <- regionComp[toConvert, stri_c(prot, prot2, inter1[toConvert, 2], sep = " ")]
    newStart <- nuc[chmatch(startID, pairPos1), pos2]
    newEnd <- nuc[chmatch(endID, pairPos1), pos2]
    inter1[toConvert, 1] <- newStart
    inter1[toConvert, 2] <- newEnd
    
    # intersection of the 2 intersections = the width of alignment/homology between the 2 regions
    regionComp$interAll <- intersection(inter1[, 1], inter1[, 2], inter2[, 1], inter2[, 2], negative = F)

    # it may be NA if the proteins didn't even align and we had some HSP coordinates as NAs
    regionComp[is.na(interAll), interAll := 0L]
    regionComp[is.na(length), length := 0L]
    regionComp[is.na(pID), pID := 0L]
    regionComp
}


ksMode <- function(ks) {
    # to investigate truncation of TE Ks due to thresholds we imposed, we determine
    # the mode of the kS within  hit group. If this mode is below the upper
    # threshold we set before, then there is little evidence that the hits
    # represent the lower tail of similarities between inherited TEs.

    # We let the hist() function choose the appropriate number of classes
    h <- hist(ks, plot = F)

    # but we modify them because the upper limit of the last class may be higher than the max Ks, so the rightmost class may artificially appear smaller
    newbreaks <- rescale(h$breaks, range(ks))
    h <- hist(ks, breaks = newbreaks, plot = F)

    # the Ks mode
    mode <- rev(h$mids)[which.max(rev(h$counts))]
    list(
        mode = mode,
        nMode = max(h$counts),
        nLast = tail(h$counts, 1),
        class = (max(ks) - min(ks)) / length(h$mids)

        # note that nLast is the size of the class of highest Ks
    )
}

# FUNCTIONS USED TO DRAW FIGURES -------------------------------------------------------

fadeTo <- function(source, dest, amount) {
    # fades source colours into destination colours, by a certain amount (from 0 to 1). 
    # alpha is not managed and any transparency is removed
   
    dest <- rep(dest, length.out = length(source))
    sourceLevels <- col2rgb(source) / 255
    destLevels <- col2rgb(dest) / 255
    delta <- destLevels - sourceLevels
    levels <- t(sourceLevels) + t(delta) * amount
    rgb(levels[, 1], levels[, 2], levels[, 3])
}

saturate <- function(col, amount) {
    # saturates colours by a certain amount from -1 (greyscale) to 1 (max saturation)
    
    amount <- rep(amount, length.out = length(col))
    amount[amount < -1] <- -1
    amount[amount > 1] <- 1
    f <- amount <= 0
    amount[f] <- amount[f] + 1
    hueSat <- rgb2hsv(col2rgb(col))
    diff <- 1 - hueSat[2, ]
    hueSat[2, f] <- hueSat[2, f] * amount[f]
    hueSat[2, !f] <- hueSat[2, !f] + amount[!f] * diff[!f]
    hueSat[hueSat > 1] <- 1
    hsv(hueSat[1, ], hueSat[2, ], hueSat[3, ])
}


linkedBarPlot <- function(height,
                          space = 1,
                          width = 1,
                          col,
                          junCol = fadeTo(col, "white", 0.5),
                          ...) {
    # links sectors between bars of a stacked bar plot
    
    b <- barplot(
        height,
        space = space,
        width = width,
        col = col,
        plot = T,
        ...
    )

    mat <- rbind(0L, height)
    nr <- nrow(mat)
    nc <- ncol(mat)

    x <- as.vector(rbind(b - width / 2, b + width / 2))
    x <- x[-c(1, length(x))]
    x <- matrix(rep(x, nr), ncol = length(x), byrow = T)

    y <- colCumsums(mat)
    
    if (nc > 2) {
        y <- y[, c(1, rep(2:(nc - 1), each = 2), nc)]
    }

    starts <- seq(1, ncol(x), 2)
    ends <- starts + 1L
    x0 <- as.vector(x[, starts])
    x1 <- as.vector(x[, ends])
    y0 <- as.vector(y[, starts])
    y1 <- as.vector(y[, ends])
    
    x <- lapply(2:length(x0), function(i) {
        c(x0[i - 1], x0[i], x1[i - 1], x1[i])
    })
    
    y <- lapply(2:length(x0), function(i) {
        c(y0[i - 1], y0[i], y1[i], y1[i - 1])
    })
    
    if (nc > 2) {
        remove <- seq(nr, length(x), nr)
    } else {
        remove <- length(x) + 1L
    }
    
    m <- Map(polygon, x[-remove], y[-remove], col = junCol, border = NA)
    args <- list(...)
    
    if (hasArg("border")) {
        linksCol <- args[["border"]]
    } else {
        linksCol <- par("col")
    }
    
    if (hasArg("lwd")) {
        linksLWD <- args[["lwd"]]
    } else {
        linksLWD <- par("lwd")
    }
    
    if (hasArg("lty")) {
        linksLTY <- args[["lty"]]
    } else {
        linksLTY <- par("lty")
    }
    
    segments(x0,
        y0,
        x1,
        y1,
        col = linksCol,
        lwd = linksLWD,
        lty = linksLTY
    )
}


xyDistance <- function(x0, y0, x1, y1) {
    # returns the euclidian distance between two points given their coordinates (arguments can be vectors)
    
    dX <- x1 - x0
    dY <- y1 - y0
    l <- sqrt(dX^2 + dY^2)
}

roundedRect <- function(x0,
                        y0,
                        x1,
                        y1,
                        l,
                        rounding = 0.5,
                        radius = NA,
                        plot = T,
                        angle = pi / 2,
                        nv = 100,
                        ...) {
    # draws a rounded rectangle. First 4 arguments are coordinates of two corners.
    # l is the length of the other rectangle side (that not defined by the two
    # corners). rounding (a number from 0 to 1) and radius (the radius of the
    # circle arc replacing a corner) control the amount of rounding. Angle is the
    # angle of the circle arc that replaces a corner and nv the number of vertices
    # in that arc. rounding, radius and angle can be defined separately for each
    # corner
    
    L <- xyDistance(x0, y0, x1, y1)
    s <- ifelse(l < L, l / 2, L / 2)
    
    if (all(!is.na(radius))) {
        rounding <- radius / s
    }
    
    rounding[rounding > 1] <- 1
    rounding[rounding < 0] <- 0
    rounding <- rep_len(rounding, 4)
    angle <- rep_len(angle, 4)
    d <- sign(y0 - y1) * sign(x1 - x0)

    # 1 if descending slope, -1 if ascending or horizontal
    d[d == 0] <- -1
    dX <- abs(x1 - x0)
    dY <- abs(y1 - y0)
    dy <- dX * l / L
    dx <- dY * l / L * d
    x2 <- x1 + dx
    y2 <- y1 + dy
    x3 <- x0 + dx
    y3 <- y0 + dy

    X <- c(x3, x0, x1, x2, x3, x0)
    Y <- c(y3, y0, y1, y2, y3, y0)
    d1 <- xyDistance(X[2:5], Y[2:5], X[1:4], Y[1:4])
    d2 <- xyDistance(X[2:5], Y[2:5], X[3:6], Y[3:6])
    x <- NULL
    y <- NULL
    
    for (corner in 1:4) {
        if (rounding[corner] == 0) {
            x <- c(x, X[corner + 1])
            y <- c(y, Y[corner + 1])
        } else {
            p1 <- s * rounding[corner] / d1[corner]
            p2 <- s * rounding[corner] / d2[corner]
            xa <- X[corner + 1] + (X[corner] - X[corner + 1]) * p1
            ya <- Y[corner + 1] + (Y[corner] - Y[corner + 1]) * p1
            xb <- X[corner + 1] + (X[corner + 2] - X[corner + 1]) * p2
            yb <- Y[corner + 1] + (Y[corner + 2] - Y[corner + 1]) * p2
            midx <- (xa + xb) / 2
            midy <- (ya + yb) / 2
            arc <- arcSegments(xa, ya, xb, yb, angle[corner], nv = nv, draw = F)[[1]]
            
            if (angle[corner] > 0) {
                x <- c(x, arc$x)
                y <- c(y, arc$y)
            }
            else {
                x <- c(x, rev(arc$x))
                y <- c(y, rev(arc$y))
            }
        }
    }
    
    if (plot) {
        polygon(x, y, ...)
    }
    
    invisible(cbind(x, y))
}



linesFromTree <- function(tree, angular = T) {
    # returns a list of lines (x, y coordinates of points) that can be used to draw a tree on a plot

    # the path to each tip
    paths <- nodepath(tree)

    # the X coordinates of nodes/tips on the plot
    X <- node.depth.edgelength(tree)

    # the Y coordinates
    Y <- node.height(tree)

    # used as we scan paths to remove branches already covered by previous paths. We start with the root node
    temp <- getMRCA(tree, tree$tip.label)

    for (i in 1:length(paths)) {

        # position of the last node already covered by previous paths
        m <- max(match(temp, paths[[i]]), na.rm = T)

        # node connexions that are specific to this path (i.e. between this node and the tip)
        nodes <- paths[[i]][m:length(paths[[i]])]
        temp <- union(temp, nodes)
        x <- X[nodes]

        # we get their coordinates
        y <- Y[nodes]
        if (angular) {
            # if this is an angular tree (not V-shaped), we add points for corners
            x <- rep(x, each = 2)
            x <- x[-length(x)]
            y <- rep(y, each = 2)
            y <- y[-1]
        }

        # we use the path list to store our coordinates
        paths[[i]] <- data.frame(x, y)
    }
    setNames(paths, tree$tip.label)
}

outlineTaxon <- function(tree,
                         node,
                         name,
                         x,
                         y,
                         xMargin = 0,
                         yMargin = 0,
                         tipPos,
                         cex = 0.7,
                         ...) {
    # names and outlines a clade in the concave tree. x and y are positions of the
    # clade label, yMargin is the margin between the base of the clade (its age on
    # the tree) and the surrounding rectangle. xMargin is the margin between the X
    # position of the tips (or the label) and the size of the rectangle, tiPos is
    # the Y position of tips, required to draw at the correct position. It is also
    # the position of the top of rectangles surrounding clades

    # age of the clade, which will be position of the basal side of the rectangle surrounding it
    age <- nodeDepth(tree, node)

    # horizontal (X) positions of the tips
    XPos <- node.height(tree)

    # horizontal range of the tips, which will set the sides of the rectangle
    xRange <- range(XPos[tipsForNode(tree, node, itself = T)])
    yM <- 15


    # the bottom of the rectangle may be pushed back down to the rectangle surrounding the label
    bottom <- max(tipPos + age + yMargin, y - (yM - 1))
    circPolygon(
        c(
            xRange[1] - xMargin,
            xRange[2] + xMargin,
            xRange[2] + xMargin,
            xRange[1] - xMargin
        ),
        c(rep(bottom, 2), c(tipPos, tipPos)),
        ...

        # draws the rectangle around the clade
    )


    # required for the rounded rectangle around the label because the different scales between X and Y 
    # coordinates on the circle messes with the rounding, also required to determine label width
    ratio <- 1 / xRatioForY(y)

    # label width in X coordinates, required to determine the size of the rectangle surrounding the label
    wordWidth <- strwidth(name) * cex / ratio

    leftX <- x - xMargin
    rightX <- x + wordWidth + xMargin
    if (abs(leftX - xRange[1] + xMargin) < 2) {
        leftX <- xRange[1] - xMargin

        # adjusts the left side of the rectangle to the left side the rectangle 
        #surrounding the clade, if it's close enough, for best visual effect
    }
    if (abs(rightX - xRange[2] - xMargin) < 2) {
        rightX <- xRange[2] + xMargin
        # same for the right side (considering that it's very unlikely that the left side was close enough in this case)
    }
    
    # we compute xy coordinates of rounded rectangle around the label
    rect <- roundedRect(
        leftX * ratio,
        y - yM,
        rightX * ratio,
        y - yM,
        strheight(name) * cex + 2 * yM,
        rounding = c(0, 0, 1, 1),
        nv = 20,
        plot = F
    )

    # we draw it after reestablishing proper X dimensions
    circPolygon(rect[, 1] / ratio, rect[, 2], ...)
    
     # draws clade label
    circText(
        x,
        y,
        labels = name,
        curved = T,
        col = grey(0.3),
        cex = cex,
        adj = 0
    )
    invisible(NULL)
}
