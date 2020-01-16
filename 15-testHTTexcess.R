

## %######################################################%##
#                                                          #
####       This stage permutates species to test        ####
####       whether some clades have more htt than       ####
####          expected by chance. This can be           ####
####          run at any point after stage 12.          ####
#                                                          #
## %######################################################%##

# The principle is to replace a species by a random one among the 307 species
# and compute the number of transfers that involves each species (or rather, each taxon)
# Hence, the distribution of htt event per species is not changed in a sense that 
# if species 1 was involved in 10 events and species 2 had 20, swapping these species
# will not affect the numbers. There would sill be 10 transfers for one species and
# 20 for another.

# as several species are found for each clade of a hit group, we only
# consider the pair of species involved in the best hit (see paper)

source("HTvFunctions.R")
tree <- read.tree("additional_files/timetree.nwk")

# we import data file provided with the paper, which is a table of hits representing HTT
retainedHits <- fread("supplementary-data4-retained_hits")



# STEP ONE, permutation of species --------------------------------------------------------------------------------
# we launch the dedicated script to perform the permutations on 20 CPUs
system('sbatch --mail-type=BEGIN,END,FAIL --cpus-per-task=20 --mem=20G --wrap="Rscript permutateSpeciesInHits.R 20"')



# STEP 2, contrasting the number of HTT per taxa in real and randomised hits.---------------------------------------

# we import the permutations generated above (= pairs of species representing real and simulated HTT)
randomHTTs <- readRDS("permutations/allPermutations_Node.308.RDS")

taxa <- fread("additional_files/namedClades.txt")

# we retrieve the species (tip numbers in the tree) composing each taxon that will be shown on the figure
sps <- tipsForNodes(tree, taxa[onPlot == T, node])

# taxon[x] below gives the taxon of species x
taxon <- sps[, node[order(tip)]]
    

# this function counts htt per taxon in a given TE super family (or batch within a larger superfamily) 
# to obtain statistics (means, quantiles) for number of simulated transfers for each taxon of interest
countHTTs <- function(superfamily) {

    # we "unfold" the species pairs in htts, separating species 1 and 2, duplicating HTTs. 
    # Effectively, if a hit involves 2 species of the same taxon, two HTTs will be counted for the taxon
    unFolded <- randomHTTs[[superfamily]][, .(sp = c(sp1, sp2), repl = c(repl, repl))]
   
    # we add a column denoting the taxon of each species
    unFolded[, taxon := taxon[sp]]

    # we count the numbers of transfers involving each taxon in each permutation (note that repl == 0 for real HTT events)
    counts <- unFolded[, .(N = .N / 2), by = .(taxon, repl)]
    
    # we retrieve the number of permutations
    nRepl <- max(randomHTTs[[superfamily]]$repl)
    
    # we need to add 0 counts for taxa not involved in HTT in each replicate
    # since they are simply missing (if we didn't do this, the stats would be biased upwards)
    # for this we need to know the taxa involved
    uTaxa <- unique(taxon)
    
    # ad we create a table indicating 0 transfer for each taxon in each replicate
    missingTaxa <- data.table(
        taxon = rep(uTaxa, nRepl),
        repl = rep(1:nRepl, each = length(uTaxa)),
        N = 0
    )
    
    # we place these rows at the end of our previous table
    counts <- rbind(counts, missingTaxa)

    # and we remove rows of 0 counts for taxa that were already present
    counts <- counts[!duplicated(data.table(taxon, repl))]
    
    # we now obtain the statistics we need, per taxon
    stats <- counts[repl > 0L, .(
        simulated = mean(N),       # mean number of transfer over permutations
        minNumber = min(N),        # min number '''
        maxNumber = max(N),        # max number '''
        qLc = quantile(N, 0.005),  # and the various quantiles to test for significance
        qRc = quantile(N, 0.995),
        qL5 = quantile(N, 0.025),
        qR5 = quantile(N, 0.975),
        superfamily
    ), by = taxon]

    # and we add number of real number of HTT events to these stats
    merge(stats, counts[repl == 0L, .(taxon, observed = N)], by = "taxon", all = T)
}


# we apply the function to all super families (batches)
permutationStats <- lapply(names(randomHTTs), countHTTs)

# we stack results in a single table
permutationStats <- rbindlist(permutationStats)

# NAs must be replaced with zeros for taxa not observed in the real transfers
permutationStats[is.na(observed), observed := 0]


# We make one plot per TE class ---------------------------------------------

# we thus add a column for TE classes
permutationStats[, class := ifelse(
    test = grepl("CMC|hAT|Mariner|Maverick|Merlin|PIF|PiggyBac|DNA",
                 superfamily),
    yes = "DNA",
    no = "RNA"
)] 

# we compute stats for randomised and observed HTT numbers of the different taxa we selected and for each TE class
# this means that we sum stats obtained over super families within classes
statsPerClass <- permutationStats[, .(
    simulated = sum(simulated),
    minNumber = sum(minNumber),
    maxNumber = sum(maxNumber),
    qLc = sum(qLc),
    qRc = sum(qRc),
    qL5 = sum(qL5),
    qR5 = sum(qR5),
    observed = sum(observed)
), by = .(class, taxon)]


# we count the number of simulated HTT per taxon, to put the taxa the most involved on top
sumSimulated <- statsPerClass[, sum(simulated), by = .(node = taxon)]
statsPerClass <- statsPerClass[order(sumSimulated[match(taxon, node), -V1], class), ]

# building the connected barplots (figure 3)

# we attribute taxa-specific colour used for the plots, the same as on figure 2
statsPerClass[, col := taxa[match(taxon, node), col]]
    
# this function generates the linked barplot figure. 
# Almost all the code adds stars for siginificance levels
linkedPlots <- function(dt, space = 2, legend = F, ...) {
    
    b <- dt[, linkedBarPlot(
        cbind(simulated, observed),
        col = col,
        junCol = fadeTo(col, "black", 0.4),
        space = space,
        main = ifelse(class[1] == "DNA", "DNA transposons", "Retrotransposons"),
        border = grey(0.3),
        ...
    )]

    # no star by default
    dt[, stars := ""]
    
    # difference at the 5% level when observed numbers are beyond the right quantiles
    dt[observed < qL5 | observed > qR5, stars := "*"]  
    
    # at the 1% level
    dt[observed < qLc | observed > qRc, stars := "**"]
    
    # at the 1/1000 level
    dt[observed < minNumber | observed > maxNumber, stars := "***"]

    # TRUE when there are significant differences
    significant <- dt$stars != ""
    
    # we will add stars in this case
    labels <- dt[significant, ]

    # below, we compute the vertical location of the stars
    y <- cumsum(c(0, dt$observed))
   
    # we need to know the mid vertical position of barplot sectors to which we want to add stars
    mids <- (y[which(significant)] + y[which(significant) + 1L]) / 2
    
    # and we adjust the vertical position of future stars to avoid collisions.
    yPos <- c(0, mids)
    for (i in 2:length(yPos)) {
        delta <- yPos[i] - yPos[i - 1L]
        if (delta < 15) {
            yPos[i] <- yPos[i] + 15 - delta
        }
    }
    
    yPos <- c(yPos[-1], sum(dt$observed))
    for (i in length(yPos):2) {
        delta <- yPos[i] - yPos[i - 1L]
        if (delta < 15) {
            yPos[i - 1L] <- yPos[i - 1L] - 15 + delta
        }
    }
    
    yPos <- yPos[-length(yPos)]
    t <- labels[, text(2 * space + 2.55, yPos, stri_c(round(observed / simulated, 2), stars), adj = c(0, 0.5))]
    
    # we add segments to link the stars to barplot sectors
    segments(2 * space + 2,
        mids,
        2 * space + 2.5,
        yPos,
        col = dt$col[significant],
        lwd = 1
    )
    
    
    # we add a legend if required (as there will be 2 plots on the figure)
    if (legend) {
        legend(
            x = 2 * space + 3.5,
            y = sum(dt$observed),
            legend = rev(taxa[match(dt$taxon, node), name]),
            bty = "n",
            fill = rev(dt$col),
            border = grey(0.3)
        )
    }
}



pdf(
    file = "Figure3.pdf",
    width = 8,
    height = 4.5
)


# in this layout, the rightmost section is there to leave room for the legend 
layout(cbind(1, 2, 3), widths = c(1, 1, 0.4))

par(
    lwd = 0.5,
    xpd = NA,
    las = 1,
    mai = c(0.6, 0.55, 0.55, 0.6)
)

# the barplots for DNA transposons
linkedPlots(statsPerClass[class == "DNA"], ylab = "number of transfers")

#and for retrotranspons
linkedPlots(statsPerClass[class == "RNA"], legend = T)

dev.off()



# for ray-finned fishes, we record results for each super family (table S1) ------------------------------------------
rayFin <- permutationStats[taxon == 312L, .(
    mean_over_simulations = round(sum(simulated), 2),
    max_over_simulations = sum(maxNumber),
    observed = sum(observed)
), by = testSplit(superfamily, ".", 1)]
writeT(rayFin, "tableS1.txt")





# ADDITIONAL ANALYSES NOT DETAILED IN THE PAPER ---------------------------------------------
# investigating if more transfers between "fishes" and tetrapods involved aquatic tetrapods 
# we first create a list indicating the main habitat of the species (species are tip numbers on the timetree)
tetrapods <- tipsForNode(tree, 375)

# we retrieve the aquatic species among tetrapods
aquatic <- c(
    tipsForNodes(tree, c(376, 472, 474, 514, 501, 463))$tip,
    grep("Trichechus", tree$tip.label)
) 

# fishes as a paraphyletic group (non-tretrapods)
fishes <- setdiff(1:length(tree$tip.label), tetrapods)

habitatEffect <- function(tetrapods, aquatic) {
    # this permutes tetrapod species among the tetrapod-fish HTT. This is fast
    # since this cannot general illegal HTTs (the randomised hit will still be
    # between fishes and tetrapods)
    
    tetraHits <- HTThits[(sp1 %in% tetrapods &
        sp2 %in% fishes) |
        (sp2 %in% tetrapods &

            # hits between a tetrapod and a non-tetrapod
            sp1 %in% fishes)]

    # always puts the tetrapod species on the left
    tetraHits[!sp1 %in% tetrapods, c("sp1", "sp2") := .(sp2, sp1)]

    # renumbers the tetrapod species to help the simulation
    tetraHits[, sp1 := match(sp1, tetrapods)]

    aquatic <- match(aquatic, tetrapods)
    l <- length(tetrapods)
    permu <- replicate(10000, sample(l))
    sp1 <- tetraHits$sp1 + rep(1:10000 - 1L, each = nrow(tetraHits)) * l
    sp1 <- permu[sp1]
    pairs <- rbind(
        data.table(sp1 = tetraHits$sp1, repl = 0L),
        data.table(sp1, repl = rep(1:10000, each = nrow(tetraHits)))
    )
    pairs[, habitat := "land"]
    pairs[sp1 %in% aquatic, habitat := "water"]
    perPerm <- pairs[repl > 0L, .N, by = .(repl, habitat)]
    simu <- perPerm[, .(
        mean = mean(N),
        q1 = quantile(N, 1 / 1000),
        q999 = quantile(N, 999 / 1000)
    ), by = habitat]
    real <- pairs[repl == 0L, .(observed = .N), by = habitat]
    merge(simu, real, by = "habitat")
}

habitatEffect(tetrapods, aquatic)

# same, but only within birds and mammals
# to increase power, we use all hit groups, not just the "independent" ones
HTThits <- retainedHits[!duplicated(hitgroup), data.table(
    sp1 = chmatch(sp1, tree$tip.label),
    sp2 = chmatch(sp2, tree$tip.label)
)]

birdsMams <- tipsForNodes(tree, c(477, 391))$tip
aquaticBM <- intersect(aquatic, birdsMams)
habitatEffect(birdsMams, aquaticBM)

# same, but ignoring birds and mammals
nonBirdsMams <- setdiff(tetrapods, birdsMams)
aquaticNBM <- setdiff(aquatic, birdsMams)
habitatEffect(nonBirdsMams, aquaticNBM)
