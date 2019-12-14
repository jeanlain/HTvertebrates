## %######################################################%##
#                                                          #
####           This stages makes the convave            ####
####        tree showing HTT events (figure 2).         ####
#                                                          #
## %######################################################%##

# This can be run at any point after step 12

source("HTvFunctions.R")


# we draw the HTTs as arcs connecting tips on a concave tree.
# we use different colors for DNA and RNA TEs.

tree <- read.tree("timetree.nwk")

# data file provided with the paper, which is a table of hits representing HTT
retainedHits <- fread("supplementary-data4-retained_hits.txt")

# we get the data ready for the plot-----------------------------

# we create a table of the best hit per transfer, the one we will show on the tree
# we thus pace hits whith highest pID on top
setorder(retainedHits, -pID)

# we extract these hits, and the columns we need
connections <- retainedHits[!duplicated(hitgroup), .(sp1, sp2, superfamily, hitgroup, independent)]

# we will show different colors for the TE class, we determine them based on superfamily name
connections[, class := ifelse(test = grepl("CMC|hAT|Mariner|Maverick|Merlin|PIF|PiggyBac", superfamily),
    yes = "DNA",
    no = "RNA"
)]


connections <- connections[
    sample(.N), # We shuffle rows to ensure that arcs will be drawn in random order,
    # to avoid visual bias between TE classes
    data.table(independent, # We retain the "independent" column to leave us the
        # option to show independent HTTs or all hit groups.
        # Species names are converted to integer that
        # correspond to tip numbers on the tree.
        tip1 = match(sp1, tree$tip.label),
        tip2 = match(sp2, tree$tip.label),
        col = ifelse(test = class == "DNA",
            yes = rgb(0.9, 0.4, 0.4, 0.4), # red
            no = rgb(0.4, 0.4, 0.9, 0.4)
        )
    )
] # blue


# rotates some nodes to avoid unnecessary long horizontal branches
for (node in c(308, 389, 390, 477, 391)) tree <- rotate(tree, node)

# this will be the position of species along the circle
xTipPos <- node.height(tree)

# we obtain the max depth of the tree (age of mrca)
treeDepth <- max(node.depth.edgelength(tree))

# see function. This uses ape to return lines (vertex coordinates) to draw curved branches
brancheLines <- linesFromTree(tree)

# Y position of the tree tips = distance from the center of the plot area
tipPos <- treeDepth * 1.7 + 10

# we locate clades younger than 120 My within which HTT was not inferred
# (excluding clades composed of just one species, as we don't need to outline those)
outlinedClades <- cladesOfAge(
    tree = tree,
    age = 120,
    withTips = T,
    names = F,
    singletons = F
)

# we import a table of taxa (clades) with latin and English names corresponding to tree nodes
# and colors used in figures 2 and 3, and offset to draw their name at the right place on the tree.
# the offset was somewhat determined by trial an error.
# Some taxa are not shown on the figures but are still in this table
taxa <- fread("namedClades.txt")

# we compute the age of the mrca of the taxa, which is use to
# find the right position to write taxa names on the tree
taxa[, age := nodeDepth(tree)[node]]

# for taxa shown on the tree but not drawn on figure 3 (the connected barplots), we attribute colors
taxa[
    onTree == T & onPlot == F,
    col := rep_len(saturate(brewer.pal(12, "Set3"), -0.5), .N)
]

# we loaf functions required to make the radial figure
source("circularPlots.R")

pdf("figure2.pdf", 7, 7) # figure 2 of the paper

# see circularPlots.R for these functions
initializePlot(tipPos + treeDepth, -1, max(xTipPos) + 2L, 180)

# we draw nice curved rectangles around some clades and name the clade --------------------
# Calling the function within the taxa table does not produce a correct pdf, hence the us of with()
for (i in which(taxa$onTree)) {
    with(taxa[i, ], outlineTaxon(tree, node, name,
        x = xTipPos[node] + xOffset,
        y = tipPos + age + 20 + yOffset,
        xMargin = 0.5, yMargin = 5,
        col = col,
        tipPos = tipPos - 3,
        cex = 0.6, border = NA
    ))
}

# we draw the timescale ----------------------------------------------
# the ages in My
ages <- seq(100, 600, 100)

# the tick-marks of the ages
circSegments(
    x0 = -2.3, # the position is just before the first tip (hence degative)
    y0 = tipPos + ages, # the "depth" of the tick is defined by the ages
    x1 = -1.7, # to create segments of appropriate length
    curved = F, # these segments need not be curved according to the tree
    lend = 1, lwd = 0.5
)

# we drow the ages themselves
circText(
    x = rep(-0.8, length(ages)),
    y = tipPos + ages,
    labels = ages / 10,
    cex = 0.5,
    correct = F
)

# we draw the concave tree -------------------------------------------------
l <- lapply(
    X = brancheLines,
    FUN = function(v) {
        circLines(
            x = v$y, # x takes column y because fundamentally, this tree is drawn horizontallt
            y = tipPos + treeDepth - v$x,
            lwd = 1, lend = 2, col = grey(0.3)
        )
    }
)

# we add colored points at the node of outlined clades ------------------------
taxa[
    onTree == T & age > 1, # age > 1 to avoid drowing points for single-tip taxa
    circPoints(
        x = xTipPos[node],
        y = tipPos + age,
        pch = 21,
        col = grey(0.3),
        bg = col,
        cex = 0.5
    )
]

# we add curved grey segments above tipes of clades younger than 120 My ---------
outlinedClades[, circSegments(
    x0 = min(xTipPos[tip]) - 0.2,
    y0 = tipPos - 20,
    max(xTipPos[tip]) + 0.2,
    col = grey(0.5),
    lwd = 3,
    lend = 1
),
by = node
]

# we show HTT event with arcs. ---------------------------------------------------
# Note that we only show "independent" transfers but we could draw
# all hit groups by removing the filter on the rows if we wanted

co <- connections[independent == T, connectionArcs(
    x0 = xTipPos[tip1],
    y0 = tipPos - 35,
    x1 = xTipPos[tip2],
    col = col,
    lwd = 1
)]


dev.off()
