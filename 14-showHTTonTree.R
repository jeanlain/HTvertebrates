#### makes the convave tree showing HTT events (figure 2). This can be run at any point after step 12

source('HTvFunctions.R')
tree = read.tree("timetree.nwk")
retainedHits = fread("gunzip -c data4-retained_hits.txt.gz")		#data file provided with the paper, which is a table of hits representing HTT

setorder(retainedHits, -pID)																		#places hits whith highest pID on top
connections = retainedHits[!duplicated(hitgroup), .(sp1,sp2, superfamily, hitgroup, independent)]	#so we can create a table of the best hit per transfer, the one we will show on the tree
connections[,class :=  ifelse(grepl("CMC|hAT|Mariner|Maverick|Merlin|PIF|PiggyBac", superfamily),"DNA","RNA")]		#attributes TE class to superfamilies in retained hits
connections = connections[sample(.N),data.table(independent, tip1 = match(sp1, tree$tip.label), tip2 = match(sp2, tree$tip.label), col = ifelse(class =="DNA", rgb(0.9,0.4,0.4,0.4), rgb(0.4,0.4,0.9,0.4)))]				#this table will be used to draw the HTT as circle arcs. Note the different colors for DNA and RNA TEs. Species names are converted to integer that correspond to tip numbers on the tree. We retain the "independent" column to leave us the option to show independent HTTs or all hit groups. sample(.N) shuffles rows to ensure that arcs will be drawn in random order, to avoid visual bias between TE classes
 
for(node in c(308, 389, 390, 477, 391)) tree = rotate(tree, node) 	#rotates some nodes to avoid unnecessary long horizontal branches
xTipPos = node.height(tree)								#this will be the position of species along the circle
treeDepth = max(node.depth.edgelength(tree))
brancheLines = linesFromTree(tree)						#see function. This returns lines (vertex coordinates) to draw curved branches
tipPos = treeDepth*1.7+10								#Y position of the tips = distance from the center of the plot 
collapsedClades = cladesOfAge(tree, 120, withTips = T, names = F, singletons = F)		#clades younger than 120 My within which HTT was not inferred (excluding clades composed of just one species, as we don't need to outline those)

taxa = fread("namedClades.txt")		#a table of taxa (clades) with latin and English names corresponding to tree nodes and colors used in figures 2 and 3, and offset to draw their name at the right place on the tree. Some taxa are not shown on the figures but are still in this table
taxa[,age := nodeDepth(tree)[node]]
taxa[onTree ==T & onPlot == F, col :=rep_len(saturate(brewer.pal(12, "Set3"),-0.5),.N)]		#for taxa shown on the tree but not drawn on figure 3 (the connected barplots), we attribute colors

source('circularPlots.R')				#functions required to make the radial figure
pdf("figure2.pdf", 7,7)					#figure 2 of the paper
initializePlot(tipPos + treeDepth, -1, max(xTipPos)+2L, 180)
for(i in which(taxa$onTree)) with(taxa[i,], outlineTaxon(tree, node, name, x=xTipPos[node]+xOffset, y=tipPos+age+20+yOffset, xMargin=0.5, yMargin=5, col=col, tipPos=tipPos-3, cex=0.6, border=NA))  #this draws nice curved rectangles around some clades. Calling the function within the data table does not produce a correct pdf, for some reason. 

times  = seq(100, 600, 100)														#ages shown on the timescale
circSegments(-2.3, tipPos + times, -1.7, curved = F, lend = 1, lwd = 0.5)		#the tick-marks of the ages
circText(rep(-0.8, length(times)), tipPos + times, labels = times/10, cex = 0.5, correct = F)	#writes the ages themselves

l =lapply(brancheLines, function(v) circLines(v$y, tipPos + treeDepth-v$x, lwd = 1, lend = 2, col = grey(0.3)))			#draws the concave tree
taxa[onTree == T & age > 1, circPoints(xTipPos[node], tipPos + age, pch = 21, col = grey(0.3), bg = col, cex =0.5)]		#adds colored points at the node of outlined clades

collapsedClades[, circSegments(min(xTipPos[tip])-0.2, tipPos-20, max(xTipPos[tip])+0.2, col = grey(0.5), lwd = 3, lend = 1), by = node]	#adds curved grey segments above tipes of clades younger than 120 My
co = connections[independent == T,connectionArcs(xTipPos[tip1], tipPos-35, xTipPos[tip2], col = col, lwd = 1)]		#shows HTT event with arcs. Note that we only show "independent" transfers but we could draw all hit groups by removing the filter on the rows

dev.off()