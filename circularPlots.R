## %######################################################%##
#                                                          #
#### these functions make circular plots, which may be  ####
####            useful to represent cicular             ####
####         genomes among other applications.          ####
#                                                          #
## %######################################################%##


# The approach is to curve the X axis and leave the Y axis perpendicular to the X
# axis. The Y axis is thus oriented toward the centre of the circle, such that
# any point at Y = 0 is at the center of the circle (Y should never be < 0). The
# Y scale is that of the plot area, but the X scale is changed and defined in a
# parameter vector (circ__, which the user is not supposed to touch but
# unknowingly defines in initializePlot()) that characterises the plot. circ__
# specifies the radius of the circle (in Y coordinates) on a normal plot,
# together with the x coordinate at noon (i.e. starting x, as in a clock), and
# the x coordinate at midnight, the angle of the circle radius at x = 0
# (noon/midnight) the starting x (as in a trigonometric circle, in degrees), the
# direction of rotation as X increases (1 is counterclockwise)


# we load the required functions, but this script was non written specifically for the analysis of HTT
source("HTvFunctions.R")

last <- function(x) {
  x[length(x)]
}

initializePlot <- function(radius,
                           xAtNoon,
                           xAtMidnigth,
                           angleAtX0 = 90,
                           dir = -1,
                           xlab = "",
                           axes = F,
                           ylab = "",
                           asp = 1,
                           bty = "n",
                           mai = rep(0, 4),
                           ...) {
  # this function initializes the plot area according to the parameters
  # described above.
  # by default, there is no margin, no axes, and no box on the plot


  circ__ <<- as.numeric(c(radius, xAtNoon, xAtMidnigth, angleAtX0, dir))
  names(circ__) <<- c("radius", "xAtNoon", "xAtMidnigth", "angleAtX0", "dir")

  par(mai = mai)

  plot(
    c(-radius, radius),
    c(-radius, radius),
    type = "n",
    xlab = xlab,
    axes = axes,
    ylab = ylab,
    asp = asp,
    bty = bty,
    ...
  )
}


angleOfPoint <- function(x) {
  # determines the angle of a point to draw, as in a trigonometric circle,
  # the origin being the center of the plot area.
  # The user is not supposed to use this function
  degreePerXUnit <- 360 / (circ__[3] - circ__[2])
  circ__[4] + (x - circ__[2]) * circ__[5] * degreePerXUnit
}


convertCoords <- function(x, y) {
  # Internal function that converts coordinates between
  # a regular plot and its circularized version.

  if (!exists("circ__")) {
    stop("initializePlot() not called yet")
  }
  x[is.na(x)] <- x[1]
  y[is.na(y)] <- y[1]
  if (length(x) < length(y)) {
    x <- rep(x, length.out = length(y))
  } else {
    y <- rep(y, length.out = length(x))
  }
  angle <- angleOfPoint(x) * 2 * pi / 360
  data.frame(x = y * cos(angle), y = y * sin(angle))
}



convertAngle <- function(x, angle) {
  # Internal function

  angleOfPoint(x) + angle - 90
}


circPoints <- function(x, y, ...) {
  # draws points on a circular plot (equivalent of the points() function).
  # Point shapes are not rotated

  coords <- convertCoords(x, y)
  points(coords$x, coords$y, ...)
}



circSegments <- function(x0,
                         y0,
                         x1 = x0,
                         y1 = y0,
                         curved = T,
                         ...) {
  # draws segments on a circular plot, similar to segments()
  # Segments will be curved unless those for which "curved" == F

  coords <- data.table(x0, y0, x1, y1)

  # we compute the number of times to subdivide segments so they appear properly curved
  div <- coords[, subdivisions(x1 - x0, (y0 + y1) / 2)]
  f <- div > 1 & curved

  for (i in which(f)) {
    # those that need to be subdivided are drawn as curved lines
    c <- coords[i, circLines(c(x0, x1), c(y0, y1), ...)]
  }

  # the others as straight segments
  coords0 <- coords[!f, convertCoords(x0, y0)]
  coords1 <- coords[!f, convertCoords(x1, y1)]

  segments(coords0$x, coords0$y, coords1$x, coords1$y, ...)
}

circLines <- function(x,
                      y,
                      curved = T,
                      draw = T,
                      ...) {
  # draws lines on a circular plot, similar to the lines() function
  # if draw is FALSE, only returns point coordinates of the lines

  coords <- subdivide(x, y, do = curved)
  coords <- convertCoords(coords$x, coords$y)

  if (draw) {
    lines(coords$x, coords$y, ...)
  } else {
    invisible(coords)
  }
}



circPolygon <- function(x,
                        y,
                        closeAtY = "n",
                        connect = F,
                        curved = T,
                        ...) {
  # draws a polygon on a circular plot, somilar to polygon()

  if (connect) {
    # if we want the polygon to be closed, that is, if the 2 ends join (likely at noon)
    x <- c(circ[2], x, circ[3])
    y <- c(mean(y[1], last(y)), y, mean(y[1], last(y)))
  }

  if (is.numeric(closeAtY)) {
    # used to automatically specify that the polygon should have
    # a bottom side at closeAtY  (so we don't have to close it manually)

    x <- c(x, last(x), x[1])
    y <- c(y, closeAtY, closeAtY)
    if (length(curved) > 1) {
      curved <- c(curved, T)
    }
  }

  coords <- subdivide(x, y, do = curved)
  coords <- convertCoords(coords$x, coords$y)
  polygon(coords$x, coords$y, ...)
}

xRatioForY <- function(y) {
  # internal function
  abs(circ__[3] - circ__[2]) / (y * 2 * pi)
}




circText <- function(x,
                     y,
                     labels,
                     rotate = T,
                     adj = par("adj"),
                     srt = par("srt"),
                     cex = par("cex"),
                     col = par("col"),
                     curved = F,
                     correct = rotate,
                     ...) {
  # displays text on a circular plot. Words can be curved to follow the circular X
  # axis. rotate indicates whether the labels should be rotated to keep the same
  # angle (specified by srt, in degrees) relative to the X axis. correct makes
  # sure that no text is written upside down. Non-default fonts are not managed
  # (never tested)

  m <- max(length(x), length(y), length(labels), length(col))
  col <- rep(col, length.out = m)
  x <- rep(x, length.out = m)
  y <- rep(y, length.out = m)
  labels <- rep(labels, length.out = m)

  coords <- convertCoords(x, y)

  if (!rotate) {
    # if we don't rotate the text, we can already draw it
    text(
      coords$x,
      coords$y,
      labels,
      srt = srt,
      adj = adj,
      col = col,
      cex = cex,
      ...
    )
  }

  if (length(adj) == 1) {
    adj <- rep(adj, 2)
  }

  adjx <- adj[1]
  adjy <- ifelse(!is.na(adj[2]), adj[2], 1 / 2)
  toReverse <- rep(F, m)

  if (correct) {
    wordBottoms <- y - strheight(labels) * cex * adjy
    wordWidths <- strwidth(labels) * cex * xRatioForY(wordBottoms)

    wordMids <- x + wordWidths * (0.5 - adjx)
    temp <- convertCoords(wordMids, wordBottoms)
    toReverse <- temp$y < 0
    adjy <- rep(adjy, m)
  }

  if (!curved & rotate) {
    adjx <- rep(adjx, m)
    adjx[toReverse] <- 1 - adjx[toReverse]
    adjy[toReverse] <- 1 - adjy[toReverse]
    angle <- convertAngle(x, srt)
    angle[toReverse] <- angle[toReverse] + 180
    adju <- Map("c", adjx, adjy)
    
    for (i in 1:length(x)) {
      text(
        coords$x[i],
        coords$y[i],
        labels = labels[i],
        srt = angle[i],
        col = col[i],
        adj = adju[[i]],
        cex = cex,
        ...
      )
    }
  }

  if (srt == 0 & curved) {
    adjy[toReverse] <- adjy[toReverse] - 1
    labels[toReverse] <- stri_reverse(labels[toReverse])
    wordBottoms <- y - strheight(labels) * cex * adjy
    wordWidths <- strwidth(labels) * cex * xRatioForY(wordBottoms)
    lets <- strsplit(labels, split = "")

    letterMids <- Map(function(lets, y) {
      mids(cumsum(
        c(0, strwidth(lets) * cex * xRatioForY(y) * -circ__[5])
      ))
    }, lets, wordBottoms)

    leftWordPos <- x - wordWidths * adjx * -circ__[5]
    letterMids <- Map("+", letterMids, leftWordPos)
    angle <- convertAngle(unlist(letterMids), srt)
    angle[rep(toReverse, nchar(labels))] <- angle[rep(toReverse, nchar(labels))] + 180
    coords <- convertCoords(unlist(letterMids), rep(wordBottoms, nchar(labels)))
    col <- rep(col, nchar(labels))
    labels <- unlist(lets)

    for (i in 1:nrow(coords)) {
      text(
        coords$x[i],
        coords$y[i],
        labels = labels[i],
        srt = angle[i],
        col = col[i],
        adj = c(0.5, 0),
        cex = cex,
        ...
      )
    }
  }
}




subdivisions <- function(diffX, y, resolution = 0.04) {
  # returns the number of subdivisions needed within segments so they appear
  # properly curved on a circular plot. Resolution is the distance between two
  # points in inches. Subdivisions only concern the X coordinates since the
  # circular plot is only curved on the x axis.

  # we compute length (in inches) of a circle of radius y units,
  # considering that circ__[1] is the max radius that takes the
  # whole plotting region (almost, but we can approximate)
  meanPerimetre <- min(par("pin")) * pi * y / circ__[1]

  inchesPerXUnit <- meanPerimetre / abs(circ__[3] - circ__[2])

  div <- as.integer(abs(diffX) * inchesPerXUnit / resolution) + 1L
}

subdivide <- function(x, y, do = T) {
  # subdivides a segmented line so that it forms a curve when displayed on a circular plot

  pos <- 2:length(x)
  diffX <- x[pos] - x[pos - 1]
  div <- subdivisions(diffX, (y[pos] + y[pos - 1]) / 2)
  div[!do] <- 1L

  subd <- function(div, x, diffX) {
    x + (1:div - 1) / div * diffX
  }

  resX <- c(unlist(Map(subd, div, x[2:length(x) - 1], diffX)), last(x))
  resY <- c(unlist(Map(subd, div, y[2:length(y) - 1], y[pos] - y[pos - 1])), last(y))
  data.frame(x = resX, y = resY)
}


connectionArcs <- function(x0,
                           y0,
                           x1,
                           y1 = y0,
                           draw = T,
                           ...) {
  # connects two points of a circular plot with an circle arc that intersects the
  # circle at perpendicular angles. x0 and y0 refer to the first point. Arguments
  # can be vectors.

  # we compute the portion of the circle that the two points span on the x axis
  portion <- abs(x1 - x0) / abs(circ__[3] - circ__[2])

  # we ensure that arcs are always drawn within the circle,
  # given how the arcSegments function works (negative arc angles may be specified)
  f <- sign((x1 - x0) * (circ__[3] - circ__[2]) * circ__[5] * (0.5 - portion))

  # when points are opposite (portion = 0.5), the arc angle is 0 (straight line),
  # when they are at the same x position (portion = 0), the arc angle is pi (or -pi).
  # This ensures that the arc always crosses the plot circle perpendicularly.
  # Note that we only consider the smaller of the two sections between the points.
  angle <- 2 * pi * (0.5 - pmin(portion, 1 - portion)) * -f
  starts <- convertCoords(x0, y0)
  ends <- convertCoords(x1, y1)

  seg <- arcSegments(
    starts$x,
    starts$y,
    ends$x,
    ends$y,
    angle = angle,
    draw = draw,
    ...
  )

  invisible(seg)
}
