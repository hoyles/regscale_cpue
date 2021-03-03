# load aggregated data
#lldat <- read.csv("IOTC-2016-DATASETS-CELongline.csv",stringsAsFactors = F)
prepCE <- function(CE) {
  CE[1,]
  table(CE$MonthStart)
  # convert these grids to centlat and centlong
  CE$Grid <- as.character(CE$Grid)
  table(CE$Grid)
  CE <- CE[CE$Grid != "F57",]
  CE <- CE[CE$Grid != "9000080",]
  CE <- CE[CE$Grid != "9000020",]
  addlat <- c(5,10,10,20,1,5)
  addlon <- c(10,20,10,20,1,5)
  addtp <- as.numeric(substring(CE$Grid,1,1))
  a <- as.numeric(substring(CE$Grid,2,2))
  CE$Grid[a == 0]
  sg <- c(1,-1)[as.numeric(substring(CE$Grid,2,2))]

  CE$latin <- as.numeric(substring(CE$Grid,3,4)) * sg
  CE$latout <- CE$latin + (addlat[addtp] * sg)
  CE$lonmin <- as.numeric(substring(CE$Grid,5,7))
  CE$lonmax <- CE$lonmin + addlon[addtp]
  CE$centlat <- (CE$latin + CE$latout)/2
  CE$centlon <- (CE$lonmin + CE$lonmax)/2
  CE$lat5 <- 5 * floor(CE$centlat/5) + 2.5
  CE$lon5 <- 5 * floor(CE$centlon/5) + 2.5
  
  colnames(CE)
  head(CE)
  table(CE$MonthS,CE$MonthE)
  # change the approach from an array of lengths across to a single column
  CE$yrqtr <-  as.factor(CE$Year + rep(c(0.125,0.375,0.625,0.875),each = 3)[CE$MonthE])
  CE$yq <- as.numeric(as.character(CE$yrqtr))
  CE$latlong <- as.factor(paste(CE$lat5,CE$lon5, sep = "_"))
  CE$latlon1 <- as.factor(paste(sg*0.5 + CE$latin,0.5+CE$lonmin, sep = "_"))
  CE$latlon2 <- as.factor(paste(sg*1   + 2*floor(CE$latin/2),1+2*floor(CE$lonmin/2), sep = "_"))
  CE$latlon3 <- as.factor(paste(sg*1.5 + 3*floor(CE$latin/3),1.5+3*floor(CE$lonmin/3), sep = "_"))
  
  CE <- setup_IO_regions(CE, regY2=TRUE, regA4 = TRUE)
  return(CE)
}

plot_patterns <- function(llmn, sp, spreg) {
  lab1 <- unlist(sapply(names(llmn),strsplit,"_"))
  nlats <- length(lab1)/3
  lats <- as.numeric(lab1[3 * seq(1:nlats) - 1])
  lons <- as.numeric(lab1[3 * seq(1:nlats)])
  pl <- tapply(llmn,list(lats,lons),mean)
  windows(width = 12, height = 8)
  latseq <- as.numeric(rownames(pl))
  lonseq <- as.numeric(colnames(pl))
  if(sp == "ALB") image(lonseq,latseq,t(pl),xlab="Longitude",ylab="Latitude",main=sp, ylim = c(-42, 0), xlim = c(18, 130))
  if(sp == "BET") image(lonseq,latseq,t(pl),xlab="Longitude",ylab="Latitude",main=sp, ylim = c(-42, 25), xlim = c(18, 130))
  if(sp == "YFT") image(lonseq,latseq,t(pl),xlab="Longitude",ylab="Latitude",main=sp, ylim = c(-42, 25), xlim = c(18, 130))
  contour(lonseq,latseq,t(pl),add = TRUE, labcex = 1)
  plot_IO(plot_title = "", uselims = c(20, 130, -50, 25), sp = spreg, newm=F, lwdm=3, axes = F, tcol = 1, mapfill = TRUE)
}

mk_wts <- function(dat, wttype, catch = NULL, sp = NULL, cell_areas = NA) {
  if (wttype == "equal")
    wts <- NULL
  if (wttype == "propn")
    wts <- catch
  # if (wttype == "area") {
  #   a <- tapply(dat$latlong, list(dat$latlong, dat$yrqtr), length)
  #   i <- match(dat$latlong, rownames(a))
  #   j <- match(dat$yrqtr, colnames(a))
  #   n <- mapply("[", list(a), i, j)
  #   wts <- 1/n
  # }
  # if (wttype == "cell_area") {
  #   areas <- cell_areas$garea[match(dat$latlong, cell_areas$latlong)]
  #   a <- tapply(dat$latlong, list(dat$latlong, dat$yrqtr), length)
  #   i <- match(dat$latlong, rownames(a))
  #   j <- match(dat$yrqtr, colnames(a))
  #   n <- mapply("[", list(a), i, j)
  #   wts <- areas/n
  # }
  if (wttype == "area") {
    a <- aggregate(dat$latlong, list(dat$latlong, dat$yrqtr), length)
    i <- match(paste(dat$latlong, dat$yrqtr), paste(a[,1], a[,2]))
    n <- a[i,3]
    wts <- 1/n
  }
  if (wttype == "cell_area") {
    areas <- cell_areas$garea[match(dat$latlong, cell_areas$latlong)]
    a <- aggregate(dat$latlong, list(dat$latlong, dat$yrqtr), length)
    i <- match(paste(dat$latlong, dat$yrqtr), paste(a[,1], a[,2]))
    n <- a[i,3]
    wts <- areas/n
  }
  if (wttype == "catch") {
    if (is.null(catch))
      catch <- tapply(dat[, sp], list(dat$latlong), sum)
    a <- tapply(dat$latlong, list(dat$latlong, dat$yrqtr), length)
    i <- match(dat$latlong, rownames(a))
    j <- match(dat$yrqtr, colnames(a))
    n <- mapply("[", list(a), i, j)
    cwts <- mapply("[", list(catch), i)/sum(catch)
    wts <- cwts/n
  }
  return(wts)
}

#' Map of the Indian Ocean.
#'
#' Function to make a map of the Indian Ocean, with regional boundaries.
#' @param plot_title Plot title.
#' @param uselims Latitudes and Longitudes for the edges of the map.
#' @param sp Region type code for the region boundaries.
#' @param newm If TRUE, create a new plot, otherwise add boundaries etc to existing plot.
#' @param lwdm Line width for boundaries.
#' @param axes If TRUE, create x and y axes.
#' @param tcol Text colour.
#' @param mapfill If TRUE, fill the land area of the map.
#' @param bgc Background colour
#'
plot_IO <- function(plot_title = "", uselims = c(20, 130, -50, 25), sp = "YFT", newm = T, lwdm = 3, axes = T, tcol = "red", mapfill = TRUE, bgc = "lightblue") {
  lims <- uselims
  if (newm) {
    plot(1, 1, yaxt = "n", xaxt = "n", type = "n", xlim = c(lims[1], lims[2]), ylim = c(lims[3], lims[4]), ylab = "", xlab = "", bg = bgc)
    polygon(c(lims[1] - 5, lims[2] + 5, lims[2] + 5, lims[1] - 5), c(lims[3] - 5, lims[3] - 5, lims[4] + 5, lims[4] + 5), col = bgc)
  }
  if (sp == "ALB") {
    lines(c(34.5, 44.2), c(-20, -20), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(49, 120), c(-20, -20), lwd = lwdm, col = "slate grey", lty = 1)
    xoffset <- 5
    yoffset <- 2.5
    text(115 + xoffset, -11 - yoffset, "N", col = tcol, cex = 1.5)
    text(115 + xoffset, -40 - yoffset, "S", col = tcol, cex = 1.5)
  }
  if (sp %in% c("YFT")) {
    lines(c(20, 120), c(-40, -40), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(20, 20), c(-40, -35), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 40), c(-40, -30), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 40), c(-10, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(60, 60), c(-30, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(75, 75), c(-15, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(100, 100), c(-5, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(110, 110), c(-10, -5), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(120, 120), c(-40, -30), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(125, 125), c(-20, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(130, 130), c(-15, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 60), c(-30, -30), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(60, 130), c(-15, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 60), c(-10, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(110, 130), c(-10, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(100, 110), c(-5, -5), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 100), c(10, 10), lwd = lwdm, col = "slate grey", lty = 1)
    text(67.5, 15, "R1", col = tcol, cex = 1.5)
    text(57.5, -2.5, "R2", col = tcol, cex = 1.5)
    text(52.5, -27.5, "R3", col = tcol, cex = 1.5)
    text(82.5, -27.5, "R4", col = tcol, cex = 1.5)
    text(85, -2.5, "R5", col = tcol, cex = 1.5)
    text(90, 15, "R6", col = tcol, cex = 1.5)
  }
  if (sp %in% c("YFT2")) {
    lines(c(20, 120), c(-40, -40), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(20, 20), c(-40, -35), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 40), c(-40, -30), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 40), c(-10, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(60, 60), c(-30, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(75, 75), c(-15, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(100, 100), c(-5, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(110, 110), c(-10, -5), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(120, 120), c(-40, -30), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(125, 125), c(-20, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(130, 130), c(-15, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 60), c(-30, -30), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(60, 130), c(-15, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 60), c(-10, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(110, 130), c(-10, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(100, 110), c(-5, -5), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 75), c(0, 0), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 100), c(10, 10), lwd = lwdm, col = "slate grey", lty = 1)
    text(62.5, 17.5, "R1", col = tcol, cex = 1.5)
    text(62.5, 6, "R2N", col = tcol, cex = 1.5)
    text(62.5, -5, "R2S", col = tcol, cex = 1.5)
    text(52.5, -25, "R3", col = tcol, cex = 1.5)
    text(82.5, -27.5, "R4", col = tcol, cex = 1.5)
    text(85, -2.5, "R5", col = tcol, cex = 1.5)
    text(90, 17.5, "R6", col = tcol, cex = 1.5)
  }
  if (sp %in% c("BET")) {
    lines(c(20, 120), c(-35, -35), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(80, 80), c(-15, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(100, 100), c(-5, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(110, 110), c(-10, -5), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(120, 120), c(-35, -30), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(125, 125), c(-20, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(130, 130), c(-15, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(47, 130), c(-15, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(35, 45), c(-20, -20), lwd = lwdm, col = "slate grey", lty = 1)
    # lines(c(45, 45), c( -20, -15), lwd = lwdm, col = 'slate grey', lty = 1)
    lines(c(110, 130), c(-10, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(100, 110), c(-5, -5), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 100), c(10, 10), lwd = lwdm, col = "slate grey", lty = 1)
    text(65, 0, "R1", col = tcol, cex = 1.5)
    text(90, 0, "R2", col = tcol, cex = 1.5)
    text(70, -30, "R3", col = tcol, cex = 1.5)
  }
  if (sp %in% c("BET3")) {
    lines(c(20, 120), c(-35, -35), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(80, 80), c(-15, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(100, 100), c(-5, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(110, 110), c(-10, -5), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(120, 120), c(-35, -30), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(125, 125), c(-20, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(130, 130), c(-15, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(47, 130), c(-15, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(35, 45), c(-20, -20), lwd = lwdm, col = "slate grey", lty = 1)
    # lines(c(45, 45), c( -20, -15), lwd = lwdm, col = 'slate grey', lty = 1)
    lines(c(110, 130), c(-10, -10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(100, 110), c(-5, -5), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 80), c(0, 0), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(40, 100), c(10, 10), lwd = lwdm, col = "slate grey", lty = 1)
    text(62.5, 6, "R1N", col = tcol, cex = 1.5)
    text(62.5, -5, "R1S", col = tcol, cex = 1.5)
    text(87.5, -2.5, "R2", col = tcol, cex = 1.5)
    text(72.5, -27.5, "R3", col = tcol, cex = 1.5)
  }
  if (sp %in% c("BETcore", "YFTcore")) {
    lines(c(20, 140), c(-35, -35), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(20, 125.5), c(-15, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(20, 100), c(10, 10), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(80, 80), c(10, -15), lwd = lwdm, col = "slate grey", lty = 1)
    xoffset <- 5
    yoffset <- 2.5
    text(115 + xoffset, -11 - yoffset, "N", col = tcol, cex = 1.5)
    text(115 + xoffset, -40 - yoffset, "S", col = tcol, cex = 1.5)
  }
  maps::map("world", yaxt = "n", xaxt = "n", add = T, resolution = 1, interior = F, fill = mapfill)
  if (axes) {
    box(lwd = 3)
    axis(1, at = seq(lims[1], lims[2], by = 10), labels = F)
    axis(2, at = seq(lims[3], lims[4], by = 5), labels = F)
    latseq <- seq(lims[3] + 10, lims[4] - 10, by = 10)
    latseq2 <- as.character(latseq)
    lonseq <- seq(lims[1] + 20, lims[2] - 10, by = 20)
    lonseq2 <- as.character(lonseq)
    latseq2[latseq < 0] <- paste(abs(latseq[latseq < 0]), "S", sep = "")
    latseq2[latseq > 0] <- paste(latseq[latseq > 0], "N", sep = "")
    lonseq2[lonseq < 180] <- paste(lonseq2[lonseq < 180], "E", sep = "")
    lonseq2[lonseq > 180] <- paste(360 - lonseq[lonseq > 180], "W", sep = "")
    axis(2, at = latseq, labels = latseq2, cex.axis = 0.75)
    axis(1, at = lonseq, labels = lonseq2, cex.axis = 0.75)
  }
  mtext(side = 3, line = 0.5, plot_title, font = 2, cex = 1.1)
}

#' Map of the Atlantic Ocean.
#'
#' Function to make a map of the Atlantic Ocean, with regional boundaries.
#' @param plot_title Plot title.
#' @param uselims Latitudes and Longitudes for the edges of the map.
#' @param sp Region type code for the region boundaries.
#' @param newm If TRUE, create a new plot, otherwise add boundaries etc to existing plot.
#' @param lwdm Line width for boundaries.
#' @param axes If TRUE, create x and y axes.
#' @param tcol Text colour.
#'
plot_AO <- function(plot_title = "", uselims = c(-100, 30, -50, 60), sp = "BET", newm = T, lwdm = 3, axes = T, tcol = "red") {
  lims <- uselims
  if (newm) {
    plot(1, 1, yaxt = "n", xaxt = "n", type = "n", xlim = c(lims[1], lims[2]), ylim = c(lims[3], lims[4]), ylab = "", xlab = "", bg = "lightblue")
    polygon(c(lims[1] - 5, lims[2] + 5, lims[2] + 5, lims[1] - 5), c(lims[3] - 5, lims[3] - 5, lims[4] + 5, lims[4] + 5), col = "lightblue")
  }
  if (sp %in% c("BET")) {
    lines(c(-80, -10), c(45, 45), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(-100, -10), c(25, 25), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(-50, 20), c(-15, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(-60, 20), c(-35, -35), lwd = lwdm, col = "slate grey", lty = 1)
    text(-30, 35, "R1", col = tcol, cex = 1.5)
    text(-20, 0, "R2", col = tcol, cex = 1.5)
    text(-10, -25, "R3", col = tcol, cex = 1.5)
  }
  maps::map("world", yaxt = "n", xaxt = "n", add = T, resolution = 1, interior = F, fill = T)
  if (axes) {
    box(lwd = 3)
    axis(1, at = seq(lims[1], lims[2], by = 10), labels = F)
    axis(2, at = seq(lims[3], lims[4], by = 5), labels = F)
    latseq <- seq(lims[3] + 10, lims[4] - 10, by = 10)
    latseq2 <- as.character(latseq)
    lonseq <- seq(lims[1] + 20, lims[2] - 10, by = 20)
    lonseq2 <- as.character(lonseq)
    latseq2[latseq < 0] <- paste(abs(latseq[latseq < 0]), "S", sep = "")
    latseq2[latseq > 0] <- paste(latseq[latseq > 0], "N", sep = "")
    lonseq2[lonseq < 180] <- paste(lonseq2[lonseq < 180], "E", sep = "")
    lonseq2[lonseq > 180] <- paste(360 - lonseq[lonseq > 180], "W", sep = "")
    axis(2, at = latseq, labels = latseq2, cex.axis = 0.75)
    axis(1, at = lonseq, labels = lonseq2, cex.axis = 0.75)
  }
  mtext(side = 3, line = 0.5, plot_title, font = 2, cex = 1.1)
}

#' Map of the Pacific Ocean.
#'
#' Function to make a map of the Pacific Ocean, with regional boundaries.
#' @param plot_title Plot title.
#' @param uselims Latitudes and Longitudes for the edges of the map.
#' @param sp Region type code for the region boundaries.
#' @param newm If TRUE, create a new plot, otherwise add boundaries etc to existing plot.
#' @param lwdm Line width for boundaries.
#' @param axes If TRUE, create x and y axes.
#' @param tcol Text colour.
#'
plot_pacific <- function(plot_title = "", uselims = c(-100, 30, -50, 60), sp = "BET", newm = T, lwdm = 3, axes = T, tcol = "red") {
  lims <- uselims
  if (newm) {
    plot(1, 1, yaxt = "n", xaxt = "n", type = "n", xlim = c(lims[1], lims[2]), ylim = c(lims[3], lims[4]), ylab = "", xlab = "", bg = "lightblue")
    polygon(c(lims[1] - 5, lims[2] + 5, lims[2] + 5, lims[1] - 5), c(lims[3] - 5, lims[3] - 5, lims[4] + 5, lims[4] + 5), col = "lightblue")
  }
  if (sp %in% c("BET")) {
    lines(c(-50, 50), c(210, 210), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(-100, -10), c(25, 25), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(-50, 20), c(-15, -15), lwd = lwdm, col = "slate grey", lty = 1)
    lines(c(-60, 20), c(-35, -35), lwd = lwdm, col = "slate grey", lty = 1)
    text(-30, 35, "R1", col = tcol, cex = 1.5)
    text(-20, 0, "R2", col = tcol, cex = 1.5)
    text(-10, -25, "R3", col = tcol, cex = 1.5)
  }
  maps::map("world", yaxt = "n", xaxt = "n", add = T, resolution = 1, interior = F, fill = T)
  if (axes) {
    box(lwd = 3)
    axis(1, at = seq(lims[1], lims[2], by = 10), labels = F)
    axis(2, at = seq(lims[3], lims[4], by = 5), labels = F)
    latseq <- seq(lims[3] + 10, lims[4] - 10, by = 10)
    latseq2 <- as.character(latseq)
    lonseq <- seq(lims[1] + 20, lims[2] - 10, by = 20)
    lonseq2 <- as.character(lonseq)
    latseq2[latseq < 0] <- paste(abs(latseq[latseq < 0]), "S", sep = "")
    latseq2[latseq > 0] <- paste(latseq[latseq > 0], "N", sep = "")
    lonseq2[lonseq < 180] <- paste(lonseq2[lonseq < 180], "E", sep = "")
    lonseq2[lonseq > 180] <- paste(360 - lonseq[lonseq > 180], "W", sep = "")
    axis(2, at = latseq, labels = latseq2, cex.axis = 0.75)
    axis(1, at = lonseq, labels = lonseq2, cex.axis = 0.75)
  }
  mtext(side = 3, line = 0.5, plot_title, font = 2, cex = 1.1)
}

#' Map of the Pacific Ocean with EPO boundaries.
#'
#' Function to make a map of the Pacific Ocean, with regional boundaries.
#' @param new Make a new plot
#' @param latlim Latitudes for the edges of the map.
#' @param lonlim Longitudes for the edges of the map.
#'
map_EPO <- function(new=FALSE, latlim = c(-40, 40), lonlim = c(140, 290)) {
  if(new) {
    doxax <- "n"
    plot(1:5, 1:5, ylim = latlim, xlim = lonlim, type = "n", xlab = "Longitude", ylab = "Latitude",xaxt=doxax)
  }
  axis(1, at = c(160, 210, 260), labels = c(-200, -150, -100))
  maps::map("world2", add = T, interior=F, fill = TRUE)
  abline(v=210, col = "slate grey", lwd = 2, lty=1)
  lines(c(250, 250), c(-70, 25), lwd = 2, col = "slate grey", lty = 1)
  lines(c(110, 280), c(-10, -10), lwd = 2, col = "slate grey", lty = 1)
  lines(c(210, 250), c(10, 10), lwd = 2, col = "slate grey", lty = 1)
  lines(c(140, 210), c(-40, -40), lwd = 2, col = "slate grey", lty = 1)
  lines(c(140, 150), c(-20, -20), lwd = 2, col = "slate grey", lty = 1)
  lines(c(140, 150), c(-15, -15), lwd = 2, col = "slate grey", lty = 1)
  lines(c(155, 160), c( -5,  -5), lwd = 2, col = "slate grey", lty = 1)
  lines(c(140, 155), c(  0,   0), lwd = 2, col = "slate grey", lty = 1)
  lines(c(140, 210), c( 10,  10), lwd = 2, col = "slate grey", lty = 1)
  lines(c(110, 140), c( 20,  20), lwd = 2, col = "slate grey", lty = 1)
  lines(c(140, 210), c( 50,  50), lwd = 2, col = "slate grey", lty = 1)
  lines(c(110, 110), c(-10,  20), lwd = 2, col = "slate grey", lty = 1)
  lines(c(120, 120), c( 20,  26), lwd = 2, col = "slate grey", lty = 1)
  lines(c(140, 140), c(-40,  20), lwd = 2, col = "slate grey", lty = 1)
  lines(c(150, 150), c(-20, -15), lwd = 2, col = "slate grey", lty = 1)
  lines(c(155, 155), c(-5,    0), lwd = 2, col = "slate grey", lty = 1)
  lines(c(160, 160), c(-10,  -5), lwd = 2, col = "slate grey", lty = 1)
  lines(c(170, 170), c(-40,  50), lwd = 2, col = "slate grey", lty = 1)
  lines(c(210, 210), c(-40,  50), lwd = 2, col = "slate grey", lty = 1)
  text(c(140,190,150,185,155,190,130,143,143), c(25,25,5,0,-30,-30,10,-4.5,-15.6), labels = 1:9, cex=1.6, font = 2, col = 4)
  text(c(230,260,230,260,230), c(0,0,-30,-30,25), labels = c(1:4,0), cex=1.6, font = 2, col = 1)
}
