#' Prepare aggregated data
#'
#' Function to prepare a file of IOTC longline data for use, such as for regional scaling analyses.
#' To set it up for different regions, change the regions parameters in the call to the setup_IO_regions function 
#' @param CE The input data file.
#'
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
  
  loc <- grep(".NO",names(CE))
  CE[,loc][is.na(CE[,loc])] <- 0
  loc <- grep(".MT",names(CE))
  CE[,loc] <- NULL
  
  
  CE <- setup_IO_regions(CE, regY2=TRUE, regB3 = TRUE)
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

#' Allocate statistical weights to each row.
#'
#' The function returns a vector indicating a statistical weight to apply to the row, based on the stratum of which the row is a member.
#' @param dat Input dataset
#' @param wttype Type of statistical weighting method; 'equal' gives equal weight to all rows, 'area' weights rows in inverse proportion to the number of rows per stratum; 'catch' weights rows in proportion to the catch in the stratum divided by the number of rows in the stratum; cell_area weights rows in proportion to the ocean area of the stratum divided by the number of rows in the stratum
#' @param catch An optional vector of catch per stratum.
#' @param sp The species to be used for catch weighting.
#' @param cell_areas An object containing the areas of all cells.
#' @return A vector indicating the statistical weight to apply to the stratum of which the row is a member.
#'
mk_wts <- function(dat, wttype, catch = NULL, sp = NULL, cell_areas = NA) {
  if (wttype == "equal")
    wts <- NULL
  if (wttype == "propn")
    wts <- catch
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

setup_IO_regions <- function(dat, regY=F, regY1=F, regY2=F, regB=F, regB1=F, regB2=F, regB3 = F, regA=F, regA1=F, regA2=F, regA3=F, regA4=F, regA5=F) {
  if(regY) {
    dat$regY <- 0
    dat <- mutate(dat,regY = replace(regY,which(lat5 >=  10 & lon5 < 80 & !is.na(lat5)),1)) %>%
      mutate(regY = replace(regY,which(lat5 <  10 & lat5 >=  -10 & lon5 >= 40 & lon5 < 75 & !is.na(lat5)),2))  %>%
      mutate(regY = replace(regY,which(lat5 <  -10 & lat5 >=  -15 & lon5 >= 60 & lon5 < 75 & !is.na(lat5)),2)) %>%
      mutate(regY = replace(regY,which(lat5 <  -10 & lat5 >= -30 & lon5 >= 20 & lon5 < 60 & !is.na(lat5)),3))  %>%
      mutate(regY = replace(regY,which(lat5 <  -30 & lat5 >= -40 & lon5 >= 20 & lon5 < 40 & !is.na(lat5)),3))  %>%
      mutate(regY = replace(regY,which(lat5 <  -15 & lat5 >= -40 & lon5 >= 60 & lon5 < 120 & !is.na(lat5)),4)) %>%
      mutate(regY = replace(regY,which(lat5 <  -30 & lat5 >= -40 & lon5 >= 40 & lon5 < 60 & !is.na(lat5)),4))  %>%
      mutate(regY = replace(regY,which(lat5 < 10 & lat5 >= -15 & lon5 >= 75 & lon5 < 100 & !is.na(lat5)),5))   %>%
      mutate(regY = replace(regY,which(lat5 < -5 & lat5 >= -15 & lon5 >= 100 & lon5 < 110 & !is.na(lat5)),5))  %>%
      mutate(regY = replace(regY,which(lat5 < -10 & lat5 >= -15 & lon5 >= 110 & lon5 < 130 & !is.na(lat5)),5)) %>%
      mutate(regY = replace(regY,which(lat5 < 30 & lat5 >= 10 & lon5 >= 80 & lon5 < 100 & !is.na(lat5)),6))
  }
  
  if(regY1) {
    dat$regY1 <- 0
    dat <- mutate(dat,regY1 = replace(regY1,which(regY %in% c(1,2)),1)) %>%
      mutate(regY1 = replace(regY1,which(regY %in% c(3)),2)) %>%
      mutate(regY1 = replace(regY1,which(regY %in% c(4)),3)) %>%
      mutate(regY1 = replace(regY1,which(regY %in% c(5,6)),4))
  }
  
  if(regY2) {
    dat$regY2 <- 0
    dat <- mutate(dat,regY2 = replace(regY2,which(lat5 >=  10 & lon5 < 80 & !is.na(lat5)),1)) %>%
      mutate(regY2 = replace(regY2,which(lat5 <  10 & lat5 >=    0 & lon5 >= 40 & lon5 < 75 & !is.na(lat5)),7))  %>%
      mutate(regY2 = replace(regY2,which(lat5 <   0 & lat5 >=  -10 & lon5 >= 40 & lon5 < 75 & !is.na(lat5)),2))  %>%
      mutate(regY2 = replace(regY2,which(lat5 <  -10 & lat5 >=  -15 & lon5 >= 60 & lon5 < 75 & !is.na(lat5)),2)) %>%
      mutate(regY2 = replace(regY2,which(lat5 <  -10 & lat5 >= -30 & lon5 >= 20 & lon5 < 60 & !is.na(lat5)),3))  %>%
      mutate(regY2 = replace(regY2,which(lat5 <  -30 & lat5 >= -40 & lon5 >= 20 & lon5 < 40 & !is.na(lat5)),3))  %>%
      mutate(regY2 = replace(regY2,which(lat5 <  -15 & lat5 >= -40 & lon5 >= 60 & lon5 < 120 & !is.na(lat5)),4)) %>%
      mutate(regY2 = replace(regY2,which(lat5 <  -30 & lat5 >= -40 & lon5 >= 40 & lon5 < 60 & !is.na(lat5)),4))  %>%
      mutate(regY2 = replace(regY2,which(lat5 < 10 & lat5 >= -15 & lon5 >= 75 & lon5 < 100 & !is.na(lat5)),5))   %>%
      mutate(regY2 = replace(regY2,which(lat5 < -5 & lat5 >= -15 & lon5 >= 100 & lon5 < 110 & !is.na(lat5)),5))  %>%
      mutate(regY2 = replace(regY2,which(lat5 < -10 & lat5 >= -15 & lon5 >= 110 & lon5 < 130 & !is.na(lat5)),5)) %>%
      mutate(regY2 = replace(regY2,which(lat5 < 30 & lat5 >= 10 & lon5 >= 80 & lon5 < 100 & !is.na(lat5)),6))
  }
  
  #regB   North of 15S and west of 80 is R1, or north of  20 and west of 45; north of 15S and east of 80 is R2; north of 35S is R3
  if(regB) {
    dat$regB <- 0
    dat <- mutate(dat,regB = replace(regB,which(lat5 <  10 & lat5 >=  -15 & lon5 >= 20 & lon5 < 80 & !is.na(lat5)),1)) %>%
      mutate(regB = replace(regB,which(lat5 <  10 & lat5 >=  -20 & lon5 >= 20 & lon5 < 46 & !is.na(lat5)),1)) %>%
      mutate(regB = replace(regB,which(lat5 < 10 & lat5 >= -15 & lon5 >= 80 & lon5 < 100 & !is.na(lat5)),2)) %>%
      mutate(regB = replace(regB,which(lat5 < -3 & lat5 >= -15 & lon5 >= 100 & lon5 < 110 & !is.na(lat5)),2)) %>%
      mutate(regB = replace(regB,which(lat5 < -7 & lat5 >= -15 & lon5 >= 110 & lon5 < 130 & !is.na(lat5)),2)) %>%
      mutate(regB = replace(regB,which(lat5 <  -20 & lat5 >= -35 & lon5 >= 20 & lon5 < 120 & !is.na(lat5)),3)) %>%
      mutate(regB = replace(regB,which(lat5 <  -15 & lat5 >= -35 & lon5 >= 46 & lon5 < 120 & !is.na(lat5)),3))
  }
  
  if(regB1) {
    dat$regB1 <- 0
    dat <- mutate(dat,regB1 = replace(regB1,which(regB1 %in% 1:3),1))  # Doesn't look right
  }
  
  if(regB2) {
    dat$regB2 <- 0
    dat <- mutate(dat,regB2 = replace(regB2,which(lat5 <  10 & lat5 >=  -15 & lon5 >= 20 & lon5 < 80 & !is.na(lat5)),1)) %>%
      mutate(regB2 = replace(regB2,which(lat5 <  10 & lat5 >=  -20 & lon5 >= 20 & lon5 < 46 & !is.na(lat5)),1)) %>%
      mutate(regB2 = replace(regB2,which(lat5 < 10 & lat5 >= -15 & lon5 >= 80 & lon5 < 100 & !is.na(lat5)),2)) %>%
      mutate(regB2 = replace(regB2,which(lat5 < -3 & lat5 >= -15 & lon5 >= 100 & lon5 < 110 & !is.na(lat5)),2)) %>%
      mutate(regB2 = replace(regB2,which(lat5 < -7 & lat5 >= -15 & lon5 >= 110 & lon5 < 130 & !is.na(lat5)),2)) %>%
      mutate(regB2 = replace(regB2,which(lat5 <  -20 & lat5 >= -35 & lon5 >= 20 & lon5 < 75 & !is.na(lat5)),3)) %>%
      mutate(regB2 = replace(regB2,which(lat5 <  -15 & lat5 >= -35 & lon5 >= 46 & lon5 < 75 & !is.na(lat5)),3)) %>%
      mutate(regB2 = replace(regB2,which(lat5 <  -15 & lat5 >= -35 & lon5 >= 75 & lon5 < 120 & !is.na(lat5)),4))
  }
  
  if(regB3) {
    dat$regB3 <- 0
    dat <- mutate(dat,regB3 = replace(regB3,which(lat5 <  10 & lat5 >=  0 & lon5 >= 20 & lon5 < 80 & !is.na(lat5)),5)) %>%
      mutate(regB3 = replace(regB3,which(lat5 <  0  & lat5 >= -15 & lon5 >= 20 & lon5 < 80 & !is.na(lat5)),1)) %>%
      mutate(regB3 = replace(regB3,which(lat5 <  0 & lat5 >=  -20 & lon5 >= 20 & lon5 < 46 & !is.na(lat5)),1)) %>%
      mutate(regB3 = replace(regB3,which(lat5 < 10 & lat5 >= -15 & lon5 >= 80 & lon5 < 100 & !is.na(lat5)),2)) %>%
      mutate(regB3 = replace(regB3,which(lat5 < -3 & lat5 >= -15 & lon5 >= 100 & lon5 < 110 & !is.na(lat5)),2)) %>%
      mutate(regB3 = replace(regB3,which(lat5 < -7 & lat5 >= -15 & lon5 >= 110 & lon5 < 130 & !is.na(lat5)),2)) %>%
      mutate(regB3 = replace(regB3,which(lat5 <  -20 & lat5 >= -35 & lon5 >= 20 & lon5 < 75 & !is.na(lat5)),3)) %>%
      mutate(regB3 = replace(regB3,which(lat5 <  -15 & lat5 >= -35 & lon5 >= 46 & lon5 < 75 & !is.na(lat5)),3)) %>%
      mutate(regB3 = replace(regB3,which(lat5 <  -15 & lat5 >= -35 & lon5 >= 75 & lon5 < 120 & !is.na(lat5)),4))
  }
  
  #regA
  if(regA) {
    dat$regA <- 0
    dat <- mutate(dat,regA = replace(regA,which(dat$lat5 < -10 & dat$lat5 >=  -20 & !is.na(dat$lat5)),1)) %>%
      mutate(regA = replace(regA,which(dat$lat5 < -20 & dat$lat5 > -40 & !is.na(dat$lat5)),2))
  }
  
  if(regA1) {
    dat$regA1 <- 0
    dat <- mutate(dat,regA1 = replace(regA1,which(dat$lat5 < -10 & dat$lat5 >=  -25 & !is.na(dat$lat5)),1)) %>%
      mutate(regA1 = replace(regA1,which(dat$lat5 < -25 & dat$lat5 > -40 & !is.na(dat$lat5)),2))
  }
  
  if(regA2) {
    dat$regA2 <- 0
    dat <- mutate(dat,regA2 = replace(regA2,which(dat$lat5 < -10 & dat$lat5 >=  -20 & dat$lon5 < 75 & !is.na(dat$lat5)),1)) %>%
      mutate(regA2 = replace(regA2,which(dat$lat5 < -10 & dat$lat5 >=  -20 & dat$lon5 >= 75 & !is.na(dat$lat5)),2)) %>%
      mutate(regA2 = replace(regA2,which(dat$lat5 < -20 & dat$lat5 > -40 & dat$lon5 < 75 & !is.na(dat$lat5)),3)) %>%
      mutate(regA2 = replace(regA2,which(dat$lat5 < -20 & dat$lat5 > -40 & dat$lon5 >= 75 & !is.na(dat$lat5)),4))
  }
  
  if(regA3) {
    dat$regA3 <- 0
    dat <- mutate(dat,regA3 = replace(regA3,which(dat$lat5 <  -10 & dat$lat5 >=  -25 & dat$lon5 < 75 & !is.na(dat$lat5)),1)) %>%
      mutate(regA3 = replace(regA3,which(dat$lat5 <  -10 & dat$lat5 >=  -25 & dat$lon5 >= 75 & !is.na(dat$lat5)),2)) %>%
      mutate(regA3 = replace(regA3,which(dat$lat5 < -25 & dat$lon5 < 75 & !is.na(dat$lat5)),3)) %>%
      mutate(regA3 = replace(regA3,which(dat$lat5 < -25 & dat$lon5 >= 75 & !is.na(dat$lat5)),4))
  }
  
  if(regA4) {
    dat$regA4 <- 0
    dat <- mutate(dat,regA4 = replace(regA4,which(dat$lat5 <  -10 & dat$lat5 >=  -25 & dat$lon5 < 75 & !is.na(dat$lat5)),1)) %>%
      mutate(regA4 = replace(regA4,which(dat$lat5 <  -10 & dat$lat5 >=  -25 & dat$lon5 >= 75 & !is.na(dat$lat5)),2)) %>%
      mutate(regA4 = replace(regA4,which(dat$lat5 < -25 & dat$lat5 > -40 & dat$lon5 < 75 & !is.na(dat$lat5)),3)) %>%
      mutate(regA4 = replace(regA4,which(dat$lat5 < -25 & dat$lat5 > -40 & dat$lon5 >= 75 & !is.na(dat$lat5)),4))
  }
  
  if(regA5) {
    dat$regA5 <- 0
    dat <- mutate(dat,regA5 = replace(regA5,which(dat$lat5 <  -15 & dat$lat5 >=  -45 & dat$lon5 > 55 & dat$lon5 < 100 & !is.na(dat$lat5)),1))
  }
  return(dat)
}

#' Diagnostic plots for Gaussian glm.
#'
#' The function produces 2 Gaussian diagnostic plots: a frequency histogram and a QQ plot.
#' @param res Fitted model for diagnostics.
#' @param ti Titles for the plots.
#'
plotdiags <- function(res, ti = "", ...) {
  hist(res, nclass = 200, freq = F, xlab = "Residuals", main = ti, ...)
  lines(-300:300, dnorm(-300:300, sd = sd(res)), col = 2)
  sdres <- res/sd(res)
  qqDist(sdres, add.median = T)
}

#' QQ plot for Gaussian glm.
#'
#' The function produces a QQ plot for a GLM.
#' @param x Residuals or standardized residuals.
#' @param standardise Standardize the residuals if TRUE.
#' @param add.median Add median line if TRUE.
#' @param ... Other qqnorm() parameters.
#'
qqDist <- function(x, standardise = F, add.median = F, ...) {
  n <- length(x)
  seq.length <- min(1000, n)
  if (standardise) {
    SEQ <- seq(1, 2 * n + 1, length = seq.length)/2
    U <- qnorm(qbeta(0.975, SEQ, rev(SEQ)))
    L <- qnorm(qbeta(0.025, SEQ, rev(SEQ)))
    if (add.median)
      M <- qnorm(qbeta(0.5, SEQ, rev(SEQ)))
  } else {
    SD <- sqrt(var(x) * (n + 1)/n)
    SEQ <- seq(1, 2 * n + 1, length = seq.length)/2
    U <- mean(x) + SD * qt(qbeta(0.975, SEQ, rev(SEQ)), n - 1)
    L <- mean(x) + SD * qt(qbeta(0.025, SEQ, rev(SEQ)), n - 1)
    if (add.median)
      M <- mean(x) + SD * qt(qbeta(0.5, SEQ, rev(SEQ)), n - 1)
  }
  X <- qnorm((SEQ - 0.25)/(n + 0.5))
  qqnorm(x, main = "", ...)
  lines(X, U, type = "l", col = 2)
  lines(X, L, type = "l", col = 2)
  if (add.median)
    lines(X, M, type = "l", col = 2)
  invisible()
}

