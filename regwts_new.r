# Regional weighting
########################################################
#projdir <- "~/../OneDrive_personal/OneDrive/Consulting/IOTC/2018_CPUE/"
projdir <- "~/../OneDrive/Consulting/IOTC/2018_CPUE/"
Rdir <- paste0(projdir, "Rfiles/")
jpdir <- paste0(projdir, "JP/")
krdir <- paste0(projdir, "KR/")
sydir <- paste0(projdir, "SY/")
twdir <- paste0(projdir, "TW/")
jointdir <- paste0(projdir, "joint/")
jntalysis_dir <- paste0(jointdir, "analyses/")
projdir17 <- "~/../OneDrive/Consulting/IOTC/2017_CPUE/"
jointdir17 <- paste0(projdir17, "joint/")
jntalysis_dir17 <- paste0(jointdir17, "analyses/")

regwt_dir <- paste0(jointdir,"regwt/")
regwt_dir2 <- paste0(jointdir,"regwt2/")
#dir.create(regwt_dir2)

setwd(regwt_dir)
setwd(regwt_dir2)

#install.packages("survival")
#install.packages("stringr")
library(stringr)
library("date")
library(splines)
library("maps")
library("mapdata")
library("maptools")
library("lunar")
library("mgcv")
library(randomForest)
library(influ)
library("nFactors")
library(plyr)
library(dplyr)
library(data.table)
library(cluster)
library(beanplot)
library(survival)
source("../../RFiles/support_functions.r")

# Generate data
# dat <- data.frame(x = 1:500,z = runif(500),k = as.factor(sample(c("a","b"),size = 500,replace = TRUE)))
# kvals <- data.frame(kn = c("a","b"),kv = c(20,30))
# dat$y = dat$x + (40*dat$z)^2 + kvals$kv[match(dat$k,kvals$kn)] + rnorm(500,0,30)
# # Fit model
# mod <- glm(y ~ x + ns(z,df = 2) + k,data = dat)
# # Create new dataset
# dat.new <- expand.grid(x = 1:3,z = seq(0.2,0.4,0.1),k = "b")
# # Predict expected values in the usual way
# predict(mod,newdata = dat.new)
# summ <- summary(mod)
# rm(mod)
# # Now, how do I predict using just the summary object and dat.new?
#
#
# # identify models
# # load models
# load(paste0(jntalysis_dir,"std_nocl_hbf/Joint_regB2_R1_lognC_novess_allyrs_summary.RData"))
# BR1 <- summ
# load(paste0(jntalysis_dir,"std_nocl_hbf/Joint_regB2_R2_lognC_novess_allyrs_summary.RData"))
# BR2 <- summ
# load(paste0(jntalysis_dir,"std_nocl_hbf/Joint_regB2_R3_lognC_novess_allyrs_summary.RData"))
# BR3 <- summ
# load(paste0(jntalysis_dir,"std_nocl_hbf/Joint_regB2_R4_lognC_novess_allyrs_summary.RData"))
# BR4 <- summ
# load(paste0(jntalysis_dir,"std_nocl_hbf/Joint_regY_R2_lognC_novess_allyrs_summary.RData"))
# YR2 <- summ
# load(paste0(jntalysis_dir,"std_nocl_hbf/Joint_regY_R5_lognC_novess_allyrs_summary.RData"))
# YR5 <- summ
# load(paste0(jntalysis_dir,"std_nocl_hbf/Joint_regY_R3_lognC_novess_allyrs_summary.RData"))
# YR3 <- summ
# load(paste0(jntalysis_dir,"std_nocl_hbf/Joint_regY_R4_lognC_novess_allyrs_summary.RData"))
# YR4 <- summ
# rm(summ)

# choose a range of years: 1990 - 2000
#sp_hbf <- BR1$coefficients[grep("ns(hbf",rownames(BR1$coefficients),fixed = TRUE)]


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
  CE$latmin <- as.numeric(substring(CE$Grid,3,4)) * sg
  CE$latmax <- CE$latmin + addlat[addtp]
  CE$lonmin <- as.numeric(substring(CE$Grid,5,7))
  CE$lonmax <- CE$lonmin + addlon[addtp]
  CE$centlat <- (CE$latmin + CE$latmax)/2
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
  CE$latlon1 <- as.factor(paste(0.5+CE$latmin,0.5+CE$lonmin, sep = "_"))
  CE$latlon2 <- as.factor(paste(1+2*floor(CE$latmin/2),1+2*floor(CE$lonmin/2), sep = "_"))
  CE$latlon3 <- as.factor(paste(1.5+3*floor(CE$latmin/3),1.5+3*floor(CE$lonmin/3), sep = "_"))

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
  if(sp == "BET") image(lonseq,latseq,t(pl),xlab="Longitude",ylab="Latitude",main=sp, ylim = c(-42, 25), xlim = c(18, 130))
  if(sp == "YFT") image(lonseq,latseq,t(pl),xlab="Longitude",ylab="Latitude",main=sp, ylim = c(-42, 25), xlim = c(18, 130))
  contour(lonseq,latseq,t(pl),add = TRUE, labcex = 1)
  plot_IO(plot_title = "", uselims = c(20, 130, -50, 25), sp = spreg, newm=F, lwdm=3, axes = F, tcol = 1)
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

## Set up data
cefile <- "~/../Google Drive/My papers/IOTC/WPTT/2017-19/IOTC-2017-WPTT19-DATA04_-_CELL.zip"
lldat <- read.csv(unzip(cefile),stringsAsFactors=F)
load(file=paste0(projdir,"cell_areas.RData"))
lldat <- prepCE(lldat)
lldat <- lldat[substring(lldat$Fleet,1,3) %in% c("JPN","KOR"),]
str(lldat)


# Data checking
table(substring(lldat$Grid,1,1), lldat$Fleet)
table(lldat$regB3, useNA = "always")
table(lldat$regY2, useNA = "always")

table(lldat$Fleet)
table(lldat$Fleet)
table(lldat$regY2)
table(lldat$lat5)
table(lldat$MonthEnd - lldat$MonthStart)
table(lldat$QualityCode)

length(unique(lldat$latlong[lldat$regB3 > 0])) # grid cells per regional structure
length(unique(lldat$latlong[lldat$regY2 > 0]))

a <- aggregate(cbind(Effort, BET.NO, YFT.NO) ~ latlong + regB3 + Year, sum, data = lldat)
a <- a[a$regB3 != 0,]
a$regB3[a$regB3  ==  4] <- 3
table(a$regB3)
windows()
aB <- with(a, tapply(Effort, list(Year, regB3), length))
plot(as.numeric(rownames(aB)), aB[,1], type = "l", ylim = c(0, 60), xlim = c(1955, 2018),
     xlab = "Year", ylab = "Number of 5 degree grid cells", main = "BET")
for(r in 2:4) lines(as.numeric(rownames(aB)), aB[,r], col = r)
legend("topright", legend = c("Region 1N", "Region 1S", "Region 2", "Region 3"), lty = 1, col = c(4,1,2,3))
savePlot("Bigeye regB3 strata", type = "png")

a <- aggregate(cbind(Effort, BET.NO, YFT.NO) ~ latlong + regY2 + Year, sum, data = lldat)
a <- a[a$regY2 != 0,]
aY <- with(a, tapply(Effort, list(Year, regY2), length))
windows()
plot(as.numeric(rownames(aY)), aY[,1], type = "l", ylim = c(0, 45), xlim = c(1955, 2018),
     xlab = "Year", ylab = "Number of 5 degree grid cells", main = "YFT")
for(r in 2:7) lines(as.numeric(rownames(aY)), aY[,r], col = r)
legend("topright", legend = c("Region 1", "Region 2N", "Region 2S", "Region 3", "Region 4",
                              "Region 5", "Region 6"), lty = 1, col = c(1,7,2:6))
savePlot("Yellowfin regY2 strata", type = "png")

windows()
names(ll)
lldat$qtr <- lldat$yq - lldat$Year
a <- aggregate(cbind(Effort, BET.NO, YFT.NO) ~ latlong + regY2 + qtr + Year, sum, data = lldat)
a <- a[a$regY2 != 0,]
windows(14, 12); par(mfrow = c(2,3))
for(r in c(2,7,5,3,4)) {
  aa <- a[a$regY2 == r,]
  aY <- with(aa, tapply(Effort, list(Year, qtr), length))
  plot(as.numeric(rownames(aY)), aY[,1], type = "l", ylim = c(0, 45), xlim = c(1955, 2018),
       xlab = "Year", ylab = "Number of 5 degree grid cells", main = paste0("R", r))
  for(q in 2:4) lines(as.numeric(rownames(aY)), aY[,q], col = q)
  legend("topright", legend = c("Q 1", "Q 2", "Q 3", "Q 4"),lty = 1, col = c(1:4))
}
savePlot("Yellowfin strata by qtr and yr", type = "png")



aY[3,]
a[a$regY2 == 6 & a$Year == 1954,]

# starts here
pds <- data.frame(st=c(1960, 1963, 1980), nd=c(1975, 1975, 2000))
pdnm <- as.character(c(6075, 6375, 8000))
allreg <- list(YFT=c(2,3,4,5,7), BET = c(1,2,3,4,5))

# select data
allwts <- list()
sp="BET"; pd = 3

for (pd in 1:3) {
  for (sp in c("YFT", "BET")) {
    ll <- lldat[lldat$yq > pds[pd,"st"] & lldat$yq < pds[pd,"nd"],]
    doreg <- with(allreg, get(sp))
    if (sp == "YFT") {
      ll$sp <- ll$YFT.NO
      ll$reg <- ll$regY2
      spreg = "YFT2"
      }
    if (sp== "BET") {
      ll$sp <- ll$BET.NO
      ll$reg <- ll$regB3
      spreg = "BET3"
    }
    ddd <- ll[ll$reg %in% doreg,]
    ddd$llxx <- as.factor(paste(ddd$reg, ddd$lat5, ddd$lon5, sep = "_"))
    # Remove where not enough effort
    a <- tapply(ddd$Effort, ddd$llxx, sum,na.rm = TRUE)
    a <- a[a > 50000]
    ddd <- ddd[is.na(match(ddd$llxx,names(a))) == F,]
    ##minimum number of quarters
    a <- tapply(ddd$llxx, ddd$llxx, length)
    a <- a[a > 5]
    ddd <- ddd[is.na(match(ddd$llxx,names(a))) == F,]
    ##small number of zero catch records
    ddd$sp[is.na(ddd$sp)] <- 0
    table(ddd$sp > 0, useNA = "always")
#    ddd <- ddd[!is.na(ddd$sp) & ddd$sp > 0,]
    mn <- 0.1 * mean(ddd$sp / ddd$Effort)
    # Calulate stat wts - by area. Weights are (cell area) / number of strata.
    wts <- mk_wts(dat = ddd, wttype = "cell_area", cell_areas = cell_areas)
    ddd <- ddd[!is.na(wts) & wts > 0,]
    wts <- wts[!is.na(wts) & wts > 0]
    wts <- wts / mean(wts, na.rm = TRUE)

    # Model
    ddd$yrqtr <- factor(ddd$yrqtr)
    ddd$llxx <- factor(ddd$llxx)
    ddd$fl <- factor(ddd$Fleet)
    ddd$yr <- as.factor(ddd$Year)
    ddd$qtr <- ddd$yq - ddd$Year
    ddd$reg_qtr <- factor(paste(ddd$reg, ddd$qtr, sep = "_"))
    model23 <- glm(log(sp/Effort + mn) ~ yrqtr + llxx, data = ddd)
    model4 <- glm(log(sp/Effort + mn) ~ yrqtr + llxx, data = ddd, weights = wts) # add statistical weights
    if (length(unique(ddd$fl)) > 1) {
      model5 <- glm(log(sp/Effort + mn) ~ yrqtr + llxx + fl, data = ddd, weights = wts) # include fleet
      model6 <- glm(log(sp/Effort + mn) ~ yr + llxx + fl + reg_qtr, data = ddd, weights = wts) # include quarterly effects
      model7 <- glm(log(sp/Effort + mn) ~ yr + llxx + fl + reg_qtr, data = ddd, weights = wts) # include quarterly effects

    } else {
      model5 <- glm(log(sp/Effort + mn) ~ yrqtr + llxx, data = ddd, weights = wts) # include fleet
      model6 <- glm(log(sp/Effort + mn) ~ yr + llxx + reg_qtr, data = ddd, weights = wts) # include quarterly effects
      model7 <- glm(log(sp/Effort + mn) ~ yr + llxx + reg_qtr, data = ddd, weights = wts) # include quarterly effects
    }
    if (length(unique(ll$fl)) > 1)
    assign(paste0(sp,"_model23"), model23)
    assign(paste0(sp,"_model4"), model4)
    assign(paste0(sp,"_model5"), model5)
    assign(paste0(sp,"_model6"), model6)
    assign(paste0(sp,"_model7"), model7)

    # Predict
    newdat <- expand.grid(llxx = levels(ddd$llxx), yrqtr = levels(ddd$yrqtr))
    newdat$fl <- levels(ddd$fl)[1]
    newdat$cpue23 <- exp(predict.glm(model23, newdata = newdat, type = "response")) - mn
    newdat$cpue4 <- exp(predict.glm(model4, newdata = newdat, type = "response")) - mn
    newdat$cpue5 <- exp(predict.glm(model5, newdata = newdat, type = "response")) - mn
    newdat$cpue5 <- exp(predict.glm(model5, newdata = newdat, type = "response")) - mn
    a <- substring(as.character(newdat$llxx), 3)
    newdat$area <- cell_areas[match(a, cell_areas$latlong), "areax"]
    newdat$area[is.na(newdat$area)] <- 0
    newdat$reg <- substring(as.character(newdat$llxx), 1, 1)

    newdat6 <- expand.grid(llxx = levels(ddd$llxx), qtr = sort(unique(ddd$qtr)), yr = levels(ddd$yr), fl = levels(ddd$fl)[1])
    newdat6$reg <- substring(as.character(newdat6$llxx), 1, 1)
    newdat6$reg_qtr <- paste(newdat6$reg, newdat6$qtr, sep = "_")
    newdat6 <- newdat6[newdat6$reg_qtr %in% unique(ddd$reg_qtr),]    # remove reg_qtr values that can't be predicted because no data were in the model
    newdat6$cpue6 <- exp(predict.glm(model6, newdata = newdat6, type = "response")) - mn
    a <- substring(as.character(newdat6$llxx), 3)
    newdat6$area <- cell_areas[match(a, cell_areas$latlong), "areax"]

#    ddd <- model7$data
    newdat7 <- expand.grid(llxx = levels(ddd$llxx), qtr = sort(unique(ddd$qtr)), yr = levels(ddd$yr), fl = levels(ddd$fl)[1])
    newdat7$reg <- substring(as.character(newdat7$llxx), 1, 1)
    newdat7$reg_qtr <- paste(newdat7$reg, newdat7$qtr, sep = "_")
    newdat7 <- newdat7[newdat7$reg_qtr %in% unique(ddd$reg_qtr),]    # remove reg_qtr values that can't be predicted because no data were in the model
    newdat7$cpue7 <- exp(predict.glm(model7, newdata = newdat7, type = "response")) - mn
    a <- substring(as.character(newdat7$llxx), 3)
    newdat7$area <- cell_areas[match(a, cell_areas$latlong), "areax"]

    llmn23 <- tapply(newdat$cpue23, newdat$llxx, mean)*3000
    regx <- substring(names(llmn23),1,1)
    wts2 <- tapply(llmn23, regx, sum)
    wts2 <- wts2 / max(wts2)
    assign(x = paste0(sp, "_wts2_",pdnm[pd]), value = wts2)

    latlongx <- substring(names(llmn23), 3)
    areax <- cell_areas[match(latlongx, cell_areas$latlong), "areax"]
    areax[is.na(areax)] <- 0
    wts3 <- tapply(llmn23 * areax, regx, sum)
    wts3 <- wts3 / max(wts3)
    assign(x = paste0(sp, "_wts3_",pdnm[pd]), value = wts3)

    llmn4 <- tapply(newdat$cpue4, newdat$llxx, mean)*3000
    regx <- substring(names(llmn4),1,1)
    wts4 <- tapply(llmn4 * areax, regx, sum)
    wts4 <- wts4 / max(wts4)
    assign(x = paste0(sp, "_wts4_",pdnm[pd]), value = wts4)

    llmn5 <- tapply(newdat$cpue5, newdat$llxx, mean)*3000
    regx <- substring(names(llmn5),1,1)
    wts5 <- tapply(llmn5 * areax, regx, sum)
    wts5 <- wts5 / max(wts5)
    assign(x = paste0(sp, "_wts5_",pdnm[pd]), value = wts5)

    llmn6 <- tapply(newdat6$cpue6, newdat6$llxx, mean)*3000
    regx6 <- substring(names(llmn6),1,1)
    latlongx6 <- substring(names(llmn6), 3)
    areax6 <- cell_areas[match(latlongx6, cell_areas$latlong), "areax"]
    areax6[is.na(areax6)] <- 0
    wts6 <- tapply(llmn6 * areax6, regx6, sum)
    wts6 <- wts6 / max(wts6)
    assign(x = paste0(sp, "_wts6_",pdnm[pd]), value = wts6)

    llmn7 <- tapply(newdat7$cpue7, newdat7$llxx, mean)*3000
    regx7 <- substring(names(llmn7),1,1)
    latlongx7 <- substring(names(llmn7), 3)
    areax7 <- cell_areas[match(latlongx7, cell_areas$latlong), "areax"]
    areax7[is.na(areax7)] <- 0
    wts7 <- tapply(llmn7 * areax7, regx7, sum)
    wts7 <- wts7 / max(wts7)
    assign(x = paste0(sp, "_wts7_",pdnm[pd]), value = wts7)

    mcalc1 <- tapply(ddd$sp/ddd$Effort, ddd$llxx, mean)*3000
    regx <- substring(names(mcalc1),1,1)
    mwts <- tapply(mcalc1, regx, sum)
    mwts <- mwts / max(mwts)
    assign(x = paste0(sp, "_wts_mean_",pdnm[pd]), value = mwts/max(mwts))

    plot_patterns(mcalc1, sp, spreg)
    savePlot(paste0("relative wt_means_",pdnm[pd],sp),type = "png")
    plot_patterns(llmn23, sp, spreg)
    savePlot(paste0("relative wt23_",pdnm[pd],sp),type = "png")
    plot_patterns(llmn4, sp, spreg)
    savePlot(paste0("relative wt4_",pdnm[pd],sp),type = "png")
    plot_patterns(llmn5, sp, spreg)
    savePlot(paste0("relative wt5_",pdnm[pd],sp),type = "png")
    plot_patterns(llmn6, sp, spreg)
    savePlot(paste0("relative wt6_",pdnm[pd],sp),type = "png")
    plot_patterns(llmn7, sp, spreg)
    savePlot(paste0("relative wt7_",pdnm[pd],sp),type = "png")
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts1")]] <- mwts/max(mwts)
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts2")]] <- wts2
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts3")]] <- wts3
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts4")]] <- wts4
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts5")]] <- wts5
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts6")]] <- wts6
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts7")]] <- wts7

    tapply(llmn23, regx, mean)
    hist(llmn23[regx>0])
    lim <- quantile(llmn23, 0.95) * 0.1
    llmn23x <- llmn23[llmn23 > lim & regx > 0]
    regx23x <- regx[llmn23 > lim & regx > 0]
    plot_patterns(llmn23x, sp, spreg)

    length(llmn23)
    length(llmn23x)
  }
  # save(BET_wts, BET_wts_mean, file = paste0("BET_wts_",pdnm[pd],".RData"))
  # save(YFT_wts, YFT_wts_mean, file = paste0("YFT_wts_",pdnm[pd],".RData"))
}
save(allwts, file = paste0("allwts.RData"))
load(file = "allwts.RData")

windows(15,15); par(mfrow=c(4,4), mar = c(4,4,3,1))
plotdiags(YFT_model23$residuals, ti = "YFT m2")
plotdiags(BET_model23$residuals, ti = "BET m2")
plotdiags(YFT_model4$residuals, ti = "YFT m4")
plotdiags(BET_model4$residuals, ti = "BET m4")
plotdiags(YFT_model5$residuals, ti = "YFT m5")
plotdiags(BET_model5$residuals, ti = "BET m5")
plotdiags(YFT_model6$residuals, ti = "YFT m6")
plotdiags(BET_model6$residuals, ti = "BET m6")
savePlot("diags_8000", type = "png")

a <- with(YFT_model6, tapply(residuals, list(data$yq, data$reg), mean))
yq <- as.numeric(rownames(a))
windows()
plot(yq, yq, type = "n", ylim = range(a, na.rm = TRUE))
for(r in c(2,3,4,5,7)) {
  lines(yq, a[,r+1], col = r - 1)
}
windows()
boxplot(YFT_model6$residuals ~ YFT_model6$data$Year)

cbind(areax, regx)
barplot(tapply(areax, regx, mean))
barplot(tapply(llmn23[areax==0], regx[areax==0], mean))
barplot(rbind(mwts,wts2,wts3,wts4,wts5,wts6), beside = T)

Anova(YFT_model6, test.statistic = "F")
Anova(BET_model6, test.statistic = "F")
yftdrop6 <- drop1(YFT_model6, test = "LRT")
betdrop6 <- drop1(BET_model6)
yftdrop5 <- drop1(YFT_model5)
betdrop5 <- drop1(BET_model5)
write.csv(yftdrop5, "yftdrop5.csv")
write.csv(betdrop5, "betdrop5.csv")

posy <- grep("YFT2",names(allwts))
posb <- grep("BET3",names(allwts))
allwtsy <- t(as.matrix(data.frame(allwts[posy])))
allwtsb <- t(as.matrix(data.frame(allwts[posb])))
write.csv(allwtsy, "allwt_yft.csv")
write.csv(allwtsb, "allwt_bet.csv")

hist(summary(model6)$coefficients[,4])
# predict relative catch rates for the same kind of effort for all quarters in 1990-2000
# average the prediction across vessels. Choose an HBF or average across. Choose a time period. Sum across areas.
str(lldat)
table(lldat$regB3)
reglist <- data.frame(index=1:8,regnumY=0:7, regnumB=c(0:5, 9,9), regY2=c(0,1,"2S",3,4,5,6,"2N"),regB3=c(0,"1S",2,3,4,"1N",9,9))


windows()
a<- with(allwts, rbind(YFT2_6075_wts6, YFT2_6375_wts6, YFT2_8000_wts6))
a <- a[,c(2,8,3,4,5,6,7)]
colnames(a) <- paste0("R", c(1,"2N", "2S", 3,4,5,6))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = c("1960-1975","1963-1975","1980-2000"), args.legend = list(x="topleft"))
savePlot("barplot_YFT6_pds", type = "png")

windows()
a<- with(allwts, rbind(BET3_6075_wts6, BET3_6375_wts6, BET3_8000_wts6))
a <- a[,c(6,2,3,4,5)]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = c("1960-1975","1963-1975","1980-2000"), args.legend = list(x="topleft"))
savePlot("barplot_BET6_pds", type = "png")

legtxt <- c("m1 means","m2 stdized","m3 areas","m4 stat m","m5 fleet","m6 qtr")

windows()
a<- with(allwts, rbind(YFT2_6375_wts1, YFT2_6375_wts2, YFT2_6375_wts3, YFT2_6375_wts4, YFT2_6375_wts5, YFT2_6375_wts6))
a <- a[,c(NA,2,8,3,4,5,6,7)]
colnames(a) <- c("", paste0("R", c(1,"2N", "2S", 3,4,5,6)))
barplot(a, beside=TRUE, names.arg = colnames(a), legend.text = legtxt, args.legend = list(x="topleft"), main = "1963 - 1975")
savePlot("barplot_YFT_wts_6375", type = "png")

windows()
a<- with(allwts, rbind(YFT2_6075_wts1, YFT2_6075_wts2, YFT2_6075_wts3, YFT2_6075_wts4, YFT2_6075_wts5, YFT2_6075_wts6))
a <- a[,c(NA,2,8,3,4,5,6,7)]
colnames(a) <- c("", paste0("R", c(1,"2N", "2S", 3,4,5,6)))
barplot(a, beside=TRUE, names.arg = colnames(a), legend.text = legtxt, args.legend = list(x="topleft"), main = "1960 - 1975")
savePlot("barplot_YFT_wts_6075", type = "png")

windows()
a<- with(allwts, rbind(YFT2_8000_wts1, YFT2_8000_wts2, YFT2_8000_wts3, YFT2_8000_wts4, YFT2_8000_wts5, YFT2_8000_wts6))
a <- a[,c(2,8,3,4,5,6,7)]
colnames(a) <- c(paste0("R", c(1,"2N", "2S", 3,4,5,6)))
barplot(a, beside=TRUE, names.arg = colnames(a), legend.text = legtxt, args.legend = list(x="topleft"), main = "1980 - 2000")
savePlot("barplot_YFT_wts_8000", type = "png")

windows()
a<- with(allwts, rbind(BET3_6075_wts1,BET3_6075_wts2,BET3_6075_wts3,BET3_6075_wts4,BET3_6075_wts5,BET3_6075_wts6))
a <- a[,c(6,2,3,4,5)]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft"), main = "1960 - 1975")
savePlot("barplot_BET_means_6075", type = "png")

windows()
a<- with(allwts, rbind(BET3_6375_wts1,BET3_6375_wts2,BET3_6375_wts3,BET3_6375_wts4,BET3_6375_wts5,BET3_6375_wts6))
a <- a[,c(6,2,3,4,5)]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft"), main = "1963 - 1975")
savePlot("barplot_BET_means_6375", type = "png")

windows()
a<- with(allwts, rbind(BET3_8000_wts1,BET3_8000_wts2,BET3_8000_wts3,BET3_8000_wts4,BET3_8000_wts5,BET3_8000_wts6))
a <- a[,c(6,2,3,4,5)]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft"), main = "1980 - 2000")
savePlot("barplot_BET_means_8000", type = "png")


dowts <- function(d1, d2, wts, yqs) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yq %in% yqs], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yq %in% yqs], na.rm = TRUE)
  dd <- data.frame(yq=d1$yq,pr=d1$pr * wts[1] + d2$pr[match(d1$yq,d2$yq)] * wts[2])
  dd$pr <- dd$pr/mean(dd$pr,na.rm=TRUE)
  plot(d1$yq,d1$pr,type="l",col=1, xlab = "Year-quarter", ylab = "Relative CPUE")
  lines(d2$yq,d2$pr,col=2)
  lines(dd$yq,dd$pr,col=3)
  return(dd)
}

dowts_3 <- function(d1, d2, d3, wts, yqs) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yq %in% yqs], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yq %in% yqs], na.rm = TRUE)
  d3$pr <- d3$pr / mean(d3$pr[d3$yq %in% yqs], na.rm = TRUE)
  dd <- data.frame(yq=d1$yq,pr=d1$pr * wts[1] + d2$pr[match(d1$yq,d2$yq)] * wts[2] + d3$pr[match(d1$yq,d3$yq)] * wts[3])
  dd$pr <- dd$pr/mean(dd$pr,na.rm=TRUE)
  plot(d1$yq,d1$pr,type="l",col=1, xlab = "Year-quarter", ylab = "Relative CPUE")
  lines(d2$yq,d2$pr,col=2)
  lines(d3$yq,d3$pr,col=3)
  lines(dd$yq,dd$pr,col=4)
  return(dd)
}

dowts_3y <- function(d1, d2, d3, wts, yrs) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yr %in% yrs], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yr %in% yrs], na.rm = TRUE)
  d3$pr <- d3$pr / mean(d3$pr[d3$yr %in% yrs], na.rm = TRUE)
  dd <- data.frame(yr=d1$yr,pr=d1$pr * wts[1] + d2$pr[match(d1$yr,d2$yr)] * wts[2] + d3$pr[match(d1$yr,d3$yr)] * wts[3])
  dd$pr <- dd$pr/mean(dd$pr,na.rm=TRUE)
  plot(d1$yr,d1$pr,type="l",col=1, xlab = "Year-quarter", ylab = "Relative CPUE")
  lines(d2$yr,d2$pr,col=2)
  lines(d3$yr,d3$pr,col=3)
  lines(dd$yr,dd$pr,col=4)
  return(dd)
}

dowts_4 <- function(d1, d2, d3, d4, wts, yqs) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yq %in% yqs], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yq %in% yqs], na.rm = TRUE)
  d3$pr <- d3$pr / mean(d3$pr[d3$yq %in% yqs], na.rm = TRUE)
  d4$pr <- d4$pr / mean(d4$pr[d4$yq %in% yqs], na.rm = TRUE)
  dd <- data.frame(yq=d1$yq,pr=d1$pr * wts[1] + d2$pr[match(d1$yq,d2$yq)] * wts[2] + d3$pr[match(d1$yq,d3$yq)] * wts[3] + d4$pr[match(d1$yq,d4$yq)] * wts[4])
  dd$pr <- dd$pr/mean(dd$pr,na.rm=TRUE)
  plot(d1$yq,d1$pr,type="l",col=1, xlab = "Year-quarter", ylab = "Relative CPUE")
  lines(d2$yq,d2$pr,col=2)
  lines(d3$yq,d3$pr,col=3)
  lines(d4$yq,d4$pr,col=3)
  lines(dd$yq,dd$pr,col=4)
  return(dd)
}

dowts_4y <- function(d1, d2, d3, d4, wts, yrs) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yr %in% yrs], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yr %in% yrs], na.rm = TRUE)
  d3$pr <- d3$pr / mean(d3$pr[d3$yr %in% yrs], na.rm = TRUE)
  d4$pr <- d4$pr / mean(d4$pr[d4$yr %in% yrs], na.rm = TRUE)
  dd <- data.frame(yr=d1$yr,pr=d1$pr * wts[1] + d2$pr[match(d1$yr,d2$yr)] * wts[2] + d3$pr[match(d1$yr,d3$yr)] * wts[3] + d4$pr[match(d1$yr,d4$yr)] * wts[4])
  dd$pr <- dd$pr/mean(dd$pr,na.rm=TRUE)
  plot(d1$yr,d1$pr,type="l",col=1, xlab = "Year-quarter", ylab = "Relative CPUE")
  lines(d2$yr,d2$pr,col=2)
  lines(d3$yr,d3$pr,col=3)
  lines(d4$yr,d4$pr,col=3)
  lines(dd$yr,dd$pr,col=4)
  return(dd)
}


# Combine indices
yqs <- seq(1980.125, 1999.875, .25)
yqs6375 <- seq(1963.125, 1975.875, .25)
bet1 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R1_dellog_boat_allyrs.csv"))
bet2 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R2_dellog_boat_allyrs.csv"))
bet1n <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R5_dellog_boat_allyrs.csv"))
bet1s <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R1_dellog_boat_allyrs.csv"))
betnovess1 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R1_dellog_novess_allyrs.csv"))
betnovess1n <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R5_dellog_novess_allyrs.csv"))
betnovess1s <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R1_dellog_novess_allyrs.csv"))
betnovess2 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R2_dellog_novess_allyrs.csv"))
bet5279_1 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R1_dellog_novess_5279.csv"))
bet5279_1n <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R5_dellog_novess_5279.csv"))
bet5279_1s <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R1_dellog_novess_5279.csv"))
bet5279_2 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R2_dellog_novess_5279.csv"))
bet79nd_1 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R1_dellog_vessid_79nd.csv"))
bet79nd_1n <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R5_dellog_vessid_79nd.csv"))
bet79nd_1s <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R1_dellog_vessid_79nd.csv"))
bet79nd_2 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R2_dellog_vessid_79nd.csv"))
# windows(12,12); par(mfrow = c(2,2))
# wts <- allwts$BET3_8000[c(5,1,2) + 1]
# wtsx <- allwts$BET3_6375[c(5,1,2) + 1]
# dd <- dowts(bet1, bet2, wts, yqs)
# write.csv(dd,"Joint_regB2_R12_dellog_boat_allyrs.csv")
# dd <- dowts(betnovess1, betnovess2, wts, yqs)
# write.csv(dd,"Joint_regB2_R12_dellog_novess_allyrs.csv")
# dd <- dowts(bet5279_1, bet5279_2, wtsx, yqs = yqs6375)
# write.csv(dd,"Joint_regB2_R12_dellog_novess_5279.csv")
# dd <- dowts(bet79nd_1, bet79nd_2, wts, yqs)
# write.csv(dd,"Joint_regB2_R12_dellog_vessid_79nd.csv")
# savePlot("Joint_regB2_R12_dellog", type = "png")

wts2 <- allwts$BET3_8000_wts6[c(5,1,2) + 1]
#wts <- allwts$BET3_8000[c(5,1,2) + 1]
wts2_early <- allwts$BET3_6375_wts6[c(5,1,2) + 1]
windows(12,12); par(mfrow = c(2,2))
dd <- dowts_3(bet1n, bet1s, bet2, wts2, yqs)
write.csv(dd,"Joint_regB3_R12_dellog_boat_allyrs.csv")
dd <- dowts_3(betnovess1n, betnovess1s, betnovess2, wts2, yqs)
write.csv(dd,"Joint_regB3_R12_dellog_novess_allyrs.csv")
dd <- dowts_3(bet5279_1n, bet5279_1s, bet5279_2, wts2_early, yqs = yqs6375)
write.csv(dd,"Joint_regB3_R12_dellog_novess_5279.csv")
dd <- dowts_3(bet79nd_1n, bet79nd_1s, bet79nd_2, wts2, yqs)
write.csv(dd,"Joint_regB3_R12_dellog_vessid_79nd.csv")
savePlot("Joint_regB3_R12_dellog", type = "png")

yrs <- seq(1980, 1999, 1)
yrs6375 <- seq(1963, 1975, 1)
bet1y <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R1_dellog_boat_allyrsyr.csv"))
bet2y <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R2_dellog_boat_allyrsyr.csv"))
bet1ny <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R5_dellog_boat_allyrsyr.csv"))
bet1sy <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R1_dellog_boat_allyrsyr.csv"))
betnovess1y <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R1_dellog_novess_allyrsyr.csv"))
betnovess1ny <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R5_dellog_novess_allyrsyr.csv"))
betnovess1sy <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R1_dellog_novess_allyrsyr.csv"))
betnovess2y <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R2_dellog_novess_allyrsyr.csv"))
bet5279_1y <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R1_dellog_novess_5279yr.csv"))
bet5279_1ny <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R5_dellog_novess_5279yr.csv"))
bet5279_1sy <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R1_dellog_novess_5279yr.csv"))
bet5279_2y <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R2_dellog_novess_5279yr.csv"))
bet79nd_1y <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R1_dellog_vessid_79ndyr.csv"))
bet79nd_1ny <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R5_dellog_vessid_79ndyr.csv"))
bet79nd_1sy <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R1_dellog_vessid_79ndyr.csv"))
bet79nd_2y <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R2_dellog_vessid_79ndyr.csv"))
wts6 <- allwts$BET3_8000_wts6[c(5,1,2) + 1]
wts6_early <- allwts$BET3_6375_wts6[c(5,1,2) + 1]
windows(12,12); par(mfrow = c(2,2))
dd <- dowts_3y(bet1ny, bet1sy, bet2y, wts6, yrs)
write.csv(dd,"Joint_regB3_R12_dellog_boat_allyrsyr.csv")
dd <- dowts_3y(betnovess1ny, betnovess1sy, betnovess2y, wts6, yrs)
write.csv(dd,"Joint_regB3_R12_dellog_novess_allyrsyr.csv")
dd <- dowts_3y(bet5279_1ny, bet5279_1sy, bet5279_2y, wts6_early, yrs = yrs6375)
write.csv(dd,"Joint_regB3_R12_dellog_novess_5279yr.csv")
dd <- dowts_3y(bet79nd_1ny, bet79nd_1sy, bet79nd_2y, wts6, yrs)
write.csv(dd,"Joint_regB3_R12_dellog_vessid_79ndyr.csv")
savePlot("Joint_regB3_R12_dellogy", type = "png")



# Combine indices
yqs <- seq(1980.125, 1999.875, .25)
yqs6375 <- seq(1963.125, 1975.875, .25)
yft5 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_boat_allyrs_yq.csv"))
yft2 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_boat_allyrs_yq.csv"))
yft3 <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_boat_allyrs_yq.csv"))
yft4 <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R4_dellog_boat_allyrs_yq.csv"))
yft2n <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_boat_allyrs_yq.csv"))
yft2s <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_boat_allyrs_yq.csv"))
yftnovess5 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_novess_allyrs_yq.csv"))
yftnovess2 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_novess_allyrs_yq.csv"))
yftnovess2n <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_novess_allyrs_yq.csv"))
yftnovess2s <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_novess_allyrs_yq.csv"))
yftnovess3 <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_novess_allyrs_yq.csv"))
yftnovess4 <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R4_dellog_novess_allyrs_yq.csv"))
yft5279_5 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_novess_5279_yq.csv"))
yft5279_2 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_novess_5279_yq.csv"))
yft5279_2n <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_novess_5279_yq.csv"))
yft5279_2s <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_novess_5279_yq.csv"))
yft5279_3 <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_novess_5279_yq.csv"))
yft5279_4 <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R4_dellog_novess_5279_yq.csv"))
yft79nd_5 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_vessid_79nd_yq.csv"))
yft79nd_2 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_vessid_79nd_yq.csv"))
yft79nd_2n <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_vessid_79nd_yq.csv"))
yft79nd_2s <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_vessid_79nd_yq.csv"))
yft79nd_3 <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_vessid_79nd_yq.csv"))
yft79nd_4 <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R4_dellog_vessid_79nd_yq.csv"))
wts6 <- allwts$YFT2_8000_wts6[c(7,2,5,3) + 1]
wts6_early <- allwts$YFT2_6375_wts6[c(7,2,3,5) + 1]
windows(12,12); par(mfrow = c(2,2))
dd <- dowts_4(yft2n, yft2s, yft5, yft3, wts6, yqs)
write.csv(dd,"Joint_regY2_R235_dellog_boat_allyrs_yq.csv")
dd <- dowts_4(yftnovess2n, yftnovess2s, yftnovess5, yftnovess3, wts6, yqs)
write.csv(dd,"Joint_regY2_R235_dellog_novess_allyrs_yq.csv")
dd <- dowts_4(yft5279_2n, yft5279_2s, yft5279_5, yft5279_3, wts6_early, yqs = yqs6375)
write.csv(dd,"Joint_regY2_R235_dellog_novess_5279_yq.csv")
dd <- dowts_4(yft79nd_2n, yft79nd_2s, yft79nd_5, yft79nd_3, wts6, yqs)
write.csv(dd,"Joint_regY2_R235_dellog_vessid_79nd_yq.csv")
savePlot("Joint_regY2_R235_dellog", type = "png")

yrs <- seq(1980, 1999, 1)
yrs6375 <- seq(1963, 1975, 1)
yft5y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_boat_allyrs_yr.csv"))
yft2y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_boat_allyrs_yr.csv"))
yft3y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_boat_allyrs_yr.csv"))
yft2ny <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_boat_allyrs_yr.csv"))
yft2sy <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_boat_allyrs_yr.csv"))
yftnovess5y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_novess_allyrs_yr.csv"))
yftnovess2y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_novess_allyrs_yr.csv"))
yftnovess2ny <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_novess_allyrs_yr.csv"))
yftnovess2sy <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_novess_allyrs_yr.csv"))
yftnovess3y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_novess_allyrs_yr.csv"))
yft5279_5y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_novess_5279_yr.csv"))
yft5279_2y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_novess_5279_yr.csv"))
yft5279_2ny <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_novess_5279_yr.csv"))
yft5279_2sy <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_novess_5279_yr.csv"))
yft5279_3y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_novess_5279_yr.csv"))
yft79nd_5y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_vessid_79nd_yr.csv"))
yft79nd_2y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_vessid_79nd_yr.csv"))
yft79nd_2ny <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_vessid_79nd_yr.csv"))
yft79nd_2sy <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_vessid_79nd_yr.csv"))
yft79nd_3y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_vessid_79nd_yr.csv"))
wts6 <- allwts$YFT2_8000_wts6[c(7,2,5,3) + 1]
wts6_early <- allwts$YFT2_6375_wts6[c(7,2,3,5) + 1]
windows(12,12); par(mfrow = c(2,2))
dd <- dowts_4y(yft2ny, yft2sy, yft5y, yft3y, wts6, yrs)
write.csv(dd,"Joint_regY2_R235_dellog_boat_allyrs_yr.csv")
dd <- dowts_4y(yftnovess2ny, yftnovess2sy, yftnovess5y, yftnovess3y, wts6, yrs)
write.csv(dd,"Joint_regY2_R235_dellog_novess_allyrs_yr.csv")
dd <- dowts_4y(yft5279_2ny, yft5279_2sy, yft5279_5y, yft5279_3y, wts6_early, yrs = yrs6375)
write.csv(dd,"Joint_regY2_R235_dellog_novess_5279_yr.csv")
dd <- dowts_4y(yft79nd_2ny, yft79nd_2sy, yft79nd_5y, yft79nd_3y, wts6, yrs)
write.csv(dd,"Joint_regY2_R235_dellog_vessid_79nd_yr.csv")
savePlot("Joint_regY2_R235_dellogyr", type = "png")


# Combine indices
yq8000 <- seq(1980.125, 1999.875, .25)
yr8000 <- 1980:1999
yq6375 <- seq(1963.125, 1974.875, .25)
yq6075 <- seq(1960.125, 1974.875, .25)

yft2n <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_novess_allyrs_yq.csv"))
yft2s <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_novess_allyrs_yq.csv"))
yft5 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_novess_allyrs_yq.csv"))
yft3 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R3_dellog_novess_allyrs_yq.csv"))
yft4 <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R4_dellog_novess_allyrs_yq.csv"))
bet1n <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R5_dellog_novess_allyrs.csv"))
bet1s <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf_spl/outputs/Joint_regB3_R1_dellog_novess_allyrs.csv"))
bet2 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R2_dellog_novess_allyrs.csv"))
bet3 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R3_dellog_novess_allyrs.csv"))
bet4 <- read.csv(paste0(jntalysis_dir17,"std_nocl_hbf/outputs/Joint_regB2_R4_dellog_novess_allyrs.csv"))
# add together to make southern BET. Add proportions based on estiates.
betsth <- dowts(bet3, bet4, allwts$BET3_8000_wts6[c(5,6)], yq8000)



betlist <- list(bet1n, bet1s, bet2, betsth)
betlist_ord <- c(5,1,2,3,4) + 1
yftlist <- list(yft2n, yft2s, yft3, yft4, yft5)
yftlist_ord <- c(7,2,3,4,5) + 1
length(betlist)
length(yftlist)

dowts_list <- function(dlist, wts, yqs) { # Note that the regions in dlist & wts must be in the same sequence
  for(i in 1:length(dlist)) {
    a <- dlist[[i]]
    a$prwt <- wts[i] * a$pr / mean(a$pr[a$yq %in% yqs], na.rm = TRUE)
    dlist[[i]] <- a
  }
  return(dlist)
}
load(file = "allwts.RData")

# Calulate means for each species and period
rng_mn <- function(x,y) mean(x$pr[x$yq %in% y], na.rm = TRUE)
brng6075 <- c(rng_mn(bet1n, yq6075),rng_mn(bet1s, yq6075),rng_mn(bet2, yq6075),rng_mn(betsth, yq6075))
brng6375 <- c(rng_mn(bet1n, yq6375),rng_mn(bet1s, yq6375),rng_mn(bet2, yq6375),rng_mn(betsth, yq6375))
brng8000 <- c(rng_mn(bet1n, yq8000),rng_mn(bet1s, yq8000),rng_mn(bet2, yq8000),rng_mn(betsth, yq8000))
yrng6075 <- c(rng_mn(yft2n, yq6075),rng_mn(yft2s, yq6075),rng_mn(yft3, yq6075),rng_mn(yft4, yq6075),rng_mn(yft5, yq6075))
yrng6375 <- c(rng_mn(yft2n, yq6375),rng_mn(yft2s, yq6375),rng_mn(yft3, yq6375),rng_mn(yft4, yq6375),rng_mn(yft5, yq6375))
yrng8000 <- c(rng_mn(yft2n, yq8000),rng_mn(yft2s, yq8000),rng_mn(yft3, yq8000),rng_mn(yft4, yq8000),rng_mn(yft5, yq8000))
brng6075; brng6375; brng8000
yrng6075; yrng6375; yrng8000

# Make adjusted series
bet8000st <- dowts_list(dlist = betlist, wts = c(allwts$BET3_8000_wts6[betlist_ord[1:3]], sum(allwts$BET3_8000_wts6[betlist_ord[4:5]])), yqs = yqs)
bet6375st <- dowts_list(dlist = betlist, wts = c(allwts$BET3_6375_wts6[betlist_ord[1:3]], sum(allwts$BET3_6375_wts6[betlist_ord[4:5]])), yqs = yqs)
yft8000st <- dowts_list(dlist = yftlist, wts = allwts$YFT2_8000_wts6[yftlist_ord], yqs = yqs)
yft6375st <- dowts_list(dlist = yftlist, wts = allwts$YFT2_6375_wts6[yftlist_ord], yqs = yqs6375)

# Barplots comparing each wt between periods
mkwt <- function(wtobj, wtnum) paste0(wtobj, "_wts", wtnum)

for (wt in 1:6) {
  a<- with(allwts, rbind(get(mkwt("YFT2_6075", wt)),get(mkwt("YFT2_6375", wt)),get(mkwt("YFT2_8000", wt))))
  a <- a[, yftlist_ord]
  a[1,] <- a[1,] / yrng6075
  a[2,] <- a[2,] / yrng6375
  a[3,] <- a[3,] / yrng8000
  a <- a / apply(a,1,max)
  colnames(a) <- paste0("R", c("2N", "2S", 3,4,5))
  windows()
  barplot(a, beside = TRUE, names.arg = colnames(a), legend = c("1960-1975","1963-1975","1980-2000"), args.legend = list(x = "topleft"), main = paste0("YFT m", wt))
  savePlot(paste0("barplot_YFTadj_pds_wts", wt), type = "png")

  a<- with(allwts, rbind(get(mkwt("BET3_6075", wt)),get(mkwt("BET3_6375", wt)),get(mkwt("BET3_8000", wt))))
  a <- a[,betlist_ord]
  sth <- apply(a[,4:5],1,sum)
  a <- cbind(a[,1:3], sth)
  a[1,] <- a[1,] / brng6075
  a[2,] <- a[2,] / brng6375
  a[3,] <- a[3,] / brng8000
  a <- a / apply(a,1,max)
  colnames(a) <- paste0("R", c("1N", "1S", 2,3))
  windows()
  barplot(a, beside = TRUE, names.arg = colnames(a), legend = c("1960-1975","1963-1975","1980-2000"), args.legend = list(x = "topleft"), main = paste0("BET m", wt))
  savePlot(paste0("barplot_BETadj_pds_wts", wt), type = "png")
}


windows(12,10); par(mfrow = c(3,2))
for(i in 1:5) {
  with(bet8000st[[i]], plot(yq, prwt, xlab = "Year-quarter", ylab = "Weighted CPUE", col= 1, type = "l", main = betlab[i]))
  with(bet6375st[[i]], lines(yq, prwt, col = 2, lty = 2))
}


dd <- dowts(bet1, bet2, BET_wts[-1], yqs)
write.csv(dd,"Joint_regB2_R12_dellog_boat_allyrs.csv")
dd <- dowts(betnovess1, betnovess2, BET_wts[-1], yqs)
write.csv(dd,"Joint_regB2_R12_dellog_novess_allyrs.csv")
dd <- dowts(bet5279_1, bet5279_2, BET_wts[-1], yqs = yqsfake)
write.csv(dd,"Joint_regB2_R12_dellog_novess_5279.csv")
dd <- dowts(bet79nd_1, bet79nd_2, BET_wts[-1], yqs)
write.csv(dd,"Joint_regB2_R12_dellog_vessid_79nd.csv")
savePlot("Joint_regB2_R12_dellog", type = "png")

windows()
a <- with(allwts, rbind(YFT2_6375_wts1, YFT2_6375_wts2, YFT2_6375_wts3, YFT2_6375_wts4, YFT2_6375_wts5, YFT2_6375_wts6))
a <- a[,c(NA,8,3,4,5,6)]
a[,-1] <- t(t(a[,-1]) / yrng6375)
a <- a / apply(a,1,max, na.rm = TRUE)
colnames(a) <- c("", paste0("R", c("2N", "2S", 3,4,5)))
barplot(a, beside=TRUE, names.arg = colnames(a), legend.text = legtxt, args.legend = list(x="topleft"), main = "YFT 1963 - 1975")
savePlot("barplotadj_YFT_wts_6375", type = "png")

windows()
a<- with(allwts, rbind(YFT2_6075_wts1, YFT2_6075_wts2, YFT2_6075_wts3, YFT2_6075_wts4, YFT2_6075_wts5, YFT2_6075_wts6))
a <- a[,c(NA,8,3,4,5,6)]
a[,-1] <- t(t(a[,-1]) / yrng6075)
a <- a / apply(a,1,max, na.rm = TRUE)
colnames(a) <- c("", paste0("R", c("2N", "2S", 3,4,5)))
barplot(a, beside=TRUE, names.arg = colnames(a), legend.text = legtxt, args.legend = list(x="topleft"), main = "YFT 1960 - 1975")
savePlot("barplotadj_YFT_wts_6075", type = "png")

windows()
a<- with(allwts, rbind(YFT2_8000_wts1, YFT2_8000_wts2, YFT2_8000_wts3, YFT2_8000_wts4, YFT2_8000_wts5, YFT2_8000_wts6))
a <- a[,c(NA,8,3,4,5,6)]
a[,-1] <- t(t(a[,-1]) / yrng8000)
a <- a / apply(a,1,max, na.rm = TRUE)
colnames(a) <- c("", paste0("R", c("2N", "2S", 3,4,5)))
barplot(a, beside=TRUE, names.arg = colnames(a), legend.text = legtxt, args.legend = list(x="topleft"), main = "YFT 1980 - 2000")
savePlot("barplotadj_YFT_wts_8000", type = "png")

windows()
a<- with(allwts, rbind(BET3_6075_wts1,BET3_6075_wts2,BET3_6075_wts3,BET3_6075_wts4,BET3_6075_wts5,BET3_6075_wts6))
a <- a[,c(6,2,3,4,5)]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- t(t(a) / brng6075)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft"), main = "BET 1960 - 1975")
savePlot("barplotadj_BET_means_6075", type = "png")

windows()
a<- with(allwts, rbind(BET3_6375_wts1,BET3_6375_wts2,BET3_6375_wts3,BET3_6375_wts4,BET3_6375_wts5,BET3_6375_wts6))
a <- a[,c(6,2,3,4,5)]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- t(t(a) / brng6375)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft"), main = "BET 1963 - 1975")
savePlot("barplotadj_BET_means_6375", type = "png")

windows()
a<- with(allwts, rbind(BET3_8000_wts1,BET3_8000_wts2,BET3_8000_wts3,BET3_8000_wts4,BET3_8000_wts5,BET3_8000_wts6))
a <- a[,c(6,2,3,4,5)]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- t(t(a) / brng8000)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft"), main = "BET 1980 - 2000")
savePlot("barplotadj_BET_means_8000", type = "png")
