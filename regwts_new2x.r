# Regional weighting
########################################################

library(mgcv)
library(here)
library(tidyverse)

repodir <- paste0(here(),"/")
projdir <- paste0(repodir,"..")

regscale_wkdir <- paste0(projdir,"/regscale_wkdir")
dir.create(regscale_wkdir)
setwd(regscale_wkdir)

source(paste0(repodir, "support_functions.r"))


windows(width=10, height = 10)
plot_IO(plot_title = "", uselims = c(20, 130, -50, 25), sp = "YFT2", newm=TRUE, lwdm=3, axes = T, tcol = 1, mapfill = TRUE, bgc = "lightgrey")
savePlot("map_YFT_regions", type = "png")
plot_IO(plot_title = "", uselims = c(20, 130, -50, 25), sp = "BET3", newm=TRUE, lwdm=3, axes = T, tcol = 1, mapfill = TRUE, bgc = "lightgrey")
savePlot("map_BET_regions", type = "png")


## Set up data
cefile <- "IOTC-2020-WPTT22-DATA04-CELongline.zip" # The IOTC 2020 data file. Data for the scaling period is unlikely to change. 
lldat <- read.csv(unzip(cefile),stringsAsFactors=F)
lldat <- prepCE(lldat)
lldat <- lldat[substring(lldat$Fleet,1,3) %in% c("JPN","KOR"),]
str(lldat)
load(file=paste0(repodir,"cell_areas.RData"))


# Data checking
table(substring(lldat$Grid,1,1), lldat$Fleet) # make sure it's getting data
table(lldat$regB3, useNA = "always") 
table(lldat$regY2, useNA = "always")
table(lldat$lat5)

length(unique(lldat$latlong[lldat$regB3 > 0])) # grid cells per regional structure
length(unique(lldat$latlong[lldat$regY2 > 0]))

# figure showing number of 5 degree grid cells per region - BET
a <- aggregate(cbind(Effort, BET.NO, YFT.NO) ~ latlong + regB3 + Year, sum, data = lldat)
a <- a[a$regB3 != 0,]
a$regB3[a$regB3  ==  4] <- 3 # for plotting, merge region 4 (southeast) into region 3 (southwest) 
table(a$regB3)
windows()
aB <- with(a, tapply(Effort, list(Year, regB3), length))
plot(as.numeric(rownames(aB)), aB[,1], type = "l", ylim = c(0, 70), xlim = c(1955, 2018),
     xlab = "Year", ylab = "Number of 5 degree grid cells", main = "BET")
for(r in 2:4) lines(as.numeric(rownames(aB)), aB[,r], col = r)
legend("topright", legend = c("Region 1N", "Region 1S", "Region 2", "Region 3"), lty = 1, col = c(4,1,2,3))
savePlot("Bigeye regB3 cells", type = "png")

# figure showing number of 5 degree grid cells per region - YFT
a <- aggregate(cbind(Effort, BET.NO, YFT.NO) ~ latlong + regY2 + Year, sum, data = lldat)
a <- a[a$regY2 != 0,]
aY <- with(a, tapply(Effort, list(Year, regY2), length))
windows()
plot(as.numeric(rownames(aY)), aY[,1], type = "l", ylim = c(0, 55), xlim = c(1955, 2018),
     xlab = "Year", ylab = "Number of 5 degree grid cells", main = "YFT")
for(r in 2:7) lines(as.numeric(rownames(aY)), aY[,r], col = r)
legend("topright", legend = c("Region 1", "Region 2N", "Region 2S", "Region 3", "Region 4",
                              "Region 5", "Region 6"), lty = 1, col = c(1,7,2:6))
savePlot("Yellowfin regY2 cells", type = "png")

# Figure showing number of 5 degree grid cells per region and qtr - YFT
names(lldat)
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
savePlot("Yellowfin cells by qtr and yr", type = "png")

# Figure showing number of 5 degree grid cells per region and qtr - YFT
a <- aggregate(cbind(Effort, BET.NO, YFT.NO) ~ latlong + regY2 + qtr + MonthStart + Year, sum, data = lldat)
a <- a[a$regY2 != 0,]
windows(14, 12); par(mfrow = c(2,3))
for(r in c(2,7,5,3,4)) {
  aa <- a[a$regY2 == r,]
  aY <- with(aa, tapply(Effort, list(Year, qtr), length))
  plot(as.numeric(rownames(aY)), aY[,1], type = "l", ylim = c(0, 90), xlim = c(1955, 2018),
       xlab = "Year", ylab = "Number of 5 degree grid cells", main = paste0("R", r))
  for(q in 2:4) lines(as.numeric(rownames(aY)), aY[,q], col = q)
  legend("topright", legend = c("Q 1", "Q 2", "Q 3", "Q 4"),lty = 1, col = c(1:4))
}
savePlot("Yellowfin strata by qtr and yr", type = "png")


# starts here
pds <- data.frame(st=c(1960, 1963, 1975, 1979, 1980), nd=c(1975, 1975, 1994, 1994, 2000)) # choose options for scaling periods
pdnm <- as.character(c(6075, 6375, 7594, 7994, 8000)) # name the alternative scaling periods
allreg <- list(YFT=c(1,2,3,4,5,6,7), BET = c(1,2,3,4,5)) # regions to include

# select data
allwts <- list()
# sp="YFT"; pd = 4 # options for testing

for (pd in 1:5) {
  for (sp in c("YFT", "BET")) {
    ll <- lldat[lldat$yq > pds[pd,"st"] & lldat$yq < pds[pd,"nd"],]
    doreg <- with(allreg, get(sp)) # get the list of regions to run
    if (sp == "YFT") {
      ll$sp <- ll$YFT.NO
      spreg = "YFT2"
      spreg2 <- "regY2"
    }
    if (sp== "BET") {
      ll$sp <- ll$BET.NO
      spreg = "BET3"
      spreg2 <- "regB3"
    }
    ll$reg <- ll[,spreg2]
    ddd <- ll[ll$reg %in% doreg,]
    ddd$llxx <- as.factor(paste(ddd$reg, ddd$lat5, ddd$lon5, sep = "_"))
    # Remove where not enough effort
    a <- tapply(ddd$Effort, ddd$llxx, sum,na.rm = TRUE)
    a <- a[a > 50000]
    ddd <- ddd[is.na(match(ddd$llxx,names(a))) == F,]
    ##minimum number of strata
    a <- tapply(ddd$llxx, ddd$llxx, length)
    a <- a[a > 6]
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
      model7 <- gam(log(sp/Effort + mn) ~ yr + te(lat5, lon5) + fl + reg_qtr, data = ddd, weights = wts) # gam

    } else {
      model5 <- glm(log(sp/Effort + mn) ~ yrqtr + llxx, data = ddd, weights = wts) # include fleet
      model6 <- glm(log(sp/Effort + mn) ~ yr + llxx + reg_qtr, data = ddd, weights = wts) # include quarterly effects
      model7 <- gam(log(sp/Effort + mn) ~ yr + te(lat5, lon5) + reg_qtr, data = ddd, weights = wts) # gam
      }
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

    # include all llxx, plus all lat*lon interactions where lat5 < 0 and lon5 > 40 & lon5 < 100
    newdat7 <- expand.grid(lat5 = sort(unique(ddd$lat5)), lon5 = sort(unique(ddd$lon5)), qtr = sort(unique(ddd$qtr)), yr = levels(ddd$yr), fl = levels(ddd$fl)[1])
    newdat7 <- setup_IO_regions(newdat7, regY2=TRUE, regB3 = TRUE)
    newdat7$reg <- newdat7[,spreg2]
    newdat7 <- newdat7[newdat7$reg %in% doreg,]
    newdat7$llxx <- paste(newdat7$reg, newdat7$lat5, newdat7$lon5, sep = "_")
    newdat7$reg_qtr <- paste(newdat7$reg, newdat7$qtr, sep = "_")
    newdat7 <- newdat7[newdat7$reg_qtr %in% sort(unique(ddd$reg_qtr)),]    # remove reg_qtr values that can't be predicted because no data were in the model
    newdat7$cpue7 <- exp(predict.gam(model7, newdata = newdat7, type = "response")) - mn
    a <- substring(as.character(newdat7$llxx), 3)
    newdat7$area <- cell_areas[match(a, cell_areas$latlong), "areax"]
    newdat7 <- newdat7[newdat7$llxx %in% newdat6$llxx | (newdat7$lat5 < 0 & newdat7$lon5 > 40 & newdat7$lon5 < 100),]

    a <- newdat7[!newdat7$llxx %in% newdat6$llxx,]
    colnames(a)[colnames(a)=="cpue7"] <- "cpue6"
    newdat8 <- rbind(newdat6, a[,colnames(newdat6)])
    colnames(newdat8)[colnames(newdat8)=="cpue6"] <- "cpue8"

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

    llmn8 <- tapply(newdat8$cpue8, newdat8$llxx, mean)*3000
    regx8 <- substring(names(llmn8),1,1)
    latlongx8 <- substring(names(llmn8), 3)
    areax8 <- cell_areas[match(latlongx8, cell_areas$latlong), "areax"]
    areax8[is.na(areax8)] <- 0
    wts8 <- tapply(llmn8 * areax8, regx8, sum)
    wts8 <- wts8 / max(wts8)
    assign(x = paste0(sp, "_wts8_",pdnm[pd]), value = wts8)

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
    plot_patterns(llmn8, sp, spreg)
    savePlot(paste0("relative wt8_",pdnm[pd],sp),type = "png")
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts1")]] <- mwts/max(mwts)
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts2")]] <- wts2
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts3")]] <- wts3
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts4")]] <- wts4
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts5")]] <- wts5
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts6")]] <- wts6
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts7")]] <- wts7
    allwts[[paste0(spreg,"_", pdnm[pd],"_wts8")]] <- wts8

    windows(height= 10, width = 8); par(mfrow = c(2,1), oma = c(0,0,2,0), mar = c(5,3,1,2)+.1)
    a <- with(model6, tapply(residuals, list(data$yq, data$reg), mean))
    yq <- as.numeric(rownames(a))
    plot(yq, yq, type = "n", xlab = "Year-quarter", ylim = c(-1.5,1.8))
    regseq <- sort(unique(ddd$reg))
    for(r in 1:length(regseq)) {
      lines(yq, a[,r], col = r)
    }
    legend("top", legend = regseq, lty = 1, col = 1:length(regseq), horiz = TRUE)
    boxplot(model6$residuals ~ model6$data$Year, xlab = "Year")
    title(main = paste(sp, " ", pds[pd,"st"],"-", pds[pd,"nd"]), outer = "TRUE")
    savePlot(paste0("resids_by_reg_",pdnm[pd],sp),type = "png")

    windows(height=15,width=8); par(mfrow=c(4,2), mar = c(4,4,3,1))
    plotdiags(model23$residuals, ti = paste(sp,"m2"))
    plotdiags(model4$residuals, ti = paste(sp,"m4"))
    plotdiags(model5$residuals, ti = paste(sp,"m5"))
    plotdiags(model6$residuals, ti = paste(sp,"m6"))
    savePlot(paste0("resid_dbns_",pdnm[pd],sp),type = "png")

    drop6 <- drop1(model6, test = "LRT")
    write.csv(drop6, paste0(sp,"drop6",pdnm[pd],sp,".csv"))
    drop5 <- drop1(model5)
    write.csv(drop5, paste0(sp,"drop5",pdnm[pd],sp,".csv"))

    graphics.off()
  }
  # save(BET_wts, BET_wts_mean, file = paste0("BET_wts_",pdnm[pd],".RData"))
  # save(YFT_wts, YFT_wts_mean, file = paste0("YFT_wts_",pdnm[pd],".RData"))
}
save(allwts, file = paste0("allwts.RData"))
load(file = "allwts.RData")


cbind(areax, regx)
barplot(tapply(areax, regx, mean))
barplot(tapply(llmn23[areax==0], regx[areax==0], mean))
barplot(rbind(mwts,wts2,wts3,wts4,wts5,wts6), beside = T)

posy <- grep("YFT2",names(allwts))
posb <- grep("BET3",names(allwts))
allwtsy <- t(as.matrix(data.frame(allwts[posy])))
allwtsb <- t(as.matrix(data.frame(allwts[posb])))
write.csv(allwtsy, "allwt_yft.csv")
write.csv(allwtsb, "allwt_bet.csv")

a5 <- rbind(cbind(tp=5, sp="BET", read.csv(paste0("BETdrop57994BET.csv"))),
            cbind(tp=5, sp="YFT", read.csv(paste0("YFTdrop57994YFT.csv"))),
            cbind(tp=6, sp="BET", read.csv(paste0("BETdrop67994BET.csv"))[,1:4]),
            cbind(tp=6, sp="YFT", read.csv(paste0("YFTdrop67994YFT.csv"))[,1:4]))
a <- a5 %>% group_by(sp, tp) %>%
      mutate(deltaAIC = AIC - min(AIC))
write.csv(a, file = "table1.csv")


hist(summary(model6)$coefficients[,4])
# predict relative catch rates for the same kind of effort for all quarters in 1990-2000
# average the prediction across vessels. Choose an HBF or average across. Choose a time period. Sum across areas.
str(lldat)
table(lldat$regB3)
reglist <- data.frame(index=1:5,regnumY=1:5, regnumB=c(1:5), regY2=c("2S",3,4,5,"2N"),regB3=c("1S",2,3,4,"1N"))


legtxt.pd <- c("1960-1975","1963-1975","1975-1994","1979-1994","1980-2000")
windows()
sp="YFT"
for (mdi in 1:8) {
  for(spr in c("YFT2")) {
    nmx <- paste0(paste(spr, pdnm, "wts", sep = "_"), mdi)
    a <- data.frame()
    for(i in 1:5) a <- rbind(a, with(allwts, get(nmx[i])))
    colnames(a) <- 1:7
    a <- a[,c(1,7,2,3,4,5,6)]
    colnames(a) <- paste0("R", c(1,"2N", "2S", 3,4,5,6))
    barplot(as.matrix(a), beside=TRUE, names.arg = colnames(a), legend = legtxt.pd, args.legend = list(x="topleft", adj = c(0,0.1)), main = paste("Model", mdi))
    savePlot(paste("barplot", sp, "wts_md", mdi, sep = "_"), type = "png")
  }
}

legtxt.pd <- c("1960-1975","1963-1975","1975-1994","1979-1994","1980-2000")
windows()
sp="BET"
for (mdi in 1:8) {
  for(spr in c("BET3")) {
    nmx <- paste0(paste(spr, pdnm, "wts", sep = "_"), mdi)
    a <- data.frame()
    for(i in 1:5) a <- rbind(a, with(allwts, get(nmx[i])))
    colnames(a) <- 1:5
    a <- a[,c(5,1,2,3,4)]
    sth <- apply(a[,4:5],1,sum)
    a <- cbind(a[,1:3], sth)
    a <- a / apply(a,1,max)
    colnames(a) <- paste0("R", c("1N", "1S", 2,3))
    barplot(as.matrix(a), beside=TRUE, names.arg = colnames(a), legend = legtxt.pd, args.legend = list(x="topleft", adj = c(0,0.1)), main = paste("Model", mdi))
    savePlot(paste("barplot", sp, "wts_md", mdi, sep = "_"), type = "png")
  }
}

# windows()
# a<- with(allwts, rbind(YFT2_6075_wts6, YFT2_6375_wts6, YFT2_7594_wts6, YFT2_7994_wts6, YFT2_8000_wts6))
# a <- a[,c(1,7,2,3,4,5,6)]
# colnames(a) <- paste0("R", c(1,"2N", "2S", 3,4,5,6))
# barplot(a, beside=TRUE, names.arg = colnames(a), legend = c("1960-1975","1963-1975","1975-1994","1979-1994","1980-2000"), args.legend = list(x="topleft"))
# savePlot("barplot_YFT6_pds", type = "png")
#
# windows()
# a<- with(allwts, rbind(BET3_6075_wts6, BET3_6375_wts6, BET3_7594_wts6, BET3_7994_wts6, BET3_8000_wts6))
# a <- a[,c(5,1,2,3,4)]
# sth <- apply(a[,4:5],1,sum)
# a <- cbind(a[,1:3], sth)
# a <- a / apply(a,1,max)
# colnames(a) <- paste0("R", c("1N", "1S", 2,3))
# barplot(a, beside=TRUE, names.arg = colnames(a), legend = c("1960-1975","1963-1975","1975-1994","1979-1994","1980-2000"), args.legend = list(x="topleft"))
# savePlot("barplot_BET6_pds", type = "png")

legtxt <- c("m1 means","m2 stdized","m3 areas","m4 stat m","m5 fleet","m6 qtr","m7 gam", "m8 merged")

windows()
sp="YFT"
for (pd in 1:5) {
  pdx <- pdnm[pd]
  for(spr in c("YFT2")) {
    nmx <- paste0(paste(spr, pdx, "wts", sep = "_"), 1:8)
    a <- data.frame()
    for(i in 1:8) a <- rbind(a, with(allwts, get(nmx[i])))
    colnames(a) <- 1:7
    a <- a[,c(1,7,2,3,4,5,6)]
    colnames(a) <- paste0("R", c(1,"2N", "2S", 3,4,5,6))
    barplot(as.matrix(a), beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft", y.intersp=0.8, adj=c(0,0.1)), main = paste(pds[pd, "st"], "-", pds[pd, "nd"]))
    savePlot(paste("barplot", sp, "wts", pdx, sep = "_"), type = "png")
  }
}

windows()
sp="BET"
for (pd in 1:5) {
  pdx <- pdnm[pd]
  for(spr in c("BET3")) {
    nmx <- paste0(paste(spr, pdx, "wts", sep = "_"), 1:8)
    a <- data.frame()
    for(i in 1:8) a <- rbind(a, with(allwts, get(nmx[i])))
    colnames(a) <- 1:5
    a <- a[,c(5,1,2,3,4)]
    sth <- apply(a[,4:5],1,sum)
    a <- cbind(a[,1:3], sth)
    a <- a / apply(a,1,max)
    colnames(a) <- paste0("R", c("1N", "1S", 2,3))
    barplot(as.matrix(a), beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft", y.intersp=0.8, adj=c(0,0.1)), main = paste(pds[pd, "st"], "-", pds[pd, "nd"]))
    savePlot(paste("barplot", sp, "wts", pdx, sep = "_"), type = "png")
  }
}


# Use the code below to prepare rescaled indices
# The functions ending in y work with annual indices, otherwise yq indices

# combines 2 indices
dowts <- function(d1, d2, wts, yq) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yq %in% yq], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yq %in% yq], na.rm = TRUE)
  dd <- data.frame(yq=d1$yq,pr=d1$pr * wts[1] + d2$pr[match(d1$yq,d2$yq)] * wts[2])
  dd$pr <- dd$pr/mean(dd$pr,na.rm=TRUE)
  plot(d1$yq,d1$pr,type="l",col=1, xlab = "Year-quarter", ylab = "Relative CPUE")
  lines(d2$yq,d2$pr,col=2)
  lines(dd$yq,dd$pr,col=3)
  return(dd)
}

# combines 3 indices
dowts_3 <- function(d1, d2, d3, wts, yq) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yq %in% yq], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yq %in% yq], na.rm = TRUE)
  d3$pr <- d3$pr / mean(d3$pr[d3$yq %in% yq], na.rm = TRUE)
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

# combines 4 indices
dowts_4 <- function(d1, d2, d3, d4, wts, yq) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yq %in% yq], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yq %in% yq], na.rm = TRUE)
  d3$pr <- d3$pr / mean(d3$pr[d3$yq %in% yq], na.rm = TRUE)
  d4$pr <- d4$pr / mean(d4$pr[d4$yq %in% yq], na.rm = TRUE)
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

# combines 5 indices
dowts_5 <- function(d1, d2, d3, d4, d5, wts, yq) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yq %in% yq], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yq %in% yq], na.rm = TRUE)
  d3$pr <- d3$pr / mean(d3$pr[d3$yq %in% yq], na.rm = TRUE)
  d4$pr <- d4$pr / mean(d4$pr[d4$yq %in% yq], na.rm = TRUE)
  d5$pr <- d5$pr / mean(d5$pr[d5$yq %in% yq], na.rm = TRUE)
  dd <- data.frame(yq=d1$yq,pr=d1$pr * wts[1] + d2$pr[match(d1$yq,d2$yq)] * wts[2] + d3$pr[match(d1$yq,d3$yq)] * wts[3] + d4$pr[match(d1$yq,d4$yq)] * wts[4] + d5$pr[match(d1$yq,d5$yq)] * wts[5])
  dd$pr <- dd$pr/mean(dd$pr,na.rm=TRUE)
  plot(d1$yq,d1$pr,type="l",col=1, xlab = "Year-quarter", ylab = "Relative CPUE")
  lines(d2$yq,d2$pr,col=2)
  lines(d3$yq,d3$pr,col=3)
  lines(d4$yq,d4$pr,col=4)
  lines(d5$yq,d5$pr,col=5)
  lines(dd$yq,dd$pr,col=6, lwd = 2)
  return(dd)
}

dowts_5y <- function(d1, d2, d3, d4, d5, wts, yrs) {
  d1$pr <- d1$pr / mean(d1$pr[d1$yr %in% yrs], na.rm = TRUE)
  d2$pr <- d2$pr / mean(d2$pr[d2$yr %in% yrs], na.rm = TRUE)
  d3$pr <- d3$pr / mean(d3$pr[d3$yr %in% yrs], na.rm = TRUE)
  d4$pr <- d4$pr / mean(d4$pr[d4$yr %in% yrs], na.rm = TRUE)
  d5$pr <- d5$pr / mean(d5$pr[d5$yr %in% yrs], na.rm = TRUE)
  dd <- data.frame(yr=d1$yr,pr=d1$pr * wts[1] + d2$pr[match(d1$yr,d2$yr)] * wts[2] + d3$pr[match(d1$yr,d3$yr)] * wts[3] + d4$pr[match(d1$yr,d4$yr)] * wts[4] + d5$pr[match(d1$yr,d5$yr)] * wts[5])
  dd$pr <- dd$pr/mean(dd$pr,na.rm=TRUE)
  plot(d1$yr,d1$pr,type="l",col=1, xlab = "Year-quarter", ylab = "Relative CPUE")
  lines(d2$yr,d2$pr,col=2)
  lines(d3$yr,d3$pr,col=3)
  lines(d4$yr,d4$pr,col=4)
  lines(d5$yr,d5$pr,col=5)
  lines(dd$yr,dd$pr,col=6, lwd = 2)
  return(dd)
}

# 
# Load bigeye indices 2017
yq <- seq(1980.125, 1999.875, .25)
yq6075 <- seq(1963.125, 1975.875, .25)
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


# get scaling factors, combine indices, and save the results
wts2 <- allwts$BET3_8000_wts6[c(5,1,2) + 1]
wts2_early <- allwts$BET3_6375_wts6[c(5,1,2) + 1]
windows(12,12); par(mfrow = c(2,2))
dd <- dowts_3(bet1n, bet1s, bet2, wts2, yq)
write.csv(dd,"Joint_regB3_R12_dellog_boat_allyrs.csv")
dd <- dowts_3(betnovess1n, betnovess1s, betnovess2, wts2, yq)
write.csv(dd,"Joint_regB3_R12_dellog_novess_allyrs.csv")
dd <- dowts_3(bet5279_1n, bet5279_1s, bet5279_2, wts2_early, yq = yq6075)
write.csv(dd,"Joint_regB3_R12_dellog_novess_5279.csv")
dd <- dowts_3(bet79nd_1n, bet79nd_1s, bet79nd_2, wts2, yq)
write.csv(dd,"Joint_regB3_R12_dellog_vessid_79nd.csv")
savePlot("Joint_regB3_R12_dellog", type = "png")

# Annual bigeye indices 2017
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


# Quarterly yellowfin indices 2018
yq <- seq(1979.125, 1993.875, .25)
yq6075 <- seq(1963.125, 1974.875, .25)
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
wts8 <- allwts$YFT2_7594_wts8[c(7,2,3,4,5)]
wts8_early <- allwts$YFT2_6375_wts8[c(7,2,3,4,5)]
windows(12,12); par(mfrow = c(2,2))
dd <- dowts_5(yft2n, yft2s, yft3, yft4, yft5, wts8, yq)
write.csv(dd,"Joint_regY2_R2345_dellog_boat_allyrs_yq.csv")
dd <- dowts_5(yftnovess2n, yftnovess2s, yftnovess3, yftnovess4, yftnovess5, wts8, yq)
write.csv(dd,"Joint_regY2_R2345_dellog_novess_allyrs_yq.csv")
dd <- dowts_5(yft5279_2n, yft5279_2s, yft5279_3, yft5279_4, yft5279_5, wts8_early, yq = yq6075)
write.csv(dd,"Joint_regY2_R2345_dellog_novess_5279_yq.csv")
dd <- dowts_5(yft79nd_2n, yft79nd_2s, yft79nd_3, yft79nd_4, yft79nd_5, wts8, yq)
write.csv(dd,"Joint_regY2_R2345_dellog_vessid_79nd_yq.csv")
savePlot("Joint_regY2_R2345_dellog", type = "png")

# Annual yellowfin indices 2018
yrs <- seq(1975, 1993, 1)
yrs6375 <- seq(1963, 1974, 1)
yft5y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_boat_allyrs_yr.csv"))
yft2y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_boat_allyrs_yr.csv"))
yft3y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_boat_allyrs_yr.csv"))
yft4y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R4_dellog_boat_allyrs_yr.csv"))
yft2ny <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_boat_allyrs_yr.csv"))
yft2sy <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_boat_allyrs_yr.csv"))
yftnovess5y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_novess_allyrs_yr.csv"))
yftnovess2y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_novess_allyrs_yr.csv"))
yftnovess2ny <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_novess_allyrs_yr.csv"))
yftnovess2sy <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_novess_allyrs_yr.csv"))
yftnovess3y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_novess_allyrs_yr.csv"))
yftnovess4y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R4_dellog_novess_allyrs_yr.csv"))
yft5279_5y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_novess_5279_yr.csv"))
yft5279_2y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_novess_5279_yr.csv"))
yft5279_2ny <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_novess_5279_yr.csv"))
yft5279_2sy <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_novess_5279_yr.csv"))
yft5279_3y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_novess_5279_yr.csv"))
yft5279_4y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R4_dellog_novess_5279_yr.csv"))
yft79nd_5y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R5_dellog_vessid_79nd_yr.csv"))
yft79nd_2y <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY_R2_dellog_vessid_79nd_yr.csv"))
yft79nd_2ny <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R7_dellog_vessid_79nd_yr.csv"))
yft79nd_2sy <- read.csv(paste0(jntalysis_dir,"cl0_hb1_hk1/outputs/Joint_regY2_R2_dellog_vessid_79nd_yr.csv"))
yft79nd_3y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R3_dellog_vessid_79nd_yr.csv"))
yft79nd_4y <- read.csv(paste0(jntalysis_dir,"cl1_hb0_hk1/outputs/Joint_regY_R4_dellog_vessid_79nd_yr.csv"))
wts8 <- allwts$YFT2_7594_wts8[c(7,2,3,4,5)]
wts8_early <- allwts$YFT2_6375_wts8[c(7,2,3,4,5)]
windows(12,12); par(mfrow = c(2,2))
dd <- dowts_5y(yft2ny, yft2sy, yft3y, yft4y, yft5y, wts8, yrs)
write.csv(dd,"Joint_regY2_R2345_dellog_boat_allyrs_yr.csv")
dd <- dowts_5y(yftnovess2ny, yftnovess2sy, yftnovess3y, yftnovess4y, yftnovess5y, wts8, yrs)
write.csv(dd,"Joint_regY2_R2345_dellog_novess_allyrs_yr.csv")
dd <- dowts_5y(yft5279_2ny, yft5279_2sy, yft5279_3y, yft5279_4y, yft5279_5y, wts8_early, yrs = yrs6375)
write.csv(dd,"Joint_regY2_R2345_dellog_novess_5279_yr.csv")
dd <- dowts_5y(yft79nd_2ny, yft79nd_2sy, yft79nd_3y, yft79nd_4y, yft79nd_5y, wts8, yrs)
write.csv(dd,"Joint_regY2_R2345_dellog_vessid_79nd_yr.csv")
savePlot("Joint_regY2_R2345_dellogyr", type = "png")

##############
# Combine indices - an alternative approach for YFT
yq6375 <- seq(1963.125, 1974.875, .25)
yq6075 <- seq(1960.125, 1974.875, .25)
yq7594 <- seq(1975.125, 1993.875, .25)
yq7994 <- seq(1979.125, 1993.875, .25)
yq8000 <- seq(1980.125, 1999.875, .25)
yr6075 <- 1960:1974
yr6375 <- 1963:1974
yr7594 <- 1975:1993
yr7994 <- 1979:1993
yr8000 <- 1980:1999

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


yftlist <- list(yft2n, yft2s, yft3, yft4, yft5)
yftlist_ord <- c(7,2,3,4,5)
length(yftlist)

dowts_list <- function(dlist, wts, yq) { # Note that the regions in dlist & wts must be in the same sequence
  for(i in 1:length(dlist)) {
    a <- dlist[[i]]
    a$prwt <- wts[i] * a$pr / mean(a$pr[a$yq %in% yq], na.rm = TRUE)
    dlist[[i]] <- a
  }
  return(dlist)
}
load(file = "allwts.RData")

# Calulate means for each species and period
rng_mn <- function(x,y) mean(x$pr[x$yq %in% y], na.rm = TRUE)

yrng6075 <- c(rng_mn(yft2n, yq6075),rng_mn(yft2s, yq6075),rng_mn(yft3, yq6075),rng_mn(yft4, yq6075),rng_mn(yft5, yq6075))
yrng6375 <- c(rng_mn(yft2n, yq6375),rng_mn(yft2s, yq6375),rng_mn(yft3, yq6375),rng_mn(yft4, yq6375),rng_mn(yft5, yq6375))
yrng7594 <- c(rng_mn(yft2n, yq7594),rng_mn(yft2s, yq7594),rng_mn(yft3, yq7594),rng_mn(yft4, yq7594),rng_mn(yft5, yq7594))
yrng7994 <- c(rng_mn(yft2n, yq7994),rng_mn(yft2s, yq7994),rng_mn(yft3, yq7994),rng_mn(yft4, yq7994),rng_mn(yft5, yq7994))
yrng8000 <- c(rng_mn(yft2n, yq8000),rng_mn(yft2s, yq8000),rng_mn(yft3, yq8000),rng_mn(yft4, yq8000),rng_mn(yft5, yq8000))


yft6375st <- dowts_list(dlist = yftlist, wts = allwts$YFT2_6375_wts8[yftlist_ord], yq = yq6075)
yft6075st <- dowts_list(dlist = yftlist, wts = allwts$YFT2_6075_wts8[yftlist_ord], yq = yq6075)
yft7594st <- dowts_list(dlist = yftlist, wts = allwts$YFT2_7594_wts8[yftlist_ord], yq = yq7594)
yft7994st <- dowts_list(dlist = yftlist, wts = allwts$YFT2_7994_wts8[yftlist_ord], yq = yq7994)
yft8000st <- dowts_list(dlist = yftlist, wts = allwts$YFT2_8000_wts8[yftlist_ord], yq = yq8000)

# Barplots comparing each wt between periods
mkwt <- function(wtobj, wtnum) paste0(wtobj, "_wts", wtnum)

for (wt in 1:6) {
  a<- with(allwts, rbind(get(mkwt("YFT2_6075", wt)),get(mkwt("YFT2_6375", wt)),get(mkwt("YFT2_7594", wt)),get(mkwt("YFT2_7994", wt)),get(mkwt("YFT2_8000", wt))))
  a <- a[, yftlist_ord]
  a[1,] <- a[1,] / yrng6075
  a[2,] <- a[2,] / yrng6375
  a[3,] <- a[3,] / yrng7594
  a[4,] <- a[4,] / yrng7994
  a[5,] <- a[5,] / yrng8000
  a <- a / apply(a,1,max)
  colnames(a) <- paste0("R", c("2N", "2S", 3,4,5))
  windows()
  barplot(a, beside = TRUE, names.arg = colnames(a), legend = c("1960-1975","1963-1975","1975-1994","1979-1994","1980-2000"), args.legend = list(x = "topleft"), main = paste0("YFT m", wt))
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
  barplot(a, beside = TRUE, names.arg = colnames(a), legend = c("1960-1975","1963-1975","1975-1994","1979-1994","1980-2000"), args.legend = list(x = "topleft"), main = paste0("BET m", wt))
  savePlot(paste0("barplot_BETadj_pds_wts", wt), type = "png")
}


windows(12,10); par(mfrow = c(3,2))
for(i in 1:5) {
  with(bet8000st[[i]], plot(yq, prwt, xlab = "Year-quarter", ylab = "Weighted CPUE", col= 1, type = "l", main = betlab[i]))
  with(bet6375st[[i]], lines(yq, prwt, col = 2, lty = 2))
}


dd <- dowts(bet1, bet2, BET_wts[-1], yq)
write.csv(dd,"Joint_regB2_R12_dellog_boat_allyrs.csv")
dd <- dowts(betnovess1, betnovess2, BET_wts[-1], yq)
write.csv(dd,"Joint_regB2_R12_dellog_novess_allyrs.csv")
dd <- dowts(bet5279_1, bet5279_2, BET_wts[-1], yq = yqfake)
write.csv(dd,"Joint_regB2_R12_dellog_novess_5279.csv")
dd <- dowts(bet79nd_1, bet79nd_2, BET_wts[-1], yq)
write.csv(dd,"Joint_regB2_R12_dellog_vessid_79nd.csv")
savePlot("Joint_regB2_R12_dellog", type = "png")

windows()
a <- with(allwts, rbind(YFT2_6375_wts1, YFT2_6375_wts2, YFT2_6375_wts3, YFT2_6375_wts4, YFT2_6375_wts5, YFT2_6375_wts8))
a <- a[,c(NA,7,3,4,5,6)]
a[,-1] <- t(t(a[,-1]) / yrng6375)
a <- a / apply(a,1,max, na.rm = TRUE)
colnames(a) <- c("", paste0("R", c("2N", "2S", 3,4,5)))
barplot(a, beside=TRUE, names.arg = colnames(a), legend.text = legtxt, args.legend = list(x="topleft"), main = "YFT 1963 - 1975")
savePlot("barplotadj_YFT_wts_6375", type = "png")

windows()
a<- with(allwts, rbind(YFT2_6075_wts1, YFT2_6075_wts2, YFT2_6075_wts3, YFT2_6075_wts4, YFT2_6075_wts5, YFT2_6075_wts8))
a <- a[,c(NA,8,3,4,5,6)]
a[,-1] <- t(t(a[,-1]) / yrng6075)
a <- a / apply(a,1,max, na.rm = TRUE)
colnames(a) <- c("", paste0("R", c("2N", "2S", 3,4,5)))
barplot(a, beside=TRUE, names.arg = colnames(a), legend.text = legtxt, args.legend = list(x="topleft"), main = "YFT 1960 - 1975")
savePlot("barplotadj_YFT_wts_6075", type = "png")

windows()
a<- with(allwts, rbind(YFT2_8000_wts1, YFT2_8000_wts2, YFT2_8000_wts3, YFT2_8000_wts4, YFT2_8000_wts5, YFT2_8000_wts8))
a <- a[,c(NA,8,3,4,5,6)]
a[,-1] <- t(t(a[,-1]) / yrng8000)
a <- a / apply(a,1,max, na.rm = TRUE)
colnames(a) <- c("", paste0("R", c("2N", "2S", 3,4,5)))
barplot(a, beside=TRUE, names.arg = colnames(a), legend.text = legtxt, args.legend = list(x="topleft"), main = "YFT 1980 - 2000")
savePlot("barplotadj_YFT_wts_8000", type = "png")

windows()
a<- with(allwts, rbind(BET3_6075_wts1,BET3_6075_wts2,BET3_6075_wts3,BET3_6075_wts4,BET3_6075_wts5,BET3_6075_wts8))
a <- a[,c(6,2,3,4,5)-1]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- t(t(a) / brng6075)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft"), main = "BET 1960 - 1975")
savePlot("barplotadj_BET_means_6075", type = "png")

windows()
a<- with(allwts, rbind(BET3_6375_wts1,BET3_6375_wts2,BET3_6375_wts3,BET3_6375_wts4,BET3_6375_wts5,BET3_6375_wts8))
a <- a[,c(6,2,3,4,5)]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- t(t(a) / brng6375)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft"), main = "BET 1963 - 1975")
savePlot("barplotadj_BET_means_6375", type = "png")

windows()
a<- with(allwts, rbind(BET3_8000_wts1,BET3_8000_wts2,BET3_8000_wts3,BET3_8000_wts4,BET3_8000_wts5,BET3_8000_wts8))
a <- a[,c(6,2,3,4,5)]
sth <- apply(a[,4:5],1,sum)
a <- cbind(a[,1:3], sth)
a <- t(t(a) / brng8000)
a <- a / apply(a,1,max)
colnames(a) <- paste0("R", c("1N", "1S", 2,3))
barplot(a, beside=TRUE, names.arg = colnames(a), legend = legtxt, args.legend = list(x="topleft"), main = "BET 1980 - 2000")
savePlot("barplotadj_BET_means_8000", type = "png")
