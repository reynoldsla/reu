#This code was written by Lyra Reynolds for a Research Experience for Undergrads funded by the National Science Foundation
#The packages and base code for these models were written by Dr. Rob Payn at Montana State University.
#The git hub access is below and is publicly available from https://github.com/robpayn. The metabc package
#runs the code through C rather than straight through R studio program. The program will print the output graphs as PDF.
#The code preforms a sensitivity analysis on alkalinity values for the Upper Clarkfork River in Montana.

#Alkalinity SA loop, using best fit values for Fix and Resp ratios (-0.32 and 1.07)
rm(list=ls())

library(metabc)
library(dictools)

setwd("C:/Users/larey/OneDrive - Montana State University/REU NOTES/RUE")


load(file = "./2014-09-12/signal.RData")

#variable name and values for the vector
alkalinity <- c(2850, 2852, 2854, 2856, 2857, 2858, 2860, 2862, 2864, 2866, 2868, 2870)

#line colors for each run in the loop
modelcolors <- c(
  "aquamarine4",
  "aquamarine3",
  "aquamarine2",
  "aquamarine1",
  "black",
  "darkslategray1",
  "darkslategray2",
  "deepskyblue",
  "deepskyblue1",
  "deepskyblue2",
  "deepskyblue3",
  "deepskyblue4"
)

#changing line widths
linewidths <- c(1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1)

#loop simulation output
modelsims <- vector(mode = "list", length = length(alkalinity))

#startloop
for(index in 1:length(alkalinity)){
  # sets up DO baseline values with OneStationMetabDO using NLS
  modelDO <- CMetabDo$new(
    type = "ForwardEuler",
    dailyGPP = 300,
    ratioDoCFix = 1,
    dailyER = 300,
    ratioDoCResp = -1,
    k600 = 18,
    airPressure = 0.84,
    par = signalOut$getVariable("par"),
    initialDO = signalOut$getVariable("do")[1],
    time = signalOut$time,
    temp = signalOut$getVariable("temp"),
    stdAirPressure = 1,
    parTotal = -1
  );
  
  y<- signalOut$getVariable("do")
  #Setting up NLS params for onestation
  runmodelDO<-function(
    params,
    metabModel) 
  {
    if(!is.na(params["ratioDicFix"])) {
      metabModel$setMetabDoParam("RatioDicCFix", params["ratioDicCix"])
    }
    if(!is.na(params["dailyER"])) {
      metabModel$setMetabParam("DailyER", params ["dailyER"])
    }
    if(!is.na(params["k600"])) {
      metabModel$setMetabParam("k600", params ["k600"])
    }
    if(!is.na(params["dailyGPP"])) {
      metabModel$setMetabParam("DailyGPP", params["dailyGPP"])
    }
    
    if(!is.na(params["ratioDoCResp"])) {
      metabModel$setMetabDoParam("RatioDoCResp", params["ratioDoCResp"])
    }
    
    if(!is.na(params["alkalinity"])) {
        metabModel$setMetabDoParam("Alkalinity", params["alkalinity"])
    }
    
    metabModel$run()
    
    return(modelDO$output$do$dox)
    
  }
  
  nlsResults<- nls(
    formula = y ~ runmodelDO(metabModel = modelDO, params = c(dailyGPP = dailyGPPest, dailyER = dailyERest, k600 = k600est)),
    start= c(dailyGPPest = 300, dailyERest = 300, k600est = 15)
  )
  
  #saving coeffcients from OneStation NLS
  coefficients<- summary(nlsResults)$coefficients
  
  runmodelDO(
    metabModel = modelDO,
    params = c(
      dailyGPP = coefficients["dailyGPPest", "Estimate"],
      dailyER = coefficients["dailyERest", "Estimate"],
      k600 = coefficients["k600est", "Estimate"]
    )
  )

  upstreamDIC <- 1e6 * 
    apply(
      X = signalIn$dataFrame$data[,c("temp", "pCO2")], 
      MARGIN = 1,
      FUN = function(row)
      {
        if (is.finite(row["temp"]) && is.finite(row["pCO2"])) {
          ce <- CarbonateEq$new(tempC = row["temp"]);
          return(ce$optDICFromfCO2TotalAlk(
            fCO2 = row["pCO2"], 
            totalAlk = alkalinity[index] * 1e-6
          )$concDIC);
        } else {
          return(NA);
        }
      }
    );
  
  #twostation do dic
  modelDoDic <- CMetabLagrangeDoDic$new(
    type = "CNOneStep",
    dailyGPP = coefficients["dailyGPPest", "Estimate"],
    dailyER = coefficients["dailyERest", "Estimate"],
    k600 = coefficients["k600est", "Estimate"],
    airPressure = 0.84,
    downstreamPAR = signalOut$getVariable("par"),
    upstreamPAR = signalIn$getVariable("par"),
    upstreamDO = signalIn$getVariable("do"),
    downstreamTime = signalOut$time,
    upstreamTime = signalIn$time,
    downstreamTemp = signalOut$getVariable("temp"),
    upstreamTemp = signalIn$getVariable("temp"),
    stdAirPressure = 1,
    parTotal = -1,
    ratioDicCFix = -0.32,
    ratioDicCResp = 1.07,
    upstreamDIC = upstreamDIC,
    pCO2air = 400,
    upstreamAlkalinity = alkalinity[index],
    downstreamAlkalinity = alkalinity[index],
    timesteps = 0
  );
  
  modelDoDic$run()
  
  modelsims[[index]] <- modelDoDic
  modelsims[[index]]$.__enclos_env__$upstreamDIC <- upstreamDIC
  
} #end of SA loop  

#finding min and max of each run of model DoDIC
minmaxmatrix <- sapply(
  X = modelsims,
  FUN = function(modelDoDic)
  {return(
    c(
      minDO = min(modelDoDic$output$do$dox, na.rm = TRUE),
      maxDO = max(modelDoDic$output$do$dox, na.rm = TRUE),
      minDOSat = min(modelDoDic$output$do$downstreamDoSat, na.rm = TRUE),
      maxDOSat = max(modelDoDic$output$do$downstreamDoSat, na.rm = TRUE),
      minpCO2 = min(modelDoDic$output$dic$pCO2, na.rm = TRUE),
      maxpCO2 = max(modelDoDic$output$dic$pCO2, na.rm = TRUE),
      minDIC = min(modelDoDic$output$dic$dic, na.rm = TRUE),
      maxDIC = max(modelDoDic$output$dic$dic, na.rm = TRUE),
      minDeltaDO = min(modelDoDic$output$do$dox - signalIn$getVariable("do"), na.rm = TRUE),
      maxDeltaDO = max(modelDoDic$output$do$dox - signalIn$getVariable("do"), na.rm = TRUE),
      minDeltaDic = min(modelDoDic$output$dic$dic - modelDoDic$.__enclos_env__$upstreamDIC, na.rm = TRUE),
      maxDeltaDic = max(modelDoDic$output$dic$dic - modelDoDic$.__enclos_env__$upstreamDIC, na.rm = TRUE),
      minDeltapCO2 = min(modelDoDic$output$dic$pCO2 - signalIn$getVariable("pCO2"), na.rm = TRUE),
      maxDeltapCO2 = max(modelDoDic$output$dic$pCO2 - signalIn$getVariable("pCO2"), na.rm = TRUE),
      minPH = min(modelDoDic$output$dic$downstreampH, na.rm = TRUE),
      maxPH = max(modelDoDic$output$dic$downstreampH, na.rm = TRUE)
    )
  )}
)

#sending to PDF
pdf(
  file = "./AltAlkalinitySA.pdf",
  height = 8,
  width = 12
)
par(
    mar = c(4, 4.5, 1, 4.5),
    layout(
      mat = matrix(
        data = c(1,2,3,4,5,5),
        nrow = 2,
        ncol = 3
      ),
      widths = c(1.1, 1, 1)
    )
)


#limits for new plot
ymin <- min(minmaxmatrix["minDIC",], signalOut$getVariable("dic"), na.rm = TRUE)
ymax <- max(minmaxmatrix["maxDIC",], signalOut$getVariable("dic"), na.rm = TRUE)
yminMass <- (ymin * 0.012)
ymaxMass <- (ymax * 0.012)
yminPAR <- min(signalOut$getVariable("par")) 
ymaxPAR <- max(signalOut$getVariable("par"))

#setting up PAR values behind lines and points
plot.new()
plot.window(
  xlim = c(min(signalOut$time),
           max(signalOut$time)),
  ylim = c(ymaxPAR, yminPAR + 0.04*(ymaxPAR-yminPAR)) #removing buffer to put PAR on edge
)
#PAR shape
polygon(
  x = signalOut$time,
  y = signalOut$getVariable("par"),
  lty = "blank",
  col = "lightgoldenrod1"
)

par(
  new = TRUE
)

#plot DIC concentrations
plot(
  x = modelDoDic$downstreamTimePOSIX,
  y = modelDoDic$output$dic$dic,
  type = "n",
  ylab = bquote(.("DIC-C Conc (") * mu * mol ~ L^-1 * .(")")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(ymin, ymax)
)

points(
  x = signalOut$time,
  y = signalOut$getVariable("dic")
)

#create legend
legend(
  x = "bottomleft",
  bty = "n",
  legend = c(alkalinity, "PAR Dist."),
  # legend = sprintf("Dic:C Fixation = %.2f", ratioDic),
  lty = c("solid"),
  col = c(modelcolors, "lightgoldenrod1"),
  lwd = c(linewidths, 13),
  title = "Alkalinity = "
)

lapply(
  X = 1:length(alkalinity),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$dic$dic,
      col = modelcolors[index],
      lwd = linewidths[index]
    )
  }
)

#secondary axis
par(
  new = TRUE
)

plot.new()
plot.window(
  xlim = c(0,1),
  ylim = c(yminMass, ymaxMass)
)

axis(
  side = 4
)

mtext(
  text = bquote(.("DIC-C Conc (") * mg ~ L^-1 * .(")")),
  side = 4,
  line = 2.5,
  cex = 0.7
)

ymin <- min(minmaxmatrix["minDeltaDic",], signalOut$getVariable("dic")-signalIn$getVariable("dic"), na.rm = TRUE)
ymax <- max(minmaxmatrix["maxDeltaDic",], signalOut$getVariable("dic")-signalIn$getVariable("dic"), na.rm = TRUE)
yminPAR <- min(signalOut$getVariable("par")) 
ymaxPAR <- max(signalOut$getVariable("par"))

#setting up PAR values behind lines and points
plot.new()
plot.window(
  xlim = c(min(signalOut$time),
           max(signalOut$time)),
  ylim = c(ymaxPAR, yminPAR + 0.04*(ymaxPAR-yminPAR)) #removing buffer to put PAR on edge
)
polygon(
  x = signalOut$time,
  y = signalOut$getVariable("par"),
  lty = "blank",
  col = "lightgoldenrod1"
)

#tells window to stay on same plot rather than moving to new plot
par(
  new = TRUE
)

plot(
  x = modelDoDic$downstreamTimePOSIX,
  y = modelDoDic$output$dic$dic-modelDoDic$.__enclos_env__$upstreamDIC,
  type = "n",
  ylab = bquote(.("Change in DIC-C Conc (") * mu * mol ~ L^-1 * .(")")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(ymin, ymax)
)

points(
  x = signalOut$time,
  y = signalOut$getVariable("dic") - signalIn$getVariable("dic")
)

lapply(
  X = 1:length(alkalinity),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$dic$dic - modelsims[[index]]$.__enclos_env__$upstreamDIC,
      col = modelcolors[index],
      lwd = linewidths[index]
    )
  }
)

#secondary axis
par(
  new = TRUE
)
plot.new()
plot.window(
  xlim = c(0,1),
  ylim = c(yminMass, ymaxMass)
)

axis(
  side = 4
)

mtext(
  text = bquote(.("Change in DIC-C Conc (") * mg ~ L^-1 * .(")")),
  side = 4,
  line = 2.5,
  cex = 0.7
)

ymin <- min(minmaxmatrix["minpCO2",], signalOut$getVariable("pCO2"), na.rm = TRUE)
ymax <- max(minmaxmatrix["maxpCO2",], signalOut$getVariable("pCO2"), na.rm = TRUE)
yminPAR <- min(signalOut$getVariable("par")) 
ymaxPAR <- max(signalOut$getVariable("par"))

#setting up PAR values behind lines and points
mar = par("mar")
mar[4] = 1.2
par(
  mar = mar
)

plot.new()
plot.window(
  xlim = c(min(signalOut$time),
           max(signalOut$time)),
  ylim = c(ymaxPAR, yminPAR + 0.04*(ymaxPAR-yminPAR)) #removing buffer to put PAR on edge
)
polygon(
  x = signalOut$time,
  y = signalOut$getVariable("par"),
  lty = "blank",
  col = "lightgoldenrod1"
)

#tells window to stay on same plot rather than moving to new plot
par(
  new = TRUE
)

plot(
  x = modelDoDic$downstreamTimePOSIX,
  y = modelDoDic$output$dic$pCO2,
  type = "n",
  ylab = bquote(.("pCO2 Conc (") * mu * atm * .(")")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(ymin, ymax)
)

points(
  x = signalOut$time,
  y = signalOut$getVariable("pCO2")
)

lapply(
  X = 1:length(alkalinity),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$dic$pCO2,
      col = modelcolors[index],
      lwd = linewidths[index]
    )
  }
)


#new plot
ymin <- min(minmaxmatrix["minDeltapCO2",], signalOut$getVariable("pCO2")-signalIn$getVariable("pCO2"), na.rm = TRUE)
ymax <- max(minmaxmatrix["maxDeltapCO2",], signalOut$getVariable("pCO2")-signalIn$getVariable("pCO2"), na.rm = TRUE)
yminPAR <- min(signalOut$getVariable("par")) 
ymaxPAR <- max(signalOut$getVariable("par"))

#setting up PAR values behind lines and points
plot.new()
plot.window(
  xlim = c(min(signalOut$time),
           max(signalOut$time)),
  ylim = c(ymaxPAR, yminPAR + 0.04*(ymaxPAR-yminPAR)) #removing buffer to put PAR on edge
)
polygon(
  x = signalOut$time,
  y = signalOut$getVariable("par"),
  lty = "blank",
  col = "lightgoldenrod1"
)

#tells window to stay on same plot rather than moving to new plot
par(
  new = TRUE
)

plot(
  x = modelDoDic$downstreamTimePOSIX,
  y = modelDoDic$output$dic$pCO2-signalIn$getVariable("pCO2"),
  type = "n",
  ylab = bquote(.("Change in pCO2 Conc (") * mu * atm * .(")")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(ymin, ymax)
)

points(
  x = signalOut$time,
  y = signalOut$getVariable("pCO2") - signalIn$getVariable("pCO2")
)

lapply(
  X = 1:length(alkalinity),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$dic$pCO2-signalIn$getVariable("pCO2"),
      col = modelcolors[index],
      lwd = linewidths[index]
    )
  }
)

#plot pH
mar = par("mar")
mar[4] = 1.3
par(
  mar = mar
)

plot.new()
plot.window(
  xlim = c(min(signalOut$time),
           max(signalOut$time)),
  ylim = c(ymaxPAR, yminPAR + 0.04*(ymaxPAR-yminPAR)) #removing buffer to put PAR on edge
)
polygon(
  x = signalOut$time,
  y = signalOut$getVariable("par"),
  lty = "blank",
  col = "lightgoldenrod1"
)

#tells window to stay on same plot rather than moving to new plot
par(
  new = TRUE
)

plot(
  x = modelDoDic$downstreamTimePOSIX,
  y = modelDoDic$output$dic$downstreampH,
  type = "n",
  ylab = bquote(.("pH")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(min(minmaxmatrix["minPH",], na.rm = TRUE),
          max(minmaxmatrix["maxPH",], na.rm = TRUE))
)

lapply(
  X = 1:length(alkalinity),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$dic$downstreampH,
      col = modelcolors[index],
      lwd = linewidths[index]
    )
  }
)

dev.off()