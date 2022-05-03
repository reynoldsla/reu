#This code was written by Lyra Reynolds for a Research Experience for Undergrads funded by the National Science Foundation
#The packages and base code for these models were written by Dr. Rob Payn at Montana State University.
#The git hub access is below and is publicly available from https://github.com/robpayn. The metabc package
#runs the code through C rather than straight through R studio program. The program will print the output graphs as PDF.


# devtools::install_github(
#   repo = "robpayn/metabc",
#   ref = "main",
#   subdir = "pkg"
# )


#Running a loop with varying values for Fixation Stoichiometric Ratio, respiration ratio is set to 1.07

rm(list=ls()) #clearing variables & environment

library(metabc) #loading library

setwd("C:/Users/larey/OneDrive - Montana State University/REU NOTES/RUE") #set working directory

load(file = "./2014-09-12/signal.RData") #stream data

#variable name and values for the vector
ratioDic <- c(-0.2, -0.26, -0.32, -0.38, -0.44, -1)
ratioDo <- -1/ratioDic

#line colors for each run in the loop
modelcolors <- c("limegreen", "blue", "black", "magenta", "red", "cyan3")

#changing line widths
linewidths <- c(1, 1, 2, 1, 1, 1)

#loop simulation output
modelsims <- vector(mode = "list", length = length(ratioDic))

#startloop
for(index in 1:length(ratioDic)){
  # sets up DO baseline values with OneStation Model using NLS
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
  #Setting up NLS parameters for OneStation model
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
  
  #Two Station DO-DIC Model
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
    ratioDicCFix = ratioDic[index],
    ratioDicCResp = 1.07,
    upstreamDIC = signalIn$getVariable("dic"),
    pCO2air = 400,
    upstreamAlkalinity = 2857,
    downstreamAlkalinity = 2857,
    timesteps = 0
  );
  
  modelDoDic$run()
  
  modelsims[[index]] <- modelDoDic
  
} #end of SA loop  

#finding min and max of each run of model DoDIC and min and max of the change over the reach
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
      minDeltaDic = min(modelDoDic$output$dic$dic - signalIn$getVariable("dic"), na.rm = TRUE),
      maxDeltaDic = max(modelDoDic$output$dic$dic - signalIn$getVariable("dic"), na.rm = TRUE),
      minDeltapCO2 = min(modelDoDic$output$dic$pCO2 - signalIn$getVariable("pCO2"), na.rm = TRUE),
      maxDeltapCO2 = max(modelDoDic$output$dic$pCO2 - signalIn$getVariable("pCO2"), na.rm = TRUE)
    )
  )}
)
# 
#Send plots to PDF
pdf(
  file = "./FixationLoop.pdf",
  height = 8,
  width = 12
)

#setting up plotting windows
par(
  mar = c(4, 4.5, 1, 4.5),
  mfcol = c(2, 3)
)

#limit values
ymin <- min(minmaxmatrix["minDO",], signalOut$getVariable("do"), na.rm = TRUE)
ymax <- max(minmaxmatrix["maxDO",], signalOut$getVariable("do"), na.rm = TRUE)
yminMass <- (ymin * 0.032)
ymaxMass <- (ymax * 0.032)
yminPAR <- min(signalOut$getVariable("par")) 
ymaxPAR <- max(signalOut$getVariable("par"))

#setting up PAR values behind lines and points
plot.new()
plot.window(
  xlim = c(min(signalOut$time),
           max(signalOut$time)),
  ylim = c(ymaxPAR, yminPAR + 0.04*(ymaxPAR-yminPAR)) #removing buffer to put PAR on edge
)
#create shape for PAR (shows photoperiod)
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

#set up plot DO from loop
plot(
  x = modelDoDic$downstreamTimePOSIX,
  y = modelDoDic$output$do$dox,
  type = "n",
  ylab = bquote(.("DO Conc (") * mu * mol ~ L^-1 * .(")")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(ymin, ymax)
)

#observed data from file load
points(
  x = signalOut$time,
  y = signalOut$getVariable("do")
)

#plot DO concentrations from the loop
lapply(
  X = 1:length(ratioDic),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$do$dox,
      col = modelcolors[index],
      lwd = linewidths[index]
    )
  }
)

#secondary axis on same plot for different units
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
  text = bquote(.("DO Conc (") * mg ~ L^-1 * .(")")),
  side = 4,
  line = 2.5,
  cex = 0.7
)

#limits for new plot
ymin <- min(minmaxmatrix["minDeltaDO",], signalOut$getVariable("do")-signalIn$getVariable("do"), na.rm = TRUE)
ymax <- max(minmaxmatrix["maxDeltaDO",], signalOut$getVariable("do")-signalIn$getVariable("do"), na.rm = TRUE)
yminPAR <- min(signalOut$getVariable("par")) 
ymaxPAR <- max(signalOut$getVariable("par"))

#setting up PAR values behind lines and points
plot.new()
plot.window(
  xlim = c(min(signalOut$time),
           max(signalOut$time)),
  ylim = c(ymaxPAR, yminPAR + 0.04*(ymaxPAR-yminPAR)) #removing buffer to put PAR on edge
)
#PAR Shape
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

#plot change in DO concentrations
plot(
  x = modelDoDic$downstreamTimePOSIX,
  y = modelDoDic$output$do$dox-signalIn$getVariable("do"),
  type = "n",
  ylab = bquote(.("Change in DO Conc (") * mu * mol ~ L^-1 * .(")")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(ymin, ymax)
)

points(
  x = signalOut$time,
  y = signalOut$getVariable("do") - signalIn$getVariable("do")
)

lapply(
  X = 1:length(ratioDic),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$do$dox-signalIn$getVariable("do"),
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
  text = bquote(.("Change in DO Conc (") * mg ~ L^-1 * .(")")),
  side = 4,
  line = 2.5,
  cex = 0.7
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
  ylab = bquote(.("DIC Conc (") * mu * mol ~ L^-1 * .(")")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(ymin, ymax)
)

points(
  x = signalOut$time,
  y = signalOut$getVariable("dic")
)

lapply(
  X = 1:length(ratioDic),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$dic$dic,
      col = modelcolors[index],
      lwd = linewidths[index]
    )
  }
)

text(
  x = signalOut$time[70],
  y = min(modelDoDic$output$dic$dic),
  adj = c(1,0),
  labels = bquote(.("DIC:C Respiration = 1.07"))
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

#plot change in DIC concentrations
plot(
  x = modelDoDic$downstreamTimePOSIX,
  y = modelDoDic$output$dic$dic-signalIn$getVariable("dic"),
  type = "n",
  ylab = bquote(.("Change in DIC Conc (") * mu * mol ~ L^-1 * .(")")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(ymin, ymax)
)

points(
  x = signalOut$time,
  y = signalOut$getVariable("dic") - signalIn$getVariable("dic")
)

lapply(
  X = 1:length(ratioDic),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$dic$dic-signalIn$getVariable("dic"),
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

#Plot pCO2 concentrations
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
  X = 1:length(ratioDic),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$dic$pCO2,
      col = modelcolors[index],
      lwd = linewidths[index]
    )
  }
)

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

#plot change in pCO2 concentration
plot(
  x = modelDoDic$downstreamTimePOSIX,
  y = modelDoDic$output$dic$pCO2-signalIn$getVariable("pCO2"),
  type = "n",
  ylab = bquote(.("Change in pCO2 Conc (") * mu * atm * .(")")),
  xlab = bquote(.("Time(")*MST*.(")")),
  ylim = c(ymin, ymax)
)

#create legend
legend(
  x = "bottomright",
  bty = "n",
  legend = c(ratioDic, "PAR Dist."),
  # legend = sprintf("Dic:C Fixation = %.2f", ratioDic),
  lty = c("solid"),
  col = c(modelcolors, "lightgoldenrod1"),
  lwd = c(linewidths, 15),
  title = "DIC:C Fixation = "
)

points(
  x = signalOut$time,
  y = signalOut$getVariable("pCO2") - signalIn$getVariable("pCO2")
)

lapply(
  X = 1:length(ratioDic),
  FUN = function(index){
    lines(
      x = modelsims[[index]]$downstreamTimePOSIX,
      y = modelsims[[index]]$output$dic$pCO2-signalIn$getVariable("pCO2"),
      col = modelcolors[index],
      lwd = linewidths[index]
    )
  }
)

dev.off()