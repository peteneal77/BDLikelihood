#
# R0 Calculator for European Covid data
#
#
#

ci=11# Selects a country
# Countries: 1-Austria 2-Belgium 3-Denmark 4-France 5-Germany 6-Italy 
#            7-Norway 8-Spain 9-Sweden 10-Switzerland 11-United_Kingdom


source("Rt_Calc.R") # Code to calculate Rt through time
source("Rt_Plot.R") # Code to plot Rt through time

library(RColorBrewer) # Colour pallette for producing R0 plots

# Run preparation file to aid with calculation of Control Measure implementations
# Sets up requirements for running preparation file.

muD=18  # Mean R -> D
alphaD=0 # Shape parameter for R -> D Gamma (alphaD, alphaD/muD). alphaD = 0 fixed infectious period.
ChooCoun=seq(1:11)    # Countries
ChCont=c(1,2,3,4,5,6,7) # Control Measures
# ChCont=c(2,5,6,7)

source("EuroCovidPrep.R") 

#

#X=read.csv("MCMC Output/Euro_FF3_1234567.csv") # Load in data

X=read.csv("MCMC Output/Euro_Covid_MCMC.csv")

NC=11 # Number of countries in analysis
NP=6 # Number of shared NPIs
NT=NP+1

Nrun=length(X[,1]) # Length of MCMC output
thinX=50 # thinning of MCMC output for R0 calculation.

thinS=seq(thinX,Nrun,thinX) # Thinned sequence of MCMC output

betaCountry=X[thinS,1:NC] # Initial value of beta (infection rate) in countries
mu=X[thinS,(NC+1)] # Removal rate
d=X[thinS,(NC+2)]  # Death rate
thetaControl=X[thinS,(NC+3):(NC+NP+2)] # Shared control measures
betaFinal=X[thinS,(NC+NP+3):(2*NC+NP+2)] # Beta final effect of control measures

#
# Dates
#

Euro_Dates=paste0(Euro_Data$Day[1],"/",Euro_Data$Month[1],"/",Euro_Data$Year[1])
for(i in 2:length(Euro_Data[,1])) Euro_Dates[i]=paste0(Euro_Data$Day[i],"/",Euro_Data$Month[i],"/",Euro_Data$Year[i])

DayX=126-DayE # Removal date of initial death (translated back) - DayE represents length

#
#
#


Tdaily1=EDeaths
Tdaily=as.matrix(Tdaily1)
Nday=length(Tdaily[,1])

Country_Control=Control_Int

#
# Setup control measures by day. 1 - Control measure in place, 0 - control measure not in place.
#

Day_Control=array(1,dim=c(NC,NT,Nday)) 
# dimension - country, control measures, day
for(i in 1:NC)
{
  for(j in 1:NT)
  {
    Day_Control[i,j,1:min(Country_Control[i,j],Nday)]=0
  }
}

#
# Selected country ci (ci=11 - United Kingdom) and associated control measures
#

DayL=DayE[ci]

DayCon=Day_Control[ci,,1:DayL]

RMat=matrix(0,ncol=length(thinS),nrow=DayL) 
# Matrix containing Rt estimates - each column uses different parameter set from the posterior

for(i in 1:length(thinS)) RMat[,i]=Rt_calc(betaCountry[i,ci],mu[i],thetaControl[i,],betaFinal[i,ci],DayCon)

#
#
#

NPI=c("Self_Isolate","Social_Distance","Public_Ban","School_Close","Lockdown","First_NPI","Final_NPI")  
# Lists NPIs

SelDates=Euro_Dates[(127-DayL):126] # Selected dates for plotting on x-axis to show affect on Rt

colQ=c("#C7E9C0","#74C476","#005A32","black") # Colours for Rt plots (Green pallette set here)

# brewer.pal(n = 8, name = "Greens")
# "#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D" "#238B45" "#005A32"

NPIX=c("First NPI","Social Distance","Lockdown and Final_NPI") 
# Text to place on graph to denote control measures enacted
NPos=c(6,2,5)

Rt_plot(RMat,Country[ci],DayL,quan=c(0.025,0.25,0.5,0.75,0.975),SelDates,colQ,NPIX,NPos) 
# Produces Rt plot for selected country

legend(x=c(40,68),y=c(3.6,2.8),legend=c("95% Credible interval","50% Credible interval","Posterior Median"),col=colQ[1:3],lty=rep(1,3),lwd=c(3,3,1))





