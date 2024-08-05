#
# R0 Calculator for Simulated data
#
# Includes comparison with true Rt
#

Country=c("Country A","Country B","Country C","Country D","Country E")

ci=4 # Selects a country from 5

source("Rt_Calc.R") # Code to calculate Rt through time
source("Rt_Plot.R") # Code to plot Rt through time

library(RColorBrewer) # Colour pallette for producing R0 plots

#

X=read.csv("MCMC Output/EuroSimS2_vary3.csv") # Load in data

NC=5 # Number of countries in analysis
NP=5 # Number of shared NPIs
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

DayE=c(55,50,45,48,50) # Days after the first death to the end of the period

Sim_Dates="Day 1"
for(i in 2:max(DayE)) Sim_Dates[i]=paste0("Day ",i)

DayX=56-DayE # Removal date of initial death (translated back) - DayE represents length

#
#
#

Country_Control=read.csv("SimControl2.csv",row.names = 1) # Read in control measures

Tdaily1=read.csv("DailyDataS2.csv",row.names = 1) # Read in data
Tdaily=as.matrix(Tdaily1)
Nday=length(Tdaily1[,1]) # Total number of days in analysis

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
# Selected country ci and associated control measures
#

DayL=DayE[ci]

DayCon=Day_Control[ci,,1:DayL]

RMat=matrix(0,ncol=length(thinS),nrow=DayL) 
# Matrix containing Rt estimates - each column uses different parameter set from the posterior

for(i in 1:length(thinS)) RMat[,i]=Rt_calc(betaCountry[i,ci],mu[i],thetaControl[i,],betaFinal[i,ci],DayCon)

#
#
#

NPI=c("Control 1","Control 2","Control 3","Control 4","Control 5","Control Final")  
# Lists NPIs

SelDates=Sim_Dates[DayX[ci]:55] # Selected dates for plotting on x-axis to show affect on Rt

colQ=c("#C7E9C0","#74C476","#005A32","black") # Colours for Rt plots (Green pallette set here)

# brewer.pal(n = 8, name = "Greens")
# "#F7FCF5" "#E5F5E0" "#C7E9C0" "#A1D99B" "#74C476" "#41AB5D" "#238B45" "#005A32"

NPIX=c()
# Text to place on graph to denote control measures enacted
NPos=c()

Rt_plot(RMat,Country[ci],DayL,quan=c(0.025,0.25,0.5,0.75,0.975),SelDates,colQ,NPIX,NPos) 
# Produces Rt plot for selected country

# Comparison with true Rt
# Plots true Rt in purple
#

R0_Country=c(4,3.5,4.5,3.8,4.2) 
theta_Control=c(0.9,0.8,0.8,0.75,0.35)
mu=0.1
beta_Final=c(-0.15, -0.01,  0.07, -0.24,0.17)

Rt_true=Rt_calc(R0_Country[ci]*mu,mu,theta_Control,beta_Final[ci],DayCon) # Calculates true Rt value

for(i in 2:(length(Rt_true)-1))
{
  lines(c(1,2),rep(Rt_true[1],2),col="purple")
  if(Rt_true[i-1]==Rt_true[i]) lines(c(i,i+1),rep(Rt_true[i],2),col="purple")
  if(Rt_true[i-1]!=Rt_true[i]) lines(c(i,i+1),rep(Rt_true[i],2),col="purple")
}



