#
#  European Covid-19 data preparation
#
#  Data contained in: Euro_Death.csv
#  Data pre-processed from:
# https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide
# Data available on many countries up to 14/12/20
# File available as: Data_COV_World.csv
#
#  Data consists of case counts and reported deaths by date from 31/12/19 to 14/12/20 (350 days) for 11 countries.
#  The 11 countries are those considered in Flaxman et al. (2020) and is extended version of:
#   https://reshare.ukdataservice.ac.uk/854380/ 
#  File: covid19-ecdc-data.csv
#  Note that there are some minor differences in counts
#
# Uncomment lines 38, 39, 166, 176 to run outside of EuroAnalysis

Euro_Data=read.csv("Euro_Death.csv") # Read in data

Euro_Death=Euro_Data[,5:15] # Reduce data down to deaths in the 11 countries

# Obtain a list of country names

Country=colnames(Euro_Death) # Names of countries

# Model: S->I->R->D
# Assumption that time in infectious (I state) is Exp (\mu).
# Main analysis based on setting \mu = 0.1

# R -> D 
# Following Flaxman et al. (2020) we take \mu_D = 18 from removal to death
# We consider two main possibilities although the code below can be altered to any 
# Gamma (\alpha_D, \alpha_D/\mu_D) distribution
# 1. Fixed length time from removal to death. \mu_D = 18 (Coded as \alpha_D = 0.)
# 2. \alpha_D = 5. This results in a coefficient of variation of 0.449 in line 
# with that chosen in Flaxman et al. (2020)

# muD=18 # Mean R -> D
# alphaD=5 # Shape parameter of Gamma distribution. alphaD = 0 (corresponds to fixed distribution)
betaD=alphaD/muD
# alphaD=0

RDeaths=matrix(0,ncol=length(Euro_Death[1,]),nrow=length(Euro_Death[,1])) 

# If alphaD>0, the cdf of the Gamma (alphaD,betaD) (with rounding to the nearest day) is used to translate
# the removal times back from the current date.
#

if (alphaD>0)
{  
  for(i in 1:length(Euro_Death[1,]))
  {
    for(ia in 1:length(Euro_Death[,1]))
    {
      if(Euro_Death[ia,i]>0)
      {
        pp=seq(1,Euro_Death[ia,i])/(Euro_Death[ia,i]+1)
        rp=round((ia-1-qgamma(pp,alphaD,betaD)),0)
        for(ib in 1:Euro_Death[ia,i])
        {
          RDeaths[rp[ib],i]=RDeaths[rp[ib],i]+1
        }
      }
    }
  }
}

if (alphaD==0)
{
  for(i in 1:length(Country))  RDeaths[1:(length(Euro_Death[,1])-muD),i]=Euro_Death[(muD+1):length(Euro_Death[,1]),i]
}

# Hu finds the day on which the first removal occurs to be in line with the death data and R->D

colnames(RDeaths)=Country # Sets the names of RDeaths

Hu=0
roD=rowSums(RDeaths)
for(i in 1:length(roD))
{
  if(roD[i]>0)
  {
    Hu=i
    break
  }
}

# Our analysis is based on data up to and including 03/05/2020. 
# Note that we base it on removal times up to this date so the analysis includes some death data from after
# 03/05/2020.

# Day 125 is 03/05/2020 and the first death is on day Hu. Thus the data we use is from day Hu to day 125.
# DDeaths stores this data.

DDeaths=RDeaths[Hu:125,] 

# For each country we want to find the date (relative to Hu) on which the first death occurs, and this is given 
# in HDeaths

HDeaths=rep(0,length(Country))
for(i in 1:length(Country))
{
  j=1
  while(DDeaths[j,i]==0)
  {
    j=j+1
  }
  HDeaths[i]=j
}

# DayH gives the number of days of data from the first death in a country to 03/05/2022
# The data in DDeaths is adjusted for each country so row 1 is the day of the first death.

DayH=length(DDeaths[,1])+1-HDeaths
for(i in 1:length(Country))
{
  if(HDeaths[i]>1)
  {
    DDeaths[,i]=c(DDeaths[HDeaths[i]:length(DDeaths[,1]),i],rep(0,length(DDeaths[,1])-DayH[i]))
  }
}


# Non-pharmaceutical interventions (NPIs)
# NPIs implemented in March 2020 as recorded in Flaxman et al. (2020)

March_Int=read.csv("Intervention_March.csv",row.names = 1) # Interventions 

# Day 62 - 1st March
# Day 125 - 3rd May

# Matrix of day of intervention as days after Hu, the first removal associated with a death. 

Int_TimeAdj=March_Int+61-Hu 
# Translating by adding on 61 for days prior to 1st March
# -Hu for date of first death.

# Adjust day of intervention to be relative to the first removal associated with a death. 
for(i in 1:length(Country))
{
  Int_TimeAdj[i,]=March_Int[i,]+61-HDeaths[i]-Hu
}

# If an intervention time is before the first death set the effect equal to 0.

Int_TimeAdj[Int_TimeAdj<0]=0

# For each country identify the day of the first and last intervention in March 2020. 
# Note that an intervention day set equal to 1000 means this intervention has not been used. 
# In the data this is only Lockdown in Sweden is not implemented.

Int_TimeAdj$First=0
Int_TimeAdj$Final=0

for(i in 1:length(Country)) 
{
  Int_TimeAdj$First[i]=min(Int_TimeAdj[i,1:5])
  Int_TimeAdj$Final[i]=max(Int_TimeAdj[i,Int_TimeAdj[i,]<1000])
}

# The final part of the preparation is to choose the countries and control measures to consider:
# Countries: 1-Austria 2-Belgium 3-Denmark 4-France 5-Germany 6-Italy 
#            7-Norway 8-Spain 9-Sweden 10-Switzerland 11-United_Kingdom
# NPIs: 1-Self_Isolate 2-Social_Distance 3-Public_Ban 4-School_Close 5-Lockdown 6-First 7-Final

# ChooCoun=seq(1:11)
# Chosen countries

EDeaths=DDeaths[,ChooCoun]
# Removal numbers per day for chosen countries

DayE=DayH[ChooCoun]
HEDeaths=HDeaths[ChooCoun]
# Relative day of first removal and number of days until 03/05/2020 respectively.

# ChCont=c(1,2,3,4,5,6,7) 
# Choose control measures (NPIs)

Control_Int=Int_TimeAdj[ChooCoun,ChCont]
# Control measures (NPIs) for chosen countries and interventions
