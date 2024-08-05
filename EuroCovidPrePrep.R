# Worldwide Covid data up to 14/12/20 (start 31/12/19)
#
# Downloaded from:
# https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide
# 
# File available as: Data_COV_World.csv

X=read.csv("Data_COV_World.csv")

# 11 European Countries for analysis 

Euro_Country=c("Austria","Belgium","Denmark","France","Germany","Italy"
              ,"Norway","Spain","Sweden","Switzerland","United_Kingdom")

# The following sorts out the dates

DateS=matrix(0,ncol=4,nrow=350)

DateS[,1]=seq(1,350)

MonthLength=c(31,29,31,30,31,30,31,31,30,31,30,14) 
# Length of Months up to 14 December in 2020

DayX=31 # 31 December 2019
for(i in 1:12) DayX=c(DayX,seq(1,MonthLength[i]))
DateS[,2]=DayX # Day
DateS[,3]=c(12,rep(seq(1,12),MonthLength))
DateS[,4]=c(2019,rep(2020,349))

#
#
#

Euro_Death=matrix(0,ncol=15,nrow=350)

Euro_Death[,1:4]=DateS

for(i in 1:11)
{
  couni=4+i # Country counter
  Z=X[X$countriesAndTerritories==Euro_Country[i],]
  for(j in 1:350) 
  {  
    if(length(Z$deaths[(Z$day==DateS[j,2])&(Z$month==DateS[j,3])])==1)
      Euro_Death[j,couni]=Z$deaths[(Z$day==DateS[j,2])&(Z$month==DateS[j,3])]
  }
}

Euro_Death=as.data.frame(Euro_Death)

colnames(Euro_Death)=c("Day_Number","Day","Month","Year",Euro_Country)

# Writes to file European deaths by day for:
# Austria Belgium Denmark France Germany Italy Norway Spain Sweden Switzerland United_Kingdom
# from 31/12/2019 to 14/12/2020
# Data tails off in the last few days so probably best not to use beyond 30/11/2020

write.csv(Euro_Death,file="Euro_Death.csv",row.names = FALSE)

