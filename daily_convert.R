#
# Daily data 
#
# This file contains two functions
# DailyCon - Generates removal (death times) given daily counts
# InterUp - Updates Nup removal times selected uniformly at random.


#
# Convert daily data to detected death times
# 
#

DailyCon=function(daily,excl)
# daily - list of daily counts
# excl - if excl = 1 delete final [incomplete] daily counts
{
  if(excl==1) daily=daily[-(length(daily))] # excludes incomplete day
  tday=length(daily) # number of days of daily data
  a=runif(daily[1],0,1)
  for(i in 2:tday) a=c(a,runif(daily[i],(i-1),i)) 
  # generate times of detections based on uniform across the days
  ao=a[order(a,decreasing=F)] # Puts times of detections in order
  acu=ao # This gives the actual times which is useful for data augmentation updates
  ao=ao-min(ao) # Put first time as zero
  lao=length(ao)
  ainter=ao[2:lao]-ao[1:(lao-1)] # Creates inter-arrival times
#  ainter # Output inter-arrival times of detected deaths
  list(ainter=ainter,acu=acu)
}



#
# Updating inter-arrival times
#

InterUp=function(TA,Tdcase,Nup)
{
  TA_new=TA # Proposed new times - start by setting equal to current times
  indi=sample(length(TA),Nup,replace=F) # Selects Nup random times to update - indi is the index
  u=runif(Nup) # Uniform draws of times over the day of death
  TA_new[indi]=Tdcase[indi]+u-1 # Times of deaths of updated individuals
  ao=TA_new[order(TA_new,decreasing=F)] # Puts times of detections in order
  acu=ao # This gives the actual times which is useful for data augmentation updates
  ao=ao-min(ao) # Put first time as zero
  lao=length(ao)
  ainter=ao[2:lao]-ao[1:(lao-1)]
  list(ainter=ainter,acu=acu)
}

