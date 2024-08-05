#
# Simulation of birth-death model with (partial) detection
#
# Parameter - inputs for simulation
#    Pday - number of days of simulation from day of first case
#    betaday - birth parameters - vector of length Pday [k-1 -> k detection] 
#    gammaday - death parameters - vector of length Pday [k-1 -> k detection] 
#    detpday - detection probability - vector of length Pday [k-1 -> k detection] 
#
# Does not allow for parameters to change as epidemic proceeds - could edit for such a situation
#
# - Daily counts [number of cases each day]
# - Daily parameter settings
# - version 2 with daz [day of the week] identifier


Sim_P_daily=function(Pday,betaday,gammaday,detpday)
{
  t=runif(1) # Random time [day 1] of initial birth/infection
  # Note that the infection rates of day 1 apply until after the first detection
  x=1 # Counts the number of individuals alive - start with 1
  count=0 # counts the number of detected deaths
  cox=1 # Parameter day
  tm=0
  death=0 # Records time of (detected) deaths
  pops=0 # Records population size prior to detected deaths
  n=1 # Number of events to occur [first event birth of initial individual]
  xt=1 # Keeps track of number of individuals over time
  tim=0 # Keeps track of the times of events
  while((x>0)&(cox<(Pday+1))) 
# Epidemic proceeds whilst there are infectives [x>0] and less than k removals have been detected
  {
    alpha=betaday[cox]
    mu=gammaday[cox]
    d=detpday[cox] # Sets parameters after count detections have been observed
    E=alpha+mu # Rate at which events occur (per individual)
    p=alpha/E # Current probability of a birth
    q=1-p
    t=t+rexp(1,(x*E)) # time at which next event occurs
    if((count==0)|((t-tm)<1)) 
    {  
    sev=sample(3,1,replace=F,prob=c(p,d*q,(1-d)*q))
    if(sev==1) # event is a birth
    {
      x=x+1
    }
    if(sev==2) # event is a detected death
    {
      count=count+1
      if(count==1) # First detection occurs on first day
      {
        t=t-floor(t) # Sets time on first day
        tm=0
      }
      x=x-1
      death[count]=t  # time of death
      pops[count]=x+1 # Population size just before death
    }
    if(sev==3) # event is an undetected death
    {
      x=x-1
    }
    n=n+1
    xt[n]=x
    tim[n]=t # Keeps track of the population size at all times
    }
    else
    {
      tm=tm+1
      cox=cox+1
#      print(x)
    }
  }
  
  # Daily deaths calculated
#  tday=ceiling(t)
  daily=rep(0,Pday)
  for(j in 1:Pday) daily[j]=sum((death>(j-1))&(death<=j))
#  while(daily[1]==0) daily=daily[-1]
  
  tinitial=death[1] # time of the first detected death
  death=death-tinitial
  tim=tim-tinitial # Resetting times so first detected death at time 0
  list(death=death,tim=tim,xt=xt,pops=pops,daily=daily) # Output information
}

#
# Version 2 designed so that we can find day of week.
#

Sim_P_daily2=function(Pday,betaday,gammaday,detpday)
{
  t=runif(1) # Random time [day 1] of initial birth/infection
  # Note that the infection rates of day 1 apply until after the first detection
  x=1 # Counts the number of individuals alive - start with 1
  count=0 # counts the number of detected deaths
  cox=1 # Parameter day
  tm=0
  death=0 # Records time of (detected) deaths
  pops=0 # Records population size prior to detected deaths
  n=1 # Number of events to occur [first event birth of initial individual]
  xt=1 # Keeps track of number of individuals over time
  tim=0 # Keeps track of the times of events
  while((x>0)&(cox<(Pday+1))) 
    # Epidemic proceeds whilst there are infectives [x>0] and less than k removals have been detected
  {
    alpha=betaday[cox]
    mu=gammaday[cox]
    d=detpday[cox] # Sets parameters after count detections have been observed
    E=alpha+mu # Rate at which events occur (per individual)
    p=alpha/E # Current probability of a birth
    q=1-p
    t=t+rexp(1,(x*E)) # time at which next event occurs
    if((count==0)|((t-tm)<1)) 
    {  
      sev=sample(3,1,replace=F,prob=c(p,d*q,(1-d)*q))
      if(sev==1) # event is a birth
      {
        x=x+1
      }
      if(sev==2) # event is a detected death
      {
        count=count+1
        if(count==1) # First detection occurs on first day
        {
          t=t-floor(t) # Sets time on first day
          tm=0
        }
        x=x-1
        death[count]=t  # time of death
        pops[count]=x+1 # Population size just before death
      }
      if(sev==3) # event is an undetected death
      {
        x=x-1
      }
      n=n+1
      xt[n]=x
      tim[n]=t # Keeps track of the population size at all times
    }
    else
    {
      tm=tm+1
      cox=cox+1
      #      print(x)
    }
  }
  
  # Daily deaths calculated
  #  tday=ceiling(t)
  daily=rep(0,Pday)
  for(j in 1:Pday) daily[j]=sum((death>(j-1))&(death<=j))
  #  while(daily[1]==0) daily=daily[-1]
  
  daz=floor(death[1])+1
  
  tinitial=death[1] # time of the first detected death
  death=death-tinitial
  tim=tim-tinitial # Resetting times so first detected death at time 0
  list(death=death,tim=tim,xt=xt,pops=pops,daily=daily,daz=daz) # Output information
}


#
# Version 3 designed so that we can find day of week. - Doesn't set initial day to Monday!!
#

Sim_P_daily3=function(Pday,betaday,gammaday,detpday)
{
  t=runif(1) # Random time [day 1] of initial birth/infection
  # Note that the infection rates of day 1 apply until after the first detection
  x=1 # Counts the number of individuals alive - start with 1
  count=0 # counts the number of detected deaths
  cox=1 # Parameter day
  tm=0
  death=0 # Records time of (detected) deaths
  pops=0 # Records population size prior to detected deaths
  n=1 # Number of events to occur [first event birth of initial individual]
  xt=1 # Keeps track of number of individuals over time
  tim=0 # Keeps track of the times of events
  while((x>0)&(cox<(Pday+1))) 
    # Epidemic proceeds whilst there are infectives [x>0] and less than k removals have been detected
  {
    alpha=betaday[cox]
    mu=gammaday[cox]
    d=detpday[cox] # Sets parameters after count detections have been observed
    E=alpha+mu # Rate at which events occur (per individual)
    p=alpha/E # Current probability of a birth
    q=1-p
    t=t+rexp(1,(x*E)) # time at which next event occurs
    if((count==0)|((t-tm)<1)) 
    {  
      sev=sample(3,1,replace=F,prob=c(p,d*q,(1-d)*q))
      if(sev==1) # event is a birth
      {
        x=x+1
      }
      if(sev==2) # event is a detected death
      {
        count=count+1
        if(count==1) # First detection occurs on first day
        {
          t=t-floor(t) # Sets time on first day
#          tm=0
        }
        x=x-1
        death[count]=t  # time of death
        pops[count]=x+1 # Population size just before death
      }
      if(sev==3) # event is an undetected death
      {
        x=x-1
      }
      n=n+1
      xt[n]=x
      tim[n]=t # Keeps track of the population size at all times
    }
    else
    {
      tm=tm+1
      cox=cox+1
      #      print(x)
    }
  }
  
  # Daily deaths calculated
  #  tday=ceiling(t)
  daily=rep(0,Pday)
  for(j in 1:Pday) daily[j]=sum((death>(j-1))&(death<=j))
  #  while(daily[1]==0) daily=daily[-1]
  
  daz=floor(death[1])+1
  
  tinitial=death[1] # time of the first detected death
  death=death-tinitial
  tim=tim-tinitial # Resetting times so first detected death at time 0
  list(death=death,tim=tim,xt=xt,pops=pops,daily=daily,daz=daz) # Output information
}




