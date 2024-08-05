#
# Rt_plot - plot through time with label option of when control measure is implemented 
#  
#

Rt_plot=function(RMat,coun,DayL,quan=c(0.025,0.25,0.5,0.75,0.975),SelDates,colQ,NPIX,NPIpos)
# RMat - Rt - Matrix, coun - Country
# DayL - Number of days from first death until the end of the period
# quan - quantiles (default = 2.5%, 25%, 50%, 75%, 97.5%)  
# SelDates - Selected dates for x-axis  
# colQ - Quantile colours (1,2,3) and text colour (4)
# NPIX - Named control measures and NPIpos - defines position (if NULL no labels)
{  
  QuanR=matrix(0,ncol=5,nrow=DayL)
  for(i in 1:DayL) QuanR[i,]=quantile(RMat[i,],quan)
  # Calculate quantiles

  plot(-1,-1,xlim=c(0,DayL),ylim=c(min(RMat),max(RMat)),
     ylab=expression(R[t]),xlab="",main=coun,xaxt='n')
  # Empty plot
  
  for(i in 1:(DayL-1))
  {
    polygon(c(i,i+1,i+1,i,i),
          c(QuanR[i,1],QuanR[i,1],QuanR[i,5],QuanR[i,5],QuanR[i,1])
          ,border=colQ[1],col=colQ[1])
  }
  # Outer quantiles - default equal-tailed 95% credible interval
  
  for(i in 1:(DayL-1))
  {
    polygon(c(i,i+1,i+1,i,i),
          c(QuanR[i,2],QuanR[i,2],QuanR[i,4],QuanR[i,4],QuanR[i,2])
          ,border=colQ[2],col=colQ[2])
  }
  # Inner quantiles - default equal-tailed 50% credible interval (inter-quartile range)
  
  for(i in 1:(DayL-1))
  {
    lines(c(i,i+1),rep(QuanR[i,3],2),col=colQ[3])
  }
  # Central quantile - default median
  
  axis(1,SelDates,at=seq(1, length(SelDates)),las=2,cex.axis=0.8) # Put in x-axis
  
  if(length(NPIpos)>0) # Labels on plots
  {
    locx=0
    for(i in 1:length(NPIpos)) locx[i]=DayL+1-sum(DayCon[NPIpos[i],])
  
    for(j in 1:length(NPIX))
    {
      text(locx[j],mean(QuanR[(locx[j]-1):locx[j],3]),NPIX[j],col=colQ[4])
    }
  }  
  
}


