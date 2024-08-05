#
# Calculation of reproduction number
#

Rt_calc=function(beta0,mu,thetaCon,betaF,DayCon)
  # Inputs beta0 - initial, mu - recovery rate, thetaCon - control measures, betaF - final beta 
  # DayCon - control measures active (1)/inactive (0)
{
  Len=length(DayCon[1,]) # Number of days to calculate R0 over
  NPX=length(DayCon[,1]) # Number of NPIs
  Rt=0 # value of Rt over time
  for(j in 1:Len) Rt[j]=(beta0/mu)*prod(thetaCon^DayCon[1:(NPX-1),j])*exp(betaF*DayCon[NPX,j])
  Rt
}
