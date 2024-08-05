#
#
#
#


source("aGSEI_functions.R")  # Gaussian approximation of birth-death likelihood
source("GSEI_functions.R")   # Birth-death likelihood
source("EGSEI_functions.R")  # Exact GSE likelihood


Abik=c(13,7,2,3,0,0,1,4,5,3,2,0,2,0,5,3,1,4,0,1,1,1,2,0,1,5,0,5,5)

#
# Find mode close to 0.09,0.0775
#


betax=seq(0.0875,0.0925,0.0005)
gammax=seq(0.0750,0.0800,0.0005)


WX90=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX90[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}

which(WX90==max(WX90))
# 26th entry row 4, column 3 (3-1)*11+4
# Max on grid at beta=0.089, gamma=0.076

betax=seq(0.0885,0.0895,0.0001)
gammax=seq(0.0755,0.0765,0.0001)


WX90=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX90[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}
which(WX90==max(WX90))
# 71st entry row 5, column 7 (7-1)*11+5
# Max on grid at beta=0.0889, gamma=0.0761
# Mode is -60.01865


#
# Find mode close to 0.16,0.004
# Initial run with beta on (0.1575,0.1625) gave max at beta=0.1625
#

betax=seq(0.16,0.165,0.0005)
gammax=seq(0.0035,0.0045,0.0001)


WX90=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX90[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}

which(WX90==max(WX90))
# 40th entry row 7, column 4 (4-1)*11+7
# Max on grid at beta=0.163, gamma=0.0038

betax=seq(0.1625,0.1635,0.0001)
gammax=seq(0.0037,0.0039,0.00002)


WX90=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX90[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}

which(WX90==max(WX90))
# 70th entry row 4, column 7 (7-1)*11+4
# Max on grid at beta=0.1628, gamma=0.00382
# Mode is -59.01405

#
# TT=95 - day 95
#

TT=95

betax=seq(0.09,0.095,0.0005)
gammax=seq(0.0775,0.0825,0.0005)


WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 95th entry row 7, column 9 (9-1)*11+7
# Max on grid at beta=0.093, gamma=0.0815

betax=seq(0.0925,0.0935,0.0001)
gammax=seq(0.081,0.082,0.0001)

WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 85th entry row 8, column 8 (8-1)*11+8
# Max on grid at beta=0.0932, gamma=0.0817
# Mode at -60.2581

betax=seq(0.1625,0.1675,0.0005)
gammax=seq(0.003,0.004,0.0001)

WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 76th entry row 10, column 6 (7-1)*11+10
# Max on grid at beta=0.167, gamma=0.0035

betax=seq(0.167,0.168,0.0001)
gammax=seq(0.0034,0.0036,0.00002)

WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 97th entry row 9, column 9 (9-1)*11+9
# Max on grid at beta=0.1678, gamma=0.00356
# Mode at -60.6742


#
# TT=94 - day 94
#

TT=94

betax=seq(0.09,0.095,0.0005)
gammax=seq(0.0775,0.0825,0.0005)


WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 83th entry row 6, column 8 (8-1)*11+6
# Max on grid at beta=0.0925, gamma=0.081

betax=seq(0.092,0.093,0.0001)
gammax=seq(0.0805,0.0815,0.0001)

WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 51th entry row 7, column 5 (5-1)*11+7
# Max on grid at beta=0.0926, gamma=0.0809
# Mode at -60.2253

betax=seq(0.1625,0.1675,0.0005)
gammax=seq(0.003,0.004,0.0001)

WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 76th entry row 10, column 6 (7-1)*11+10
# Max on grid at beta=0.167, gamma=0.0035

betax=seq(0.1665,0.1675,0.0001)
gammax=seq(0.0035,0.0037,0.00002)

WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 68th entry row 2, column 7 (7-1)*11+2
# Max on grid at beta=0.1666, gamma=0.00362
# Mode at -60.35149



#
# TT=93 - day 93
#

TT=93

betax=seq(0.09,0.095,0.0005)
gammax=seq(0.0775,0.0825,0.0005)


WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 60th entry row 5, column 6 (6-1)*11+5
# Max on grid at beta=0.092, gamma=0.08

betax=seq(0.0915,0.0925,0.0001)
gammax=seq(0.0795,0.0805,0.0001)

WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 48th entry row 4, column 5 (5-1)*11+4
# Max on grid at beta=0.0918, gamma=0.0799
# Mode at -60.1862

betax=seq(0.1625,0.1675,0.0005)
gammax=seq(0.003,0.004,0.0001)

WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 83th entry row 6, column 8 (8-1)*11+6
# Max on grid at beta=0.165, gamma=0.0037

betax=seq(0.165,0.166,0.0001)
gammax=seq(0.0036,0.0038,0.00002)

WX=matrix(0,nrow=length(betax),ncol=length(gammax))

Et0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,TT,120)
  }
  print(i)
}

which(WX==max(WX))
# 43th entry row 10, column 4 (4-1)*11+10
# Max on grid at beta=0.1659, gamma=0.00366
# Mode at -60.02418

#
# Refined search for modes using the approximate likelihood.
#
#
#


betax=seq(0.0925,0.0975,0.0005)
gammax=seq(0.0750,0.0800,0.0005)

AX=matrix(0,nrow=length(betax),ncol=length(gammax))

for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    AX[i,j]=alike_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}
which(AX==max(AX))
# 72th entry row 6, column 7 (7-1)*11+6
# Max on grid at beta=0.095, gamma=0.078


betax=seq(0.0945,0.0955,0.0001)
gammax=seq(0.0775,0.0785,0.0001)

AX=matrix(0,nrow=length(betax),ncol=length(gammax))

for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    AX[i,j]=alike_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}
which(AX==max(AX))
# 47th entry row 3, column 5 (5-1)*11+3
# Max on grid at beta=0.0947, gamma=0.0779
# -60.09181  R0=1.2156

betax=seq(0.1375,0.1425,0.0005)
gammax=seq(0.004,0.005,0.0001)

AX=matrix(0,nrow=length(betax),ncol=length(gammax))

for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    AX[i,j]=alike_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}
which(AX==max(AX))
# 76th entry row 10, column 7 (7-1)*11+10
# Max on grid at beta=0.142, gamma=0.0046

betax=seq(0.14,0.141,0.0001)
gammax=seq(0.0045,0.0047,0.00002)

AX=matrix(0,nrow=length(betax),ncol=length(gammax))

for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    AX[i,j]=alike_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}
which(AX==max(AX))
# 29th entry row 7, column 3 (3-1)*11+7
# Max on grid at beta=0.1406, gamma=0.00454
#  -59.9182; R0 = 30.969

