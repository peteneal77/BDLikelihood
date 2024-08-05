#
# Code for computing likelihood contour plots in the paper.
#
# Exact and approximate likelihoods calculated for T=90 and T=46.99
#
#

# Load in functions for computing log-likelihood

source("aGSEI_functions.R")  # Gaussian approximation of birth-death likelihood
source("GSEI_functions.R")   # Birth-death likelihood
source("EGSEI_functions.R")  # Exact GSE likelihood

# Inter-removal times for Abakiliki data

Abik=c(13,7,2,3,0,0,1,4,5,3,2,0,2,0,5,3,1,4,0,1,1,1,2,0,1,5,0,5,5)

# Parameter values alpha (beta) and mu (gamma) at which the likelihood is computed

betax=seq(0.0025,0.25,0.0025)
gammax=c(seq(0.0005,0.012,0.0005),seq(0.0125,0.2,0.0025))

# Exact likelihood for T=90
# Pre-calculated in "Analysis/Abak_grid_E_90.csv"

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
Et1=proc.time()

# Et1-Et0 gives a measure of the time taken.

# Load in pre-calculated data.

WX90=read.csv(file="Analysis/Abak_grid_E_90.csv")
WX90=as.matrix(WX90)

# Contour plot using grayscale

gg=gray.colors(n=23,start=1,end=0,gamma=1)

filled.contour(betax,gammax,exp(WX90),xlab=expression(alpha),ylab=expression(mu),
               levels=seq(0,2.3*10^-26,10^-27),col=gg,main="Exact: T=90")

# Figure produced: Analysis/Abak_T90_E.png

filled.contour(betax,gammax,exp(WX90),xlab=expression(alpha),ylab=expression(mu),main="Exact: T=90",levels=seq(0,2.3*10^-26,10^-27)) 

# Figure produced: Analysis/Abak_T90_E_col.png

# Exact likelihood for T=46.99
# Pre-calculated in "Analysis/Abak_grid_E_46b.csv"

WX46b=matrix(0,nrow=length(betax),ncol=length(gammax))

for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    WX46b[i,j]=Elike_GSEI(betax[i],gammax[j],Abik,46.99,120)
  }
  print(i)
}

# Load in pre-calculated data.

WX46b=read.csv(file="Analysis/Abak_grid_E_46b.csv")
WX46b=as.matrix(WX46b)

# Contour plot using grayscale

gg=gray.colors(n=23,start=1,end=0,gamma=1)

filled.contour(betax,gammax,exp(WX46b),xlab=expression(alpha),ylab=expression(mu),
               levels=seq(0,3.45*10^-13,1.5*10^-14),col=gg,main="Exact: T=46.99")

# Analysis/Abak_T46b_E.png

filled.contour(betax,gammax,exp(WX46b),xlab=expression(alpha),ylab=expression(mu),
               main="Exact: T=46.99",levels=seq(0,3.45*10^-13,1.5*10^-14))


# Analysis/Abak_T46b_E_col.png

# Approximate likelihood calculations using the Gaussian approximation of the birth-death process.

AX90=matrix(0,nrow=length(betax),ncol=length(gammax))

At0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    AX90[i,j]=alike_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}
At1=proc.time()

# At1-At0 gives a measure of the time taken.

# Approximate likelihood calculations using the birth-death process.

BX90=matrix(0,nrow=length(betax),ncol=length(gammax))

Bt0=proc.time()
for(i in 1:length(betax))
{
  for(j in 1:length(gammax))
  {
    BX90[i,j]=like_GSEI(betax[i],gammax[j],Abik,90,120)
  }
  print(i)
}
Bt1=proc.time()

# Bt1-Bt0 gives a measure of the time taken.

# Contour plot using grayscale (Gaussian approximationg of birth-death process)

gg=gray.colors(n=23,start=1,end=0,gamma=1)

filled.contour(betax,gammax,exp(AX90),xlab=expression(alpha),ylab=expression(mu),
               levels=seq(0,2.3*10^-26,10^-27),col=gg,main="Approximate likelihood")

# Analysis/Abak_T90_A.png

filled.contour(betax,gammax,exp(AX90),xlab=expression(alpha),ylab=expression(mu),main="Approximate likelihood",
               levels=seq(0,2.3*10^-26,10^-27)) 

# Figure produced: Analysis/Abak_T90_A_col.png


# Contour plot using grayscale (Birth-death process)

gg=gray.colors(n=23,start=1,end=0,gamma=1)

filled.contour(betax,gammax,exp(BX90),xlab=expression(alpha),ylab=expression(mu),
               levels=seq(0,2.3*10^-26,10^-27),col=gg,main="Birth-death likelihood")

# Analysis/Abak_T90_B.png

filled.contour(betax,gammax,exp(BX90),xlab=expression(alpha),ylab=expression(mu),main="Birth-death likelihood",
               levels=seq(0,2.3*10^-26,10^-27)) 

# Figure produced: Analysis/Abak_T90_B_col.png