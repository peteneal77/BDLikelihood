#
# Code to compute posterior means, variances and correlations using numerical integration over grid.
#

X90=read.csv("Abak_grid_E_90.csv")

betax=seq(0.0025,0.25,0.0025)
gammax=c(seq(0.0005,0.012,0.0005),seq(0.0125,0.2,0.0025)) # 25 0.0005-0.0125, 75 after

TX90=as.matrix(X90)

AA90=exp(TX90-max(TX90))

for(i in 1:100)
{
  for(j in 1:25)
  {
    AA90[i,j]=AA90[i,j]*dgamma(betax[i],10,500/6)*dgamma(gammax[j],10,100)*0.2 # for finer grid on first points
  }
  for(j in 26:100)
  {
    AA90[i,j]=AA90[i,j]*dgamma(betax[i],10,500/6)*dgamma(gammax[j],10,100)
  }
}

sum(rowSums(AA90)*betax)/sum(AA90)
sum(colSums(AA90)*gammax)/sum(AA90)

alpha_bar=sum(rowSums(AA90)*betax)/sum(AA90)
mu_bar=sum(colSums(AA90)*gammax)/sum(AA90)

alpha_sd=sqrt(sum(rowSums(AA90)*betax^2)/sum(AA90)-alpha_bar^2)
mu_sd=sqrt(sum(colSums(AA90)*gammax^2)/sum(AA90)-mu_bar^2)

Ebetagamma=AA90/sum(AA90)
for(i in 1:100)
{
  for(j in 1:100)
  {
    Ebetagamma[i,j]=Ebetagamma[i,j]*betax[i]*gammax[j]
  }
}

alpha_mu_cor=(sum(Ebetagamma)-alpha_bar*mu_bar)/(alpha_sd*mu_sd)

alpha_bar
mu_bar
alpha_sd
mu_sd
alpha_mu_cor

#> alpha_bar
#[1] 0.1071805
#> mu_bar
#[1] 0.09236288
#> alpha_sd
#[1] 0.02305163
#> mu_sd
#[1] 0.02076402
#> alpha_mu_cor
#[1] 0.4541561


#
#
#

X46b=read.csv("Abak_grid_E_46b.csv")

betax=seq(0.0025,0.25,0.0025)
gammax=c(seq(0.0005,0.012,0.0005),seq(0.0125,0.2,0.0025)) # 25 0.0005-0.0125, 75 after

TX46b=as.matrix(X46b)

AA46b=exp(TX46b-max(TX46b))

for(i in 1:100)
{
  for(j in 1:25)
  {
    AA46b[i,j]=AA46b[i,j]*dgamma(betax[i],10,500/6)*dgamma(gammax[j],10,100)*0.2 # for finer grid on first points
  }
  for(j in 26:100)
  {
    AA46b[i,j]=AA46b[i,j]*dgamma(betax[i],10,500/6)*dgamma(gammax[j],10,100)
  }
}

alpha_bar=sum(rowSums(AA46b)*betax)/sum(AA46b)
mu_bar=sum(colSums(AA46b)*gammax)/sum(AA46b)

alpha_sd=sqrt(sum(rowSums(AA46b)*betax^2)/sum(AA46b)-alpha_bar^2)
mu_sd=sqrt(sum(colSums(AA46b)*gammax^2)/sum(AA46b)-mu_bar^2)

Ebetagamma=AA46b/sum(AA46b)
for(i in 1:100)
{
  for(j in 1:100)
  {
    Ebetagamma[i,j]=Ebetagamma[i,j]*betax[i]*gammax[j]
  }
}

alpha_mu_cor=(sum(Ebetagamma)-alpha_bar*mu_bar)/(alpha_sd*mu_sd)

alpha_bar
mu_bar
alpha_sd
mu_sd
alpha_mu_cor

#> alpha_bar
#[1] 0.1103015
#> mu_bar
#[1] 0.08388515
#> alpha_sd
#[1] 0.02565957
#> mu_sd
#[1] 0.02446432
#> alpha_mu_cor
#[1] 0.2700584