#
# Table 7 and 8 computation
#

A=read.csv("MCMC Output/EuroSimS2_vary3.csv")

growth_mean_0=0 # Initial growth rate mean
growth_sd_0=0 # Initial growth rate standard deviation
growth_mean_t=0 # Final growth rate mean
growth_sd_t=0 # Final growth rate standard deviation

R_mean_0=0 # Initial R_0 mean
R_sd_0=0 # Initial R_0  standard deviation
R_mean_t=0 # Final R_t mean
R_sd_t=0 # Final R_t standard deviation

NPI_effect=A[,8]*A[,9]*A[,10]*A[,11]*A[,12] # Effects of all NPIs

for(i in 1:5)
{
  growth_mean_0[i]=mean(A[,i]-A[,6])
  growth_sd_0[i]=sd(A[,i]-A[,6])
  growth_mean_t[i]=mean(A[,i]*NPI_effect*exp(A[,(i+12)])-A[,6])
  growth_sd_t[i]=sd(A[,i]*NPI_effect*exp(A[,(i+12)])-A[,6])
  R_mean_0[i]=mean(A[,i]/A[,6])
  R_sd_0[i]=sd(A[,i]/A[,6])
  R_mean_t[i]=mean(A[,i]*NPI_effect*exp(A[,(i+12)])/A[,6])
  R_sd_t[i]=sd(A[,i]*NPI_effect*exp(A[,(i+12)])/A[,6])
}

round(growth_mean_0,4)
round(growth_sd_0,4)
round(growth_mean_t,4)
round(growth_sd_t,4)

round(R_mean_0,4)
round(R_sd_0,4)
round(R_mean_t,4)
round(R_sd_t,4)

#
# Simulation values
#

alpha=c(0.4,0.35,0.45,0.38,0.42)
mu=0.1
zeta=c(0.9,0.8,0.8,0.75,0.35)
xi=c(-0.15,-0.01,0.07,-0.24,0.17)

npis=prod(zeta)

pre_growth=alpha-mu
post_growth=alpha*npis*exp(xi)-mu
pre_R0=alpha/mu
post_Rt=alpha*npis*exp(xi)/mu

pre_growth
round(post_growth,4)
pre_R0
round(post_Rt,4)