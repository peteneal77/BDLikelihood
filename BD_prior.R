#
# Calculation of prior density
#
# Distributions: 1 - uniform, 2 - gamma, 3 - beta, 4 - normal
#

priorden=function(theta,theta_prior)
{
  logp=0 
  logp=logp+sum(log(dunif(theta[theta_prior[,1]==1],theta_prior[theta_prior[,1]==1,2],theta_prior[theta_prior[,1]==1,3])))
  logp=logp+sum(log(dgamma(theta[theta_prior[,1]==2],theta_prior[theta_prior[,1]==2,2],theta_prior[theta_prior[,1]==2,3])))
  logp=logp+sum(log(dbeta(theta[theta_prior[,1]==3],theta_prior[theta_prior[,1]==3,2],theta_prior[theta_prior[,1]==3,3])))
  logp=logp+sum(log(dnorm(theta[theta_prior[,1]==4],theta_prior[theta_prior[,1]==4,2],theta_prior[theta_prior[,1]==4,3])))
  logp
}

