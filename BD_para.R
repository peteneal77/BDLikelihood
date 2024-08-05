#
# Parameter convert
#


# Scenario 2 - time homogeneous with changepoint in alpha
#
# theta=(alpha1,alpha2,mu,d,L)

para2=function(theta,N)
{
  alpha=c(rep(theta[1],round(theta[5],0)),rep(theta[2],(N-round(theta[5],0))))
  mu=rep(theta[3],N)
  d=rep(theta[4],N)
  list(alpha=alpha,mu=mu,d=d)
}

# Scenario 3 - time-inhomogeneous with decaying birth rate
#
# theta=(beta,phi,gamma,delta)

para3=function(theta,N)
{
  alpha=theta[1]*theta[2]^(seq(0,(N-1))/N)
  mu=rep(theta[3],N)
  d=rep(theta[4],N)
  list(alpha=alpha,mu=mu,d=d)
}

# Scenario 4 - time homogeneous with changepoint in all parameters
#
# theta=(beta1,gamma1,delta1,beta2,gamma2,delta2,L)

para4=function(theta,N)
{
  alpha=c(rep(theta[1],theta[7]),rep(theta[4],(N-theta[7])))
  mu=c(rep(theta[2],theta[7]),rep(theta[5],(N-theta[7])))
  d=c(rep(theta[3],theta[7]),rep(theta[6],(N-theta[7])))
  list(alpha=alpha,mu=mu,d=d)
}

# Scenario 5 -  Time-inhomogeneous model with changing rates in all parameters
# beta decreasing, gamma and delta increasing
#
# theta=(beta,phi_beta,gamma,phi_gamma,delta,phi_delta)

para5=function(theta,N)
{
  alpha=theta[1]*theta[2]^(seq(0,(N-1))/N)
  mu=theta[3]*theta[4]^(seq(N,1,-1)/N)
  d=theta[5]*theta[6]^(seq(N,1,-1)/N)
  list(alpha=alpha,mu=mu,d=d)
}