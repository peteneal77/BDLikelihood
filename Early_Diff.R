# Exploring approximation of E[R_{k-1} r_k (\tilde{t}_k)^{R_{k-1}}] by
# \xi_k^{(1)} (\varphi_k) for k small.
#

# k=4

# h3 (= h_3 (\tilde{t}_3)) can take values on the range 0 to 1
# r4 (= r_4 (\tilde{t}_4)) can take values on the range q to 1
# varphi4 = log(r4)
#
# Given h3: P(R_3 = 2) = (1-h3)  and P(R_3=3) = h3
# eta3 = 2 + h3 and sigma3 = h3 (1-h3)
#

# LogDif4 computes the maximum difference between log E[R_{k-1} r_k (\tilde{t}_k)^{R_{k-1}}] 
# and log \xi_k^{(1)} (\varphi_k) given q is the lower bound for r4.

LogDif4=function(q)
{
  h3=seq(0,1,0.001)

  eta3=2+h3
  sigma3=h3*(1-h3)

  r4=seq(q,1,0.001)
  varphi4=log(r4)

  E4=matrix(0,nrow=length(h3),ncol=length(r4))
  A4=matrix(0,nrow=length(h3),ncol=length(r4))

  for(i in 1:length(h3))
  {
    E4[i,]=(1-h3[i])*(2*r4^2)+h3[i]*(3*r4^3)
    A4[i,]=(eta3[i]+varphi4*sigma3[i])*exp(varphi4*eta3[i]+0.5*sigma3[i]*varphi4^2)
  }
  
  c(min(log(E4)-log(A4)),max(log(E4)-log(A4)))
}

LogDif4(0.4)
LogDif4(0.3)
LogDif4(0.2)
LogDif4(0.1)

# > LogDif4(0.4)
# [1] -0.003166808  0.006652274
# > LogDif4(0.3)
#[1] -0.003166808  0.006697162
#> LogDif4(0.2)
#[1] -0.006982917  0.019722893
#> LogDif4(0.1)
#[1] -0.06514778  0.09233568
# > exp(LogDif4(0.1))
#[1] 0.936929 1.096733
#> exp(LogDif4(0.4))
#[1] 0.9968382 1.0066744

# k=5

# h3 (= h_3 (\tilde{t}_3)) can take values on the range 0 to 1
# phi(t_4) can take values on the range 0 to 1.
# r5 (= r_4 (\tilde{t}_4)) can take values on the range q to 1
# varphi4 = log(r4)
#
# Given h3: P(R_3 = 2) = (1-h3)  and P(R_3=3) = h3
# eta3 = 2 + h3 and sigma3 = h3 (1-h3)
#

# LogDif4 computes the maximum difference between log E[R_{k-1} r_k (\tilde{t}_k)^{R_{k-1}}] 
# and log \xi_k^{(1)} (\varphi_k) given q is the lower bound for r4.

Mmat4=function(phi,q)
{
  r=q+(1-q)*phi
  h=phi/r
  S=seq(0,2)
  M=matrix(0,nrow=2,ncol=3)
  for(i in 1:2) M[i,]=(i+1)*dbinom(S,i,h)*r^i
  M
}


LogDif5=function(q)
{
  h3=seq(0,1,0.005)
  phi4=seq(0,1,0.005)
  h4=phi4/(q+(1-q)*phi4)
  varphi4=log(q+(1-q)*phi4)
  
  eta3=2+h3
  sigma3=h3*(1-h3)
  
  r5=seq(q,1,0.005)
  varphi5=log(r5)
  
  E5=array(0,dim=c(length(h3),length(phi4),length(r5)))
  A5=array(0,dim=c(length(h3),length(phi4),length(r5)))
 
  for(i in 1:length(h3))
  {
    B3=c(1-h3[i],h3[i])
    for(k in 1:length(phi4))
    {
      M=Mmat4(phi4[k],q)
      B4=B3%*%M
      B4=B4/sum(B4)
      E5[i,k,]=B4[1]*(2*r5^2)+B4[2]*(3*r5^3)+B4[3]*(4*r5^4)
      
      eta4=2+h4[k]*(-1+eta3[i]+varphi4[k]*sigma3[i]+sigma3[i]/(eta3[i]+varphi4[k]*sigma3[i]))
      sigma4=h4[k]*(1-h4[k])*(-1+eta3[i]+varphi4[k]*sigma3[i]+sigma3[i]/(eta3[i]+varphi4[k]*sigma3[i]))+
        h4[k]^2*sigma3[i]*(1-sigma3[i]/(eta3[i]+varphi4[k]*sigma3[i]))
      A5[i,k,]=(eta4+varphi5*sigma4)*exp(varphi5*eta4+0.5*sigma4*varphi5^2)
    }
  }

  c(min(log(E5)-log(A5)),max(log(E5)-log(A5)))
}

LogDif5(0.4)
LogDif5(0.3)
LogDif5(0.2)
LogDif5(0.1)

# > LogDif5(0.4)
#[1] -0.004564235  0.011677936
#> LogDif5(0.3)
#[1] -0.004683583  0.013991649
#> LogDif5(0.2)
#[1] -0.01631615  0.05541796
#> LogDif5(0.1)
#[1] -0.1202969  0.2275937

#
# The following computes the maximum difference for k=4 and k=5 combined as a function of q,
# h3 (phi3), h4 (phi4), r5 (varphi5)
# 
#

Mmat4=function(phi,q)
{
  r=q+(1-q)*phi
  h=phi/r
  S=seq(0,2)
  M=matrix(0,nrow=2,ncol=3)
  for(i in 1:2) M[i,]=(i+1)*dbinom(S,i,h)*r^i
  M
}


LogDif45=function(q)
{
  h3=seq(0,1,0.005)
  phi4=seq(0,1,0.005)
  h4=phi4/(q+(1-q)*phi4)
  r4=q+(1-q)*phi4
  varphi4=log(q+(1-q)*phi4)
  
  eta3=2+h3
  sigma3=h3*(1-h3)
  
  r5=seq(q,1,0.005)
  varphi5=log(r5)
  
  E5=array(0,dim=c(length(h3),length(phi4),length(r5)))
  A5=array(0,dim=c(length(h3),length(phi4),length(r5)))
  
  E4=array(0,dim=c(length(h3),length(phi4)))
  A4=array(0,dim=c(length(h3),length(phi4)))
  
  for(i in 1:length(h3))
  {
    B3=c(1-h3[i],h3[i])
    for(k in 1:length(phi4))
    {
      M=Mmat4(phi4[k],q)
      B4=B3%*%M
      B4=B4/sum(B4)
      E4[i,k]=B3[1]*(2*r4[k]^2)+B3[2]*(3*r4[k]^3)
      E5[i,k,]=(B4[1]*(2*r5^2)+B4[2]*(3*r5^3)+B4[3]*(4*r5^4))*E4[i,k]
      
      eta4=2+h4[k]*(-1+eta3[i]+varphi4[k]*sigma3[i]+sigma3[i]/(eta3[i]+varphi4[k]*sigma3[i]))
      sigma4=h4[k]*(1-h4[k])*(-1+eta3[i]+varphi4[k]*sigma3[i]+sigma3[i]/(eta3[i]+varphi4[k]*sigma3[i]))+
        h4[k]^2*sigma3[i]*(1-sigma3[i]/(eta3[i]+varphi4[k]*sigma3[i]))
      A4[i,k]=(eta3[i]+varphi4[k]*sigma3[i])*exp(varphi4[k]*eta3[i]+0.5*sigma3[i]*varphi4[k]^2)
      A5[i,k,]=(eta4+varphi5*sigma4)*exp(varphi5*eta4+0.5*sigma4*varphi5^2)*A4[i,k]
    }
  }
  
  c(min(log(E5)-log(A5)),max(log(E5)-log(A5)))
}

LogDif45(0.4)
LogDif45(0.3)
LogDif45(0.2)
LogDif45(0.1)


#> LogDif45(0.4)
#[1] -0.005295876  0.015580769
#> LogDif45(0.3)
#[1] -0.006002601  0.015523085
#> LogDif45(0.2)
#[1] -0.01433809  0.05541796
#> LogDif45(0.1)
#[1] -0.1400411  0.2259259


