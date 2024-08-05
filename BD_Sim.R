#
# Simulation Birth-Death Process
#

# Inputs: N - Process is observed until N^th death (or extinction)
# alpha, mu, d - vectors of length N defining birth rates, death rates and detection probabilities
#
#

BD=function(N,alpha,mu,d)
{
  y=1
  detect=0
  tim=0
  t=0
  while((detect<N)&(y>0))
  {
    pro=c(alpha[detect+1],mu[detect+1])
    tim=tim+rexp(1,y*sum(pro))
    dec=sample(2,1,prob=pro,replace=T)
    if(dec==1) y=y+1
    if(dec==2)
    {
      y=y-1
      see=rbinom(1,1,d[detect+1])
      if(see==1)
      {
        detect=detect+1
        t[detect]=tim
      }
    }
  }
  inter=0
  if(detect>1) inter=t[2:detect]-t[1:(detect-1)]
  inter # inter-arrival times of detected deaths - returns 0 if 0/1 detected deaths
}

