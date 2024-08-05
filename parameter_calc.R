#
# Function for computing birth, death and detection parameters for Birth-Death process
#

Control_para=function(theta,theta_Control,Day_Con,TD)  
  # Inputs - model parameters theta=(beta,mu,d)
  # theta_Control - effects of control measures [probability or allow greater than 1 - control measures]
  # Day_Con - When the control measures are introduced [simple introduction scheme]
  # Tdaily - Daily counts
{
  NT=length(theta_Control)  # Total number of controls
  K=sum(TD) # Total number of cases
  Pday=length(TD) # Total number of days
  DayK=matrix(rep(theta_Control,Pday),ncol=Pday,nrow=NT)
  betaday=theta[1]*exp(colSums(Day_Con*log(DayK))) 
  # Note initial infection rate is theta[1]
  alpha=c(rep(betaday,TD),betaday[Pday]) # Birth rate
  mu=rep(theta[2],(K+1)) # Death rate
  detprob=rep(theta[3],(K+1))# Detection probability
  list(alpha=alpha,mu=mu,detprob=detprob)
}




