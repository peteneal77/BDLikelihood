# MCMC algorithm for Covid model with control measures
# - Multiple independent countries
#
# 24/10/22 - Code adapted to include input Covariance matrix input

Covid_MCMC=function(Nrun,Tdata,DayZ,Country_Control,
                     theta,theta_fix,theta_prob,theta_int
                     ,prior,sigma,Nup)
  # Nrun - Number of iterations per run (all but last are burn-in)  
  # Tdata - daily data [matrix]
  # DayZ - Number of days of data on number of cases for each country
  # Country_Control - Day after which control measures are introduced in each country
  # theta - model parameters (initial values)
  # theta_fix - 1 if parameter changes; 0 if parameter fixed [theta_fix - same length theta]
  # theta_prob - probability parameters (detection) - parameters which are probabilities
  # prior - priors
  # sigma - RWM standard deviation - be smarter.
  # Initial Sigma = sigma I and Sigma updated after each burn-in run. 
  # Nup - Number of death times to update
{
  Nloop=length(Nrun) # Number of cycles of the MCMC algorithm, last run stored.
  
  CouN=length(DayZ) # Total number of countries
  ConT=length(Country_Control[1,]) # The number of control measures
  
  Pday=max(DayZ) # Maximum number of days
  
  beta=theta[1:CouN] # The first CouN parameters are beta for each country
  gamma=theta[CouN+1] # Common recovery rate
  d=theta[CouN+2] # Detection (death) probability same across countries
  theta_Control=theta[(CouN+3):(CouN+1+ConT)] # Control measure parameters (subtract 1 for country effect)
  beta_Final=theta[(CouN+ConT+2):(2*CouN+ConT+1)] # Final control measure in each country
  
  Mpara=sum(theta_fix)  # Number of parameters to update 
  lambda=2.4^2/Mpara
  # Scaling weight for covariance matrix  
  
  SIGMA=as.matrix(sigma)
  if(length(SIGMA[,1])==1) SIGMA=diag(rep(sigma,length(theta))) # Initial proposal variance-covariance matrix
  SIGMA1=0.1*diag(diag(SIGMA)) # Allows for small moves to stop MCMC getting stuck, independent moves in each component
  
  if(length(Nup)==1) Nup=rep(Nup,CouN) 
  # Allows for a single value to be inserted for number of times to update - Nup=0 means no update
  
  # Create actual times and inter-arrival times
  KN=colSums(Tdaily) # Number of cases in each country
  KX=max(KN) # maximum number of cases in a country
  T_Actual=matrix(0,ncol=CouN,nrow=KX)  # Actual times of events
  T_Inter=matrix(0,ncol=CouN,nrow=(KX-1)) # Inter-arrival times of events 
  # Above matrices needed for updates and likelihood calculations
  for(ia in 1:CouN)
  {
    TX=DailyCon(Tdaily[,ia],0)
    T_Inter[1:(KN[ia]-1),ia]=TX$ainter 
    T_Actual[1:KN[ia],ia]=TX$acu 
  }
  Tnew_Inter=T_Inter # This will be used for updates of the inter-arrival times
  Tnew_Actual=T_Actual 
  
  # Day of the cases
  T_Cases=matrix(0,ncol=CouN,nrow=KX)
  for(ia in 1:CouN) T_Cases[1:KN[ia],ia]=rep(seq(1,DayZ[ia]),Tdaily[1:DayZ[ia],ia])
  
  
  Day_Control=array(1,dim=c(CouN,ConT,Pday)) 
  # dimension - country, control measures, day
  for(i in 1:CouN)
  {
    for(j in 1:ConT)
    {
      Day_Control[i,j,1:min(Country_Control[i,j],Pday)]=0
    }
  }
  
  
  # Calculate Likelihood 
  # LC - likelihood contribution from each Country - LL = Total [Note all on log scale]  
  LC=0  
  for(ia in 1:CouN)
  {
    Para=Control_para(c(beta[ia],gamma,d),c(theta_Control,exp(beta_Final[ia])),Day_Control[ia,,1:DayZ[ia]],Tdaily[1:DayZ[ia],ia])
    # Sets parameters
    alpha=Para$alpha  # alpha_k
    mu=Para$mu # mu_k = gamma (fixed)
    detprob=Para$detprob # detprob_k = d (fixed)
    LC[ia]=alike_BD(alpha,mu,detprob,T_Inter[(1:(KN[ia]-1)),ia],(DayZ[ia]-T_Actual[KN[ia],ia])) 
    # Calculates likelihood in country 
  }
  LL=sum(LC) # LL - log-likelihood over all countries
  
  # Prior - LP (log)
  theta_ip=theta_prob+theta_int # if entry 0 - then neither prob or integer
  LP=sum(dgamma(theta[theta_ip==0],prior[theta_ip==0,1],prior[theta_ip==0,2],log=T))
  LP=LP+sum(dbeta(theta[theta_prob==1],prior[theta_prob==1,1],prior[theta_prob==1,2],log=T))
  LP=LP+sum(dnorm(theta[theta_prob==2],prior[theta_prob==2,1],prior[theta_prob==2,2],log=T))
  # theta_prob - records whether a parameter is positive (0), probability (1) or on the whole real line (2)
  
  # Posterior (log) equal to sum of LL and LP.
  
  #
  # MCMC loops
  #
  
  for(il in 1:Nloop)
  {
    print(il) # Shows which loop the MCMC is on
    output=matrix(0,ncol=length(theta),nrow=Nrun[il])
    Acc=0 # Parameter acceptance number
    AcT=rep(0,CouN) # Time change acceptance number
    coxy=0
    for(i in 1:Nrun[il])
    {
      coxy=coxy+1
      if(coxy==100) # 1000
      {
        print(i) # Counter to provide an indication of how far through the MCMC run is.
        coxy=0
      }
      #   Update death times
      for(ia in 1:CouN) # Consider each country in turn [this could be parallelised]
      {  
        # Update removal times
        if(Nup[ia]>0) # Only update any of the times if Nup[ia] positive
        {
          TXnew=InterUp(T_Actual[1:KN[ia],ia],T_Cases[1:KN[ia],ia],Nup[ia]) # Propose new removal times
          Tnew_Inter[1:(KN[ia]-1),ia]=TXnew$ainter
          Tnew_Actual[1:KN[ia],ia]=TXnew$acu
          Para=Control_para(c(beta[ia],gamma,d),c(theta_Control,exp(beta_Final[ia])),
                            Day_Control[ia,,1:DayZ[ia]],Tdaily[1:DayZ[ia],ia])
          # Update parameters based on new removal times.
          alpha=Para$alpha
          mu=Para$mu
          detprob=Para$detprob
          # Calculate likelihood
          Lnew=alike_BD(alpha,mu,detprob,Tnew_Inter[(1:(KN[ia]-1)),ia],
                        (DayZ[ia]-Tnew_Actual[KN[ia],ia]))
          u=log(runif(1))
          if(u<(Lnew-LC[ia])) # No effect of prior when updating infection times.
          {
            T_Inter[1:(KN[ia]-1),ia]=TXnew$ainter
            T_Actual[1:KN[ia],ia]=TXnew$acu 
            AcT[ia]=AcT[ia]+1
            LC[ia]=Lnew
          }
        }
      }
      
      theta=c(beta,gamma,d,theta_Control,beta_Final) # Overkill but confirms the parameters
  
      
      #   Update parameters      
      thetanew=theta+theta_fix*mvrnorm(1,rep(0,length(theta)),SIGMA)
      
      tt=runif(1)
      if(tt<0.05) thetanew=theta+theta_fix*mvrnorm(1,rep(0,length(theta)),SIGMA1) # Small perturbation
      # SIGMA1 is used to avoid getting stuck.
      
      # Set parameters to their new values given thetanew
      # This assists with understanding what each parameter represents
      betanew=thetanew[1:CouN] # The first CouN parameters are R0 for each country
      gammanew=thetanew[CouN+1]
      dnew=thetanew[CouN+2]
      theta_Controlnew=thetanew[(CouN+3):(CouN+1+ConT)]
      beta_Finalnew=thetanew[(CouN+2+ConT):(2*CouN+1+ConT)] # Effect of final control measure
      # Control measure parameters
      #      thetanew[theta_int==1]=round(thetanew[theta_int==1]) 
      #    No integer "parameters" so this can be omitted      
      
      LPnew=sum(dgamma(thetanew[theta_ip==0],prior[theta_ip==0,1],prior[theta_ip==0,2],log=T))
      LPnew=LPnew+sum(dbeta(thetanew[theta_prob==1],prior[theta_prob==1,1],prior[theta_prob==1,2],log=T))
      LPnew=LPnew+sum(dnorm(thetanew[theta_prob==2],prior[theta_prob==2,1],prior[theta_prob==2,2],log=T))
      # New prior weights
      
      if ((min(thetanew[theta_prob<2])>0)&(max(thetanew[theta_prob==1])<1)) 
      # More general put in conditions on theta, this ensures positive values are positive and probabilities are between 0 and 1 
      { 
        LCnew=0  
        for(ia in 1:CouN)
        {
          Para=Control_para(c(betanew[ia],gammanew,dnew),c(theta_Controlnew,exp(beta_Finalnew[ia]))
                            ,Day_Control[ia,,1:DayZ[ia]],Tdaily[1:DayZ[ia],ia])
          alpha=Para$alpha
          mu=Para$mu
          detprob=Para$detprob
          LCnew[ia]=alike_BD(alpha,mu,detprob,T_Inter[(1:(KN[ia]-1)),ia],(DayZ[ia]-T_Actual[KN[ia],ia]))
        }
        u=log(runif(1))
        if(u<(sum(LCnew)+LPnew-sum(LC)-LP)) # Compare posteriors with new parameters versus old parameters
        {
          # Parameter updates
          theta=thetanew
          beta=betanew
          beta_Final=beta_Finalnew
          gamma=gammanew
          d=dnew
          theta_Control=theta_Controlnew
          Acc=Acc+1
          LC=LCnew
          LP=LPnew
        }
      }
      
#      print(LCnew)
      
      output[i,]=theta
      #      print(c(theta[theta_fix==1],Q$Liklog))
    }
    Vop=rep(0,length(theta))
    Nlow=round(0.5*Nrun[il],0)+1
    # Use last half of MCMC run to choose proposal variance for theta 
    # 
    for(j in 1:length(theta)) if(theta_fix[j]==1) Vop[j]=var(output[Nlow:Nrun[il],j])
    SIGMA[theta_fix==1,theta_fix==1]=lambda*(0.95*var(output[Nlow:Nrun[il],theta_fix==1])
                                             +0.05*diag(Vop[theta_fix==1]))
    SIGMA1=diag(0.1*Vop) # Diagonal matrix
    # Print output from run
    print(Acc)
    print(SIGMA)
    print(AcT)
    print(colMeans(output))
  }
  
  #
  # MCMC output
  #
  output
  
}