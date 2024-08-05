#
# Section 5 Covid-19 plots
#

X=read.csv("MCMC Output/Euro_Covid_MCMC.csv") # Main MCMC output with mu=0.1

country=c("Austria", "Belgium", "Denmark","France", "Germany", "Italy", "Norway", "Spain", "Sweden", "Swiss", "UK")

#
# Figure 3
#

# a

R0=X[,1:11]*10 # Computes R_0 values since mu=0.1 fixed
R0=as.data.frame(R0)
colnames(R0)=country

boxplot(R0,las=2,main="R0: Prior intervention",ylab=expression(R[0]))

# b

Rt=R0

Rt[,9]=Rt[,9]/X[,18] # Simplest way of discounting absence of lockdown in Sweden

for(i in 1:6)
{
  for(j in 1:11) Rt[,j]=X[,(13+i)]*Rt[,j] # Effects of NPI
}

for(j in 1:11) Rt[,j]=exp(X[,(19+j)])*Rt[,j]  # Effect of last NPI - country specific

boxplot(Rt,las=2,main="Rt: Post intervention",ylab=expression(R[t]))

#
# Figure 4
#

X8=read.csv("MCMC Output/Euro_FF3_8_1234567.csv")  # mu=1/8
X12=read.csv("MCMC Output/Euro_FF3_12_1234567.csv") # mu=1/12

ci=11 # Country index - ci=11 - UK

# ci=5 # Germany

# a: Pre-intervention

pre=matrix(0,ncol=3,nrow=50000)

pre[,1]=X8[,ci]-X8[,12]
pre[,2]=X[,ci]-X[,12]
pre[,3]=X12[,ci]-X12[,12]

pre=as.data.frame(pre)

#colnames(pre)=c(expression(paste(mu,"=","1/8")),"1/10","1/12")

boxplot(pre,main=paste0("Pre-intervention growth rate ",country[ci]),ylab="Growth rate",xaxt="n")

axis(1,at=c(1,2,3),labels=c(expression(paste(mu,"=","1/8")),expression(paste(mu,"=","1/10")),expression(paste(mu,"=","1/12"))))
# b: Post-intervention

post=matrix(0,ncol=3,nrow=50000)

post[,1]=X8[,ci]
post[,2]=X[,ci]
post[,3]=X12[,ci]

if(ci==9)  # Simplest way of discounting absence of lockdown in Sweden
{
  post[,1]=post[,1]/X8[,18]
  post[,2]=post[,2]/X[,18]
  post[,3]=post[,3]/X12[,18]
}

for(i in 1:6)
{
  post[,1]=X8[,(13+i)]*post[,1]
  post[,2]=X[,(13+i)]*post[,2]
  post[,3]=X12[,(13+i)]*post[,3]
}

post[,1]=exp(X8[,(19+ci)])*post[,1]
post[,2]=exp(X[,(19+ci)])*post[,2]
post[,3]=exp(X12[,(19+ci)])*post[,3]

post[,1]=post[,1]-X8[,12]
post[,2]=post[,2]-X[,12]
post[,3]=post[,3]-X12[,12]

post=as.data.frame(post)

colnames(post)=c("1/8","1/10","1/12")

boxplot(post,main=paste0("Post-intervention growth rate ",country[ci]),ylab="Growth rate",xaxt="n")

axis(1,at=c(1,2,3),labels=c(expression(paste(mu,"=","1/8")),expression(paste(mu,"=","1/10")),expression(paste(mu,"=","1/12"))))
