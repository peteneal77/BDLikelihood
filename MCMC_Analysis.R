X12=read.csv("Euro_12_FF3_1234567.csv")
X10=read.csv("MCMC Output/Euro_FF3_1234567.csv")
X8=read.csv("MCMC Output/Euro_FF3_8_1234567.csv")

Y=read.csv("MCMC Output/Euro5_muD25_FF3_1234567.csv")

boxplot(X10[,1:11]*10,x)

# Country set - Country[11]="UK"

muM=13/2

R0=muM*X10[,1:11]
colnames(R0)=Country
boxplot(R0,las=2)

Rt=muM*X10[,1:11]*exp(X10[,20:30])
for(i in 1:11)
{
  Rt[,i]=Rt[,i]*X10[,14]*X10[,15]*X10[,16]*X10[,17]*X10[,18]*X10[,19]
}

Rt[,9]=Rt[,9]/X10[,18]

colnames(Rt)=Country
boxplot(Rt,las=2)

#
#
#

CountryX=Country
CountryX[11]="United Kingdom"

ci=6

g0_8=X8[,ci]-X8[,12]
gt_8=X8[,ci]*X8[,14]*X8[,15]*X8[,16]*X8[,17]*X8[,18]*X8[,19]*exp(X8[,(19+ci)])-X8[,12]
if(ci==9) gt_8=X8[,ci]*X8[,14]*X8[,15]*X8[,16]*X8[,17]*X8[,19]*exp(X8[,(19+ci)])-X8[,12]

g0_10=X10[,ci]-X10[,12]
gt_10=X10[,ci]*X10[,14]*X10[,15]*X10[,16]*X10[,17]*X10[,18]*X10[,19]*exp(X10[,(19+ci)])-X10[,12]
if(ci==9) gt_10=X10[,ci]*X10[,14]*X10[,15]*X10[,16]*X10[,17]*X10[,19]*exp(X10[,(19+ci)])-X10[,12]

g0_12=X12[,ci]-X12[,12]
gt_12=X12[,ci]*X12[,14]*X12[,15]*X12[,16]*X12[,17]*X12[,18]*X12[,19]*exp(X12[,(19+ci)])-X12[,12]
if(ci==9)  gt_12=X12[,ci]*X12[,14]*X12[,15]*X12[,16]*X12[,17]*X12[,19]*exp(X12[,(19+ci)])-X12[,12]

G0=matrix(c(g0_8,g0_10,g0_12),ncol=3,byrow=F)
Gt=matrix(c(gt_8,gt_10,gt_12),ncol=3,byrow=F)

colnames(G0)=c("mu=1/8","mu=1/10","mu=1/12")
colnames(Gt)=c("mu=1/8","mu=1/10","mu=1/12")
tit=paste("Pre-intervention growth rate ",Country[ci])
boxplot(G0,las=1,main=tit)
tit2=paste("Post-intervention growth rate ",Country[ci])
boxplot(Gt,las=1,main=tit2)

#
#

ci=6

g0_8=X8[,ci]-X8[,12]
gt_8=X8[,ci]*X8[,14]*X8[,15]*X8[,16]*X8[,17]*X8[,18]*X8[,19]*exp(X8[,(19+ci)])-X8[,12]
if(ci==9) gt_8=X8[,ci]*X8[,14]*X8[,15]*X8[,16]*X8[,17]*X8[,19]*exp(X8[,(19+ci)])-X8[,12]

g0_10=X10[,ci]-X10[,12]
gt_10=X10[,ci]*X10[,14]*X10[,15]*X10[,16]*X10[,17]*X10[,18]*X10[,19]*exp(X10[,(19+ci)])-X10[,12]
if(ci==9) gt_10=X10[,ci]*X10[,14]*X10[,15]*X10[,16]*X10[,17]*X10[,19]*exp(X10[,(19+ci)])-X10[,12]

g0_12=X12[,ci]-X12[,12]
gt_12=X12[,ci]*X12[,14]*X12[,15]*X12[,16]*X12[,17]*X12[,18]*X12[,19]*exp(X12[,(19+ci)])-X12[,12]
if(ci==9)  gt_12=X12[,ci]*X12[,14]*X12[,15]*X12[,16]*X12[,17]*X12[,19]*exp(X12[,(19+ci)])-X12[,12]

g0_y=Y[,ci]-Y[,12]
gt_y=Y[,ci]*Y[,14]*Y[,15]*Y[,16]*Y[,17]*Y[,18]*Y[,19]*exp(Y[,(19+ci)])-Y[,12]
if(ci==9)  gt_y=Y[,ci]*Y[,14]*Y[,15]*Y[,16]*Y[,17]*Y[,19]*exp(Y[,(19+ci)])-Y[,12]


G0=matrix(c(g0_8,g0_10,g0_12,g0_y),ncol=4,byrow=F)
Gt=matrix(c(gt_8,gt_10,gt_12,gt_y),ncol=4,byrow=F)

colnames(G0)=c("mu=1/8","mu=1/10","mu=1/12","muD")
colnames(Gt)=c("mu=1/8","mu=1/10","mu=1/12","muD")
tit=paste("Pre-intervention growth rate ",Country[ci])
boxplot(G0,las=1,main=tit)
tit2=paste("Post-intervention growth rate ",Country[ci])
boxplot(Gt,las=1,main=tit2)

