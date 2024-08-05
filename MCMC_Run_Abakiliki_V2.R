#
# Running MCMC code for analysis of Abakiliki data
#

# Load in MCMC algorithm and required functions:
source("MCMC_IncGSE.R")  # Treats observed data as removal times
source("MCMC_daily_IncGSE.R")  # Treats observed data as aggregated data

# Load in Abakiliki data - inter-removal times.
Abik=c(13,7,2,3,0,0,1,4,5,3,2,0,2,0,5,3,1,4,0,1,1,1,2,0,1,5,0,5,5)

Abik_rem=rep(0,101)
Abik_rem[1]=1
count=1
for(i in 1:29) 
{
  count=count+Abik[i]
  Abik_rem[count]=Abik_rem[count]+1
}

# All MCMC output available in MCMC folder.

# T=46.99,47,80,90,100,1000(infinity) - burn=3,bits=1000,nits=10000
# T=46.99,47,90 - burn=3,bits=5000,nits=50000 (Longer runs)

QB47=MCMC_IncGSE(Abik,47,like_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QA47=MCMC_IncGSE(Abik,47,alike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QE47=MCMC_IncGSE(Abik,47,Elike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)




AB47=MCMC_daily_IncGSE(Abik_rem[1:47],47,like_GSEI,120,0.1,0.087,3,1000,10000,0.01,5,c(1,1),nub=10,lamb=500/6,nug=10,lamg=100)
AA47=MCMC_daily_IncGSE(Abik_rem[1:47],47,alike_GSEI,120,0.1,0.087,3,1000,10000,0.01,5,c(1,1),nub=10,lamb=500/6,nug=10,lamg=100)
AE47=MCMC_daily_IncGSE(Abik_rem[1:47],47,Elike_GSEI,120,0.1,0.087,3,1000,10000,0.01,5,c(1,1),nub=10,lamb=500/6,nug=10,lamg=100)

#write.csv(QB47,row.names=FALSE,file="MCMC/Abak_47_BD.csv")
#write.csv(QA47,row.names=FALSE,file="MCMC/Abak_47_approx.csv")
#write.csv(QE47,row.names=FALSE,file="MCMC/Abak_47_exact.csv")

#write.csv(AB47,row.names=FALSE,file="MCMC/Abak_D47_BD.csv")
#write.csv(AA47,row.names=FALSE,file="MCMC/Abak_D47_approx.csv")
#write.csv(AE47,row.names=FALSE,file="MCMC/Abak_D47_exact.csv")


#QA47=read.csv(file="MCMC/Abak_47_approx.csv")
#QB47=read.csv(file="MCMC/Abak_47_BD.csv")
#QE47=read.csv(file="MCMC/Abak_47_exact.csv")

#QA47=as.matrix(QA47)
#QB47=as.matrix(QB47)
#QE47=as.matrix(QE47)

Res47=matrix(0,ncol=5,nrow=3)

colnames(Res47)=c("Mean alpha","Mean mu","sd alpha","sd mu","cor(alpha,mu)")

Res47[1,1:2]=colMeans(QA47)
Res47[1,3]=sd(QA47[,1])
Res47[1,4]=sd(QA47[,2])
Res47[1,5]=cor(QA47[,1],QA47[,2])
Res47[2,1:2]=colMeans(QB47)
Res47[2,3]=sd(QB47[,1])
Res47[2,4]=sd(QB47[,2])
Res47[2,5]=cor(QB47[,1],QB47[,2])
Res47[3,1:2]=colMeans(QE47)
Res47[3,3]=sd(QE47[,1])
Res47[3,4]=sd(QE47[,2])
Res47[3,5]=cor(QE47[,1],QE47[,2])

#      Mean alpha    Mean mu   sd alpha      sd mu cor(alpha,mu)
#[1,]  0.1183999 0.08046293 0.02820822 0.02356408     0.2751180
#[2,]  0.1170105 0.07922061 0.02724181 0.02379063     0.2452757
#[3,]  0.1141952 0.07982052 0.02645023 0.02313379     0.2782876
# Approx, BD, Exact

QB46b=MCMC_IncGSE(Abik,46.99,like_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QA46b=MCMC_IncGSE(Abik,46.99,alike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QE46b=MCMC_IncGSE(Abik,46.99,Elike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)


Abik_rem_b=c(Abik_rem[1:46],0)

AB46b=MCMC_daily_IncGSE(Abik_rem_b,46.99,like_GSEI,120,0.1,0.087,3,5000,50000,0.01,5,c(1,1),nub=10,lamb=500/6,nug=10,lamg=100)
AA46b=MCMC_daily_IncGSE(Abik_rem_b,46.99,alike_GSEI,120,0.1,0.087,3,5000,50000,0.01,5,c(1,1),nub=10,lamb=500/6,nug=10,lamg=100)
AE46b=MCMC_daily_IncGSE(Abik_rem_b,46.99,Elike_GSEI,120,0.1,0.087,3,5000,50000,0.01,5,c(1,1),nub=10,lamb=500/6,nug=10,lamg=100)

write.csv(AB46b,row.names=FALSE,file="MCMC/Abak_D46bL_BD.csv")
write.csv(AA46b,row.names=FALSE,file="MCMC/Abak_D46bL_approx.csv")
write.csv(AE46b,row.names=FALSE,file="MCMC/Abak_D46bL_exact.csv")

#write.csv(QB46b,row.names=FALSE,file="MCMC/Abak_46b_BD.csv")
#write.csv(QA46b,row.names=FALSE,file="MCMC/Abak_46b_approx.csv")
#write.csv(QE46b,row.names=FALSE,file="MCMC/Abak_46b_exact.csv")

Res46b=matrix(0,ncol=5,nrow=3)

colnames(Res46b)=c("Mean alpha","Mean mu","sd alpha","sd mu","cor(alpha,mu)")

Res46b[1,1:2]=colMeans(QA46b)
Res46b[1,3]=sd(QA46b[,1])
Res46b[1,4]=sd(QA46b[,2])
Res46b[1,5]=cor(QA46b[,1],QA46b[,2])
Res46b[2,1:2]=colMeans(QB46b)
Res46b[2,3]=sd(QB46b[,1])
Res46b[2,4]=sd(QB46b[,2])
Res46b[2,5]=cor(QB46b[,1],QB46b[,2])
Res46b[3,1:2]=colMeans(QE46b)
Res46b[3,3]=sd(QE46b[,1])
Res46b[3,4]=sd(QE46b[,2])
Res46b[3,5]=cor(QE46b[,1],QE46b[,2])

#     Mean alpha    Mean mu   sd alpha      sd mu cor(alpha,mu)
#[1,]  0.1125013 0.08282005 0.02622035 0.02423797     0.2548407
#[2,]  0.1145931 0.08550293 0.02773972 0.02552127     0.2462794
#[3,]  0.1105847 0.08359836 0.02592695 0.02363639     0.2593803

QB90=MCMC_IncGSE(Abik,90,like_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QA90=MCMC_IncGSE(Abik,90,alike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QE90=MCMC_IncGSE(Abik,90,Elike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)


AB90=MCMC_daily_IncGSE(Abik_rem[1:90],90,like_GSEI,120,0.1,0.087,3,5000,50000,0.01,5,c(1,1),nub=10,lamb=500/6,nug=10,lamg=100)
AA90=MCMC_daily_IncGSE(Abik_rem[1:90],90,alike_GSEI,120,0.1,0.087,3,5000,50000,0.01,5,c(1,1),nub=10,lamb=500/6,nug=10,lamg=100)
AE90=MCMC_daily_IncGSE(Abik_rem[1:90],90,Elike_GSEI,120,0.1,0.087,3,5000,50000,0.01,5,c(1,1),nub=10,lamb=500/6,nug=10,lamg=100)

write.csv(AB90,row.names=FALSE,file="MCMC/Abak_D90L_BD.csv")
write.csv(AA90,row.names=FALSE,file="MCMC/Abak_D90L_approx.csv")
write.csv(AE90,row.names=FALSE,file="MCMC/Abak_D90L_exact.csv")


#write.csv(QB90,row.names=FALSE,file="MCMC/Abak_90_BD.csv")
#write.csv(QA90,row.names=FALSE,file="MCMC/Abak_90_approx.csv")
#write.csv(QE90,row.names=FALSE,file="MCMC/Abak_90_exact.csv")

#QA90=read.csv(file="MCMC/Abak_90_approx.csv")
#QB90=read.csv(file="MCMC/Abak_90_BD.csv")
#QE90=read.csv(file="MCMC/Abak_90_exact.csv")

#QA90=as.matrix(QA90)
#QB90=as.matrix(QB90)
#QE90=as.matrix(QE90)

Res90=matrix(0,ncol=5,nrow=3)

colnames(Res90)=c("Mean alpha","Mean mu","sd alpha","sd mu","cor(alpha,mu)")

Res90[1,1:2]=colMeans(QA90)
Res90[1,3]=sd(QA90[,1])
Res90[1,4]=sd(QA90[,2])
Res90[1,5]=cor(QA90[,1],QA90[,2])
Res90[2,1:2]=colMeans(QB90)
Res90[2,3]=sd(QB90[,1])
Res90[2,4]=sd(QB90[,2])
Res90[2,5]=cor(QB90[,1],QB90[,2])
Res90[3,1:2]=colMeans(QE90)
Res90[3,3]=sd(QE90[,1])
Res90[3,4]=sd(QE90[,2])
Res90[3,5]=cor(QE90[,1],QE90[,2])

#    Mean alpha    Mean mu   sd alpha      sd mu cor(alpha,mu)
#[1,]  0.1111400 0.09259390 0.02375487 0.02084274     0.3972052
#[2,]  0.1109329 0.09214457 0.02359599 0.02053128     0.3917148
#[3,]  0.1079846 0.09311511 0.02296936 0.02143423     0.4588312


QB80=MCMC_IncGSE(Abik,80,like_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QA80=MCMC_IncGSE(Abik,80,alike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QE80=MCMC_IncGSE(Abik,80,Elike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)

#write.csv(QB80,row.names=FALSE,file="MCMC/Abak_80_BD.csv")
#write.csv(QA80,row.names=FALSE,file="MCMC/Abak_80_approx.csv")
#write.csv(QE80,row.names=FALSE,file="MCMC/Abak_80_exact.csv")

#QA80=read.csv(file="MCMC/Abak_80_approx.csv")
#QB80=read.csv(file="MCMC/Abak_80_BD.csv")
#QE80=read.csv(file="MCMC/Abak_80_exact.csv")

#QA80=as.matrix(QA80)
#QB80=as.matrix(QB80)
#QE80=as.matrix(QE80)

Res80=matrix(0,ncol=5,nrow=3)

colnames(Res80)=c("Mean alpha","Mean mu","sd alpha","sd mu","cor(alpha,mu)")

Res80[1,1:2]=colMeans(QA80)
Res80[1,3]=sd(QA80[,1])
Res80[1,4]=sd(QA80[,2])
Res80[1,5]=cor(QA80[,1],QA80[,2])
Res80[2,1:2]=colMeans(QB80)
Res80[2,3]=sd(QB80[,1])
Res80[2,4]=sd(QB80[,2])
Res80[2,5]=cor(QB80[,1],QB80[,2])
Res80[3,1:2]=colMeans(QE80)
Res80[3,3]=sd(QE80[,1])
Res80[3,4]=sd(QE80[,2])
Res80[3,5]=cor(QE80[,1],QE80[,2])

QB100=MCMC_IncGSE(Abik,100,like_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QA100=MCMC_IncGSE(Abik,100,alike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QE100=MCMC_IncGSE(Abik,100,Elike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)

#write.csv(QB100,row.names=FALSE,file="MCMC/Abak_100_BD.csv")
#write.csv(QA100,row.names=FALSE,file="MCMC/Abak_100_approx.csv")
#write.csv(QE100,row.names=FALSE,file="MCMC/Abak_100_exact.csv")

#QA100=read.csv(file="MCMC/Abak_100_approx.csv")
#QB100=read.csv(file="MCMC/Abak_100_BD.csv")
#QE100=read.csv(file="MCMC/Abak_100_exact.csv")

#QA100=as.matrix(QA100)
#QB100=as.matrix(QB100)
#QE100=as.matrix(QE100)

Res100=matrix(0,ncol=5,nrow=3)

colnames(Res100)=c("Mean alpha","Mean mu","sd alpha","sd mu","cor(alpha,mu)")

Res100[1,1:2]=colMeans(QA100)
Res100[1,3]=sd(QA100[,1])
Res100[1,4]=sd(QA100[,2])
Res100[1,5]=cor(QA100[,1],QA100[,2])
Res100[2,1:2]=colMeans(QB100)
Res100[2,3]=sd(QB100[,1])
Res100[2,4]=sd(QB100[,2])
Res100[2,5]=cor(QB100[,1],QB100[,2])
Res100[3,1:2]=colMeans(QE100)
Res100[3,3]=sd(QE100[,1])
Res100[3,4]=sd(QE100[,2])
Res100[3,5]=cor(QE100[,1],QE100[,2])

#> Res80
#Mean alpha    Mean mu   sd alpha      sd mu cor(alpha,mu)
#[1,]  0.1071654 0.08308673 0.02352621 0.02163635     0.4226526
#[2,]  0.1085472 0.08493041 0.02437631 0.02180280     0.4404398
#[3,]  0.1040188 0.08414777 0.02243796 0.02059952     0.4897506
#> Res100
#Mean alpha    Mean mu   sd alpha      sd mu cor(alpha,mu)
#[1,]  0.1128736 0.09500838 0.02401796 0.02059150     0.3799843
#[2,]  0.1125586 0.09559278 0.02343980 0.02015779     0.3570123
#[3,]  0.1096849 0.09640843 0.02352108 0.02108499     0.4280151

QB1000=MCMC_IncGSE(Abik,1000,like_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QA1000=MCMC_IncGSE(Abik,1000,alike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
QE1000=MCMC_IncGSE(Abik,1000,Elike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)

#write.csv(QB1000,row.names=FALSE,file="MCMC/Abak_1000_BD.csv")
#write.csv(QA1000,row.names=FALSE,file="MCMC/Abak_1000_approx.csv")
#write.csv(QE1000,row.names=FALSE,file="MCMC/Abak_1000_exact.csv")

#QA1000=read.csv(file="MCMC/Abak_1000_approx.csv")
#QB1000=read.csv(file="MCMC/Abak_1000_BD.csv")
#QE1000=read.csv(file="MCMC/Abak_1000_exact.csv")

#QA1000=as.matrix(QA1000)
#QB1000=as.matrix(QB1000)
#QE1000=as.matrix(QE1000)

Res1000=matrix(0,ncol=5,nrow=3)

colnames(Res1000)=c("Mean alpha","Mean mu","sd alpha","sd mu","cor(alpha,mu)")

Res1000[1,1:2]=colMeans(QA1000)
Res1000[1,3]=sd(QA1000[,1])
Res1000[1,4]=sd(QA1000[,2])
Res1000[1,5]=cor(QA1000[,1],QA1000[,2])
Res1000[2,1:2]=colMeans(QB1000)
Res1000[2,3]=sd(QB1000[,1])
Res1000[2,4]=sd(QB1000[,2])
Res1000[2,5]=cor(QB1000[,1],QB1000[,2])
Res1000[3,1:2]=colMeans(QE1000)
Res1000[3,3]=sd(QE1000[,1])
Res1000[3,4]=sd(QE1000[,2])
Res1000[3,5]=cor(QE1000[,1],QE1000[,2])

#> Res1000
#     Mean alpha    Mean mu   sd alpha      sd mu cor(alpha,mu)
#[1,]  0.1128001 0.09596812 0.02515912 0.02080596     0.3725554
#[2,]  0.1137305 0.09606721 0.02471504 0.02060236     0.3665421
#[3,]  0.1094680 0.09648065 0.02297605 0.02038819     0.4218273

QB1000a=MCMC_IncGSE(Abik,1000,like_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=1,lamb=0.001,nug=1,lamg=0.001)
QA1000a=MCMC_IncGSE(Abik,1000,alike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=1,lamb=0.001,nug=1,lamg=0.001)
QE1000a=MCMC_IncGSE(Abik,1000,Elike_GSEI,120,0.1,0.087,3,1000,10000,0.01,nub=1,lamb=0.001,nug=1,lamg=0.001)

#write.csv(QB1000a,row.names=FALSE,file="MCMC/Abak_1000a_BD.csv")
#write.csv(QA1000a,row.names=FALSE,file="MCMC/Abak_1000a_approx.csv")
#write.csv(QE1000a,row.names=FALSE,file="MCMC/Abak_1000a_exact.csv")

#QA1000=read.csv(file="MCMC/Abak_1000a_approx.csv")
#QB1000=read.csv(file="MCMC/Abak_1000a_BD.csv")
#QE1000=read.csv(file="MCMC/Abak_1000a_exact.csv")

#QA1000a=as.matrix(QA1000a)
#QB1000a=as.matrix(QB1000a)
#QE1000a=as.matrix(QE1000a)

Res1000a=matrix(0,ncol=5,nrow=3)

colnames(Res1000a)=c("Mean alpha","Mean mu","sd alpha","sd mu","cor(alpha,mu)")

Res1000a[1,1:2]=colMeans(QA1000a)
Res1000a[1,3]=sd(QA1000a[,1])
Res1000a[1,4]=sd(QA1000a[,2])
Res1000a[1,5]=cor(QA1000a[,1],QA1000a[,2])
Res1000a[2,1:2]=colMeans(QB1000a)
Res1000a[2,3]=sd(QB1000a[,1])
Res1000a[2,4]=sd(QB1000a[,2])
Res1000a[2,5]=cor(QB1000a[,1],QB1000a[,2])
Res1000a[3,1:2]=colMeans(QE1000a)
Res1000a[3,3]=sd(QE1000a[,1])
Res1000a[3,4]=sd(QE1000a[,2])
Res1000a[3,5]=cor(QE1000a[,1],QE1000a[,2])

#> Res1000a
#Mean alpha   Mean mu   sd alpha      sd mu cor(alpha,mu)
#[1,]  0.1233623 0.1056000 0.03744951 0.03164268     0.5737152
#[2,]  0.1246990 0.1060180 0.03802188 0.03198264     0.5773341
#[3,]  0.1174258 0.1061479 0.03455501 0.03090267     0.6098327

Abik=c(13,7,2,3,0,0,1,4,5,3,2,0,2,0,5,3,1,4,0,1,1,1,2,0,1,5,0,5,5)

RB47=MCMC_IncGSE(Abik,47,like_GSEI,120,0.1,0.087,3,5000,50000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
RA47=MCMC_IncGSE(Abik,47,alike_GSEI,120,0.1,0.087,3,5000,50000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
RE47=MCMC_IncGSE(Abik,47,Elike_GSEI,120,0.1,0.087,3,5000,50000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)

#write.csv(RB47,row.names=FALSE,file="MCMC/Abak_47L_BD.csv")
#write.csv(RA47,row.names=FALSE,file="MCMC/Abak_47L_approx.csv")
#write.csv(RE47,row.names=FALSE,file="MCMC/Abak_47L_exact.csv")

RB90=MCMC_IncGSE(Abik,90,like_GSEI,120,0.1,0.087,3,5000,50000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
RA90=MCMC_IncGSE(Abik,90,alike_GSEI,120,0.1,0.087,3,5000,50000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
RE90=MCMC_IncGSE(Abik,90,Elike_GSEI,120,0.1,0.087,3,5000,50000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)

#write.csv(RB90,row.names=FALSE,file="MCMC/Abak_90L_BD.csv")
#write.csv(RA90,row.names=FALSE,file="MCMC/Abak_90L_approx.csv")
#write.csv(RE90,row.names=FALSE,file="MCMC/Abak_90L_exact.csv")

RB46b=MCMC_IncGSE(Abik,46.99,like_GSEI,120,0.1,0.087,3,5000,50000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
RA46b=MCMC_IncGSE(Abik,46.99,alike_GSEI,120,0.1,0.087,3,5000,50000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)
RE46b=MCMC_IncGSE(Abik,46.99,Elike_GSEI,120,0.1,0.087,3,5000,50000,0.01,nub=10,lamb=500/6,nug=10,lamg=100)

#write.csv(RB46b,row.names=FALSE,file="MCMC/Abak_46bL_BD.csv")
#write.csv(RA46b,row.names=FALSE,file="MCMC/Abak_46bL_approx.csv")
#write.csv(RE46b,row.names=FALSE,file="MCMC/Abak_46bL_exact.csv")

QA46b=RA46b
QB46b=RB46b
QE46b=RE46b

Res46b=matrix(0,ncol=5,nrow=3)

colnames(Res46b)=c("Mean alpha","Mean mu","sd alpha","sd mu","cor(alpha,mu)")

Res46b[1,1:2]=colMeans(QA46b)
Res46b[1,3]=sd(QA46b[,1])
Res46b[1,4]=sd(QA46b[,2])
Res46b[1,5]=cor(QA46b[,1],QA46b[,2])
Res46b[2,1:2]=colMeans(QB46b)
Res46b[2,3]=sd(QB46b[,1])
Res46b[2,4]=sd(QB46b[,2])
Res46b[2,5]=cor(QB46b[,1],QB46b[,2])
Res46b[3,1:2]=colMeans(QE46b)
Res46b[3,3]=sd(QE46b[,1])
Res46b[3,4]=sd(QE46b[,2])
Res46b[3,5]=cor(QE46b[,1],QE46b[,2])

#Mean alpha    Mean mu   sd alpha      sd mu cor(alpha,mu)
#[1,]  0.1149950 0.08458420 0.02759912 0.02517710     0.2511495
#[2,]  0.1142582 0.08372964 0.02649531 0.02432794     0.2424608
#[3,]  0.1106943 0.08405296 0.02575210 0.02467165     0.2640295
