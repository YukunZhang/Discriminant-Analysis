set.seed(2014)
library(clusterGeneration)
library(psych)
library(mvtnorm)
library(compositions)   #for multivariate lognormal
library(rrcov)
library(evd)
library(fMultivar)
library(Matrix)
#setting=read.csv(file="simulation-settings.csv",header=TRUE)
#setting=read.csv(file="/files3/home/yukzhang/simulation-settings2.csv",header=TRUE)
setting=commandArgs(trailing=TRUE)
# setting=list(c(2,7,'CS',0.3,32,48,1,9,1,"norm",1))
setting=strsplit(setting[[1]],',')
tryCatch({
  p=as.numeric(setting[[1]][1])  # number of time points
  q=as.numeric(setting[[1]][2])# number of measurements
  pi=3.14159265
  converge=0.00001
  nsim=1000 #number of simulations
  claserr=matrix(0,nsim,2)
  one=matrix(rep(1,p),p,1)   #increasing structured mean
  one1=matrix(rep(1,q),q,1)  
  c1=diag(c(0,rep(1,p-2),0))
  c11=diag(c(0,rep(1,q-2),0))
  
  c2=matrix(0,nrow=p,ncol=p)
  c21=matrix(0,nrow=q,ncol=q)
  c2[row(c2)==col(c2)]=0
  c21[row(c21)==col(c21)]=0
  c2[abs(row(c2)-col(c2))==1]=1
  c21[abs(row(c21)-col(c21))==1]=1
  rou1=as.numeric(setting[[1]][4]) #correlation of repeated measures for population 1  rou can be 0.5 or 0.8
  # rou2=0.5
  if(setting[[1]][3]=="AR")
  {  V1=(1-rou1^2)*solve((diag(p)+rou1^2*c1-rou1*c2))  }  else #AR(1)correlation structure
    if(setting[[1]][3]=="CS")
    {  V1=(1-rou1)*diag(p)+rou1*one%*%t(one)}      #CS correlation structure
  # V2=(1-rou2^2)*solve((diag(p)+rou2^2*c1-rou2*c2))    #AR(1)correlation structure
  if (setting[[1]][11]=="CS"){sigma1=(1-rou1)*diag(q)+rou1*one1%*%t(one1)
  }else 
    if (setting[[1]][11]=="AR"){sigma1=(1-rou1^2)*solve((diag(q)+rou1^2*c11-rou1*c21)) }
  
  sigma1=sigma1*60
  omega1=V1%x%sigma1
  omega1log=V1%x%log(sigma1)
  #omega1=round(omega1,2)
  
  #sigma2=genPositiveDefMat("eigen",dim=q)
  #omega2=V2%x%sigma2$Sigma
  omega2=omega1/as.numeric(setting[[1]][7])
  omega2log=omega1log/as.numeric(setting[[1]][7])
  
  #omega2=round(omega2,2)
  
  # mu1=matrix(c(1,1.1,1.2),p,1)%x%del1
  # mu2=matrix(c(1,1.1,1.2),p,1)%x%del2
  
  if (q==3) {del2=c(20,20,20)
             if(as.numeric(setting[[1]][8])==1)
             {    
               del1=c(25,25,25) 
             }else
               if (as.numeric(setting[[1]][8])==2)
               {   del1=c(30,30,30) }else
                 if (  as.numeric(setting[[1]][8])==3)
                 {   del1=c(23,25,27) }else
                   if (as.numeric(setting[[1]][8])==4)
                   {del1=c(27,29,31)}else
                     if( as.numeric(setting[[1]][8])==5)
                     {del1=c(25,23,25)}else
                       if(as.numeric(setting[[1]][8])==6)
                       {del1=c(30,25,30)}
  }else 
    if (q==7){del2=rep(20,7)
              if(as.numeric(setting[[1]][8])==1)
              {    
                del1=rep(25,7)
              }else
                if (as.numeric(setting[[1]][8])==2)
                {   del1=rep(30,7) }else
                  if (  as.numeric(setting[[1]][8])==3)
                  {   del1=c(25,26,27,28,29,30,31) }else
                    if (  as.numeric(setting[[1]][8])==4)
                    {del1=c(30,31,32,33,34,35,36)}else
                      if(as.numeric(setting[[1]][8])==5)
                      {del1=c(25,24,23,22,23,24,25)}else
                        if(as.numeric(setting[[1]][8])==6)
                        {del1=c(30,29,28,27,28,29,30)}
              
    }  
  
  if(p==2&as.numeric(setting[[1]][9])==1)
  {
    Wvm=matrix(c(1,1),p,1) #define within-variable means
  }
  if(p==2&as.numeric(setting[[1]][9])==2)
  {
    Wvm=matrix(c(1,0.5),p,1) #define within-variable means
  }
  if(p==3&as.numeric(setting[[1]][9])==1)
  {
    Wvm=matrix(c(1,1,1),p,1) #define within-variable means
  }
  if(p==3&as.numeric(setting[[1]][9])==2)
  {
    Wvm=matrix(c(1,0.5,1),p,1) #define within-variable means
  }
  
  mu1=Wvm%x%del1
  mu2=Wvm%x%del2
  x1=c(mu1)
  x2=c(mu2)
  Math.cbrt <- function(x) {
    sign(x) * abs(x)^(1/3)
  }