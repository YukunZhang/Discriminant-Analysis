####trimmed####

####conventional####
rm(list = ls(all = TRUE))
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
setting=strsplit(setting[[1]],',')
  tryCatch({
    p=as.numeric(setting[[1]][1])  # number of time points
    q=as.numeric(setting[[1]][2])# number of measurements
    pi=3.14159265
    converge=0.00001
    nsim=1000 #number of simulations
    claserr=matrix(0,nsim,2)
    one=matrix(rep(1,p),p,1)   #increasing structured mean
    c1=diag(c(0,rep(1,p-2),0))
    c2=matrix(0,nrow=p,ncol=p)
    c2[row(c2)==col(c2)]=0
    c2[abs(row(c2)-col(c2))==1]=1
    rou1=as.numeric(setting[[1]][4]) #correlation of repeated measures for population 1  rou can be 0.5 or 0.8
    # rou2=0.5
    if(setting[[1]][3]=="AR")
    {  V1=(1-rou1^2)*solve((diag(p)+rou1^2*c1-rou1*c2))  }  else #AR(1)correlation structure
      if(setting[[1]][3]=="CS")
      {  V1=(1-rou1)*diag(p)+rou1*one%*%t(one)}      #CS correlation structure
    # V2=(1-rou2^2)*solve((diag(p)+rou2^2*c1-rou2*c2))    #AR(1)correlation structure
    if (q==7&rou1==0.7){sigma1=matrix(c(1,0.65,0.66,0.7,0.72,0.65,0.75,
                                        0.65,1,0.7,0.7,0.67,0.72,0.74,
                                        0.66,0.7,1,0.71,0.75,0.63,0.74,
                                        0.7,0.7,0.71,1,0.7,0.65,0.73,
                                        0.72,0.67,0.75,0.7,1,0.68,0.74,
                                        0.65,0.72,0.63,0.65,0.68,1,0.75,
                                        0.75,0.74,0.74,0.73,0.74,0.75,1),nrow=7)
    }else 
      if (q==3&rou1==0.7){sigma1=matrix(c(1,0.65,0.66,
                                          0.65,1,0.7,
                                          0.66,0.7,1),nrow=3)
      }else
        if(q==3 &rou1==0.3){sigma1=matrix(c(1,0.15,0.3,
                                            0.15,1,0.45,
                                            0.3,0.45,1),nrow=3)} else
                                              if(q==7&rou1==0.3){sigma1=matrix(c(1,0.35,0.2,0.2,0.3,0.25,0.3,
                                                                                 0.35,1,0.3,0.3,0.37,0.32,0.34,
                                                                                 0.2,0.3,1,0.31,0.35,0.32,0.34,
                                                                                 0.2,0.3,0.31,1,0.3,0.25,0.33,
                                                                                 0.3,0.37,0.35,0.3,1,0.28,0.34,
                                                                                 0.25,0.32,0.32,0.25,0.28,1,0.35,
                                                                                 0.3,0.34,0.34,0.33,0.34,0.35,1
                                              ),nrow=7)}
    
    
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
    } else 
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
    for (t in 1:nsim)
    {
      n1=as.numeric(setting[[1]][5])
      n2=as.numeric(setting[[1]][6])
      n=n1+n2
      n1sum=n1
      n2sum=n2
      nsum=n1+n2
      
      if (setting[[1]][10]=="norm"){
        a=mvrnorm(n1, mu1, omega1)
        b=mvrnorm(n2, mu2, omega2)
      }else
        if (setting[[1]][10]=="T")
        { a=rmvt(n1, sigma = omega1, df = (n1-p),delta=mu1)
          b=rmvt(n2, sigma = omega2, df = (n1-p),delta=mu2)
        }else
          if (setting[[1]][10]=="lognormal")
          {   a=rlnorm.rplus(n1,c(log(mu1)),forceSymmetric(omega1log) ) #log normal
              b=rlnorm.rplus(n2,c(log(mu2)),forceSymmetric(omega2log))
          } else 
            if (setting[[1]][10]=="cauchy")
            {
              
              a=sn::rmsc(n1, x1,forceSymmetric(omega1),c(rep(4,p*q)), dp=NULL)
              b=sn::rmsc(n2, x2, forceSymmetric(omega2),c(rep(4,p*q)), dp=NULL)
            }
      
      
      data=rbind(a,b)
      datasum=data
      
      amu=t(matrix(apply(a,2,sum)/(n1)))
      bmu=t(matrix(apply(b,2,sum)/(n2)))
      asummu=amu
      bssummu=bmu
      asum=a
      bsum=b
      
      d1=matrix(amu,q,p)
      d2=matrix(bmu,q,p)
      
      ####to calculate amu and bmu####
      
      
      
      c1=diag(c(0,rep(1,p-2),0))
      c2=matrix(0,nrow=p,ncol=p)
      one=matrix(rep(1,p),p,1)   #increasing structured mean
      j0=matrix(0,q,q)
      if (p==2)
      {
        j1=cbind(diag(q),j0)
        j2=cbind(j0,diag(q))}else{
          j1=cbind(diag(q),j0,j0)
          j2=cbind(j0,j0,diag(q))
          j3=cbind(j0,j0,diag(q))
        }
      
      ####trimmed####
      S=matrix(0,p*q,p*q)
      S1=matrix(0,p*q,p*q)
      i=0
      for (i in 1:n1)
      {S=t(a[i,]-amu)%*%(a[i,]-amu)
       # print(S)
       S1=S+S1
      }
      per=0.05*n1
      S1=S1/(n1-1)
      i=0
      loglik=matrix(0,n1,1)
      for (i in 1:n1)
      {
        loglik[i,]=-p*q/2*log(2*pi)-q/2*log(det(S1))-0.5*((a[i,]-amu)%*%solve(S1)%*%t(a[i,]-amu))
      }
      delete=c(rank(loglik)[(n1-per):n1],rank(loglik)[1:per])
      #delete=c(rank(loglik)[(n1-19):n1])
      
      a=a[-delete,]
      n1=nrow(a)
      amu=t(matrix(apply(a,2,sum)/(n1)))
      
      w1=matrix(0,p*q,p*q)
      for (i in 1:n1)
      {
        subb=a[i,]                                                                                                                                                                                                                                                                       
        w1=w1+t(subb-amu)%*%(subb-amu)
      }
      
      
      a1=matrix(0,n1,p*q)
      gp1bar=matrix(0,1,p*q)
      vi1=0
      
      ####when p=2####
      if(p==2){
        for(i in 1:q)
        {
          
          a1[,(1+(i-1)*p):(i*p)]=a[,c(i,(i+q))]   #measurement 1
          gp1bar[,(1+(i-1)*p):(i*p)]=t(matrix(apply(a1[,(1+(i-1)*p):(i*p)],2,sum)/n))
          gp1a=sweep(a1[,(1+(i-1)*p):(i*p)],MARGIN=2,gp1bar[,(1+(i-1)*p):(i*p)],FUN="-")
          gp1ass=apply(gp1a^2,2,sum)
          gp1an=sweep(gp1a,MARGIN=2,sqrt(gp1ass),FUN="/")
          vi1=vi1+t(gp1an)%*%gp1an
        }
      }
      if(p==3){ 
        for(i in 1:q)
        {
          
          a1[,(1+(i-1)*p):(i*p)]=a[,c(i,(i+q),(i+2*q))]   #measurement 1
          gp1bar[,(1+(i-1)*p):(i*p)]=t(matrix(apply(a1[,(1+(i-1)*p):(i*p)],2,sum)/n))
          gp1a=sweep(a1[,(1+(i-1)*p):(i*p)],MARGIN=2,gp1bar[,(1+(i-1)*p):(i*p)],FUN="-")
          gp1ass=apply(gp1a^2,2,sum)
          gp1an=sweep(gp1a,MARGIN=2,sqrt(gp1ass),FUN="/")
          vi1=vi1+t(gp1an)%*%gp1an
        }
      }
      vi=vi1/q
      
      
      #### Calculating the initial estimate Ve of V. rest is the common estimate of rho in both the classes####
      rest=(sum(vi)-tr(vi))/(p*(p-1))
      ve=((1-rest^2)^(-1))*(diag(p)+rest^2*c1-rest*c2)
      #print (ve)
      #### Calculating the mle's mlv and mlsig of V and sigma respectivaly in population 2####
      mlsige=matrix(0,q,q)
      iter=0
      maxab=1
      
      while ((maxab>converge)&(iter<100))
      {
        ive=solve(ve)
        try=diag(1,p)
        try1=t(try)%*%ive%*%try
        
        
        if (p==3)
        {delta1=(sum(try1[1,])*d1[,1]+sum(try1[2,])*d1[,2]+sum(try1[3,])*d1[,3])/sum(try1)
        }else
        {delta1=(sum(try1[1,])*d1[,1]+sum(try1[2,])*d1[,2])/sum(try1)
        }
        
        sigma2e=matrix(0,q,q)
        
        #  detw=det(w)
        k1=(t(amu)-one%x%delta1)%*%t(t(amu)-one%x%delta1)
        
        s1=w1+n1*k1
        sigma2e=matrix(0,q,q)
        for (j in 1:n1)
        {
          subb=a[j,]
          subba=subb-amu
          
          sub1=j1%*%t(subba)
          sub2=j2%*%t(subba)
          if (p==3){
            sub3=j3%*%t(subba)
            sigma2e=sigma2e+try1[1,1]*(sub1%*%t(sub1))+try1[2,1]*(sub1%*%t(sub2))+try1[3,1]*(sub1%*%t(sub3))+try1[1,2]*(sub2%*%t(sub1))+try1[2,2]*(sub2%*%t(sub2))+try1[3,2]*(sub2%*%t(sub3))+try1[1,3]*(sub3%*%t(sub1))+try1[2,3]*(sub3%*%t(sub2))+try1[3,3]*(sub3%*%t(sub3))
          }else
          {          sigma2e=sigma2e+try1[1,1]*(sub1%*%t(sub1))+try1[2,1]*(sub1%*%t(sub2))+try1[1,2]*(sub2%*%t(sub1))+try1[2,2]*(sub2%*%t(sub2))
          }
        }
        
        mlsig=sigma2e/(n1*p)
        absig=abs(tr(mlsige-mlsig))
        mlsige=mlsig
        imlsig=solve(mlsig)
        
        f1=tr((diag(p)%x%imlsig)%*%s1)
        
        f2=tr((c1%x%imlsig)%*%s1)
        f3=tr((c2%x%imlsig)%*%s1)
        ####solving the cubic equation ####
        s=p-1
        pp=(f3)/(-2*n1*s*q)
        qq=(2*n1*s*q-2*f1-2*f2)/(-2*n1*s*q)
        rr=f3/(-2*n1*s*q)
        
        aa=(1/3)*(3*qq-pp^2)
        bb=(1/27)*(2*pp^3 -9*pp*qq +27*rr)
        discrim=(bb^2)/4+(aa^3)/27
        if (discrim>0){
          ss1=((bb^2)/4+(aa^3)/27)^(1/2)
          ss2=Math.cbrt(-bb/2+ss1)
          ar=-bb/2-ss1
          if(ar<0){ar1=-ar
                   ar2=Math.cbrt(ar1)
                   ss3=-ar2}else{
                     ss3=Math.cbrt(ar)}
          
          ro2=ss2+ss3-pp/3
          mlv=(1-ro2^2)*solve((diag(p)+ro2^2*c1-ro2*c2)) 
          ve=mlv}
        else{
          dm=sqrt(-discrim)
          r=sqrt((bb^2)/4+dm^2)
          th=atan(-2*dm/bb)
          pplq=2*(Math.cbrt(r))*cos(th/3)
          pmiq=2*(Math.cbrt(r))*sin(th/3)
          rt1=pplq-pp/3
          rt2=-0.5*pplq-pp/3-0.5*pmiq*sqrt(3)
          rt3=-0.5*pplq-pp/3+0.5*pmiq*sqrt(3)
          if(rt1>1){rt1=0}
          if(rt2>1){rt2=0}
          if(rt3>1){rt3=0}
          ro2=max(rt1,rt2,rt3)
          mlv=(1-ro2^2)*solve((diag(p)+ro2^2*c1-ro2*c2)) 
          ve=mlv
        }
        iter=iter+1
        abr1=abs(rest-ro2)
        rest=ro2
        maxab=max(absig,abr1)
      }
      
      aro2=ro2
      amlsig=mlsig
      av=mlv
      
      c1=diag(c(0,rep(1,p-2),0))
      c2=matrix(0,nrow=p,ncol=p)
      one=matrix(1,p,1)
      j0=matrix(0,q,q)
      if (p==2)
      {
        j1=cbind(diag(q),j0)
        j2=cbind(j0,diag(q))}else{
          j1=cbind(diag(q),j0,j0)
          j2=cbind(j0,j0,diag(q))
          j3=cbind(j0,j0,diag(q))
        }
      #gpbar=t(matrix(apply(gp,2,sum)/n))
      # gp1bar=t(matrix(apply(gp1,2,sum)/n2))
      # gp2bar=t(matrix(apply(gp2,2,sum)/n2))
      # gp3bar=t(matrix(apply(gp3,2,sum)/n2))
      
      # w=matrix(0,p*q,p*q)
      # for (i in 1:n)
      # {
      #   subb=gp[i,]
      #   w=w+t(subb-gpbar)%*%(subb-gpbar)
      # }
      # detw=det(w)
      # 
      # n1=n2=9
      per=0.05*n2
      S=matrix(0,p*q,p*q)
      S1=matrix(0,p*q,p*q)
      for (i in 1:n2)
      {S=t(b[i,]-bmu)%*%(b[i,]-bmu)
       S1=S+S1
      }
      S1=S1/(n2-1)
      i=0
      loglik=matrix(0,n2,1)
      
      for (i in 1:n2)
      {
        loglik[i,]=-p*q/2*log(2*pi)-1/2*log(det(S1))-0.5*((b[i,]-bmu)%*%solve(S1)%*%t(b[i,]-bmu))
      }
      #delete=c(rank(loglik)[(n2-n2*per+1):n2],rank(loglik)[1:(n2*per)])
      #delete=c(rank(loglik)[1:(n2*per*2)])
      #delete=c(rank(loglik)[(n2-19):n2])
      delete=c(rank(loglik)[(n2-per):n2],rank(loglik)[1:per])
      b=b[-delete,]
      n2=nrow(b)
      bmu=t(matrix(apply(b,2,sum)/(n2)))
      
      w2=matrix(0,p*q,p*q)
      for (i in 1:n2)
      {
        subb=b[i,]
        w2=w2+t(subb-bmu)%*%(subb-bmu)
      }
      b1=matrix(0,n2,p*q)
      gp1bar=matrix(0,1,p*q)
      vi1=0
      ####when p=2####
      if(p==2){
        for(i in 1:q)
        {
          
          b1[,(1+(i-1)*p):(i*p)]=b[,c(i,(i+q))]   #measurement 1
          gp1bar[,(1+(i-1)*p):(i*p)]=t(matrix(apply(b1[,(1+(i-1)*p):(i*p)],2,sum)/n))
          gp1a=sweep(b1[,(1+(i-1)*p):(i*p)],MARGIN=2,gp1bar[,(1+(i-1)*p):(i*p)],FUN="-")
          gp1ass=apply(gp1a^2,2,sum)
          gp1an=sweep(gp1a,MARGIN=2,sqrt(gp1ass),FUN="/")
          vi1=vi1+t(gp1an)%*%gp1an
        }
      }
      if(p==3){ 
        for(i in 1:q)
        {
          
          b1[,(1+(i-1)*p):(i*p)]=b[,c(i,(i+q),(i+2*q))]   #measurement 1
          gp1bar[,(1+(i-1)*p):(i*p)]=t(matrix(apply(b1[,(1+(i-1)*p):(i*p)],2,sum)/n))
          gp1a=sweep(b1[,(1+(i-1)*p):(i*p)],MARGIN=2,gp1bar[,(1+(i-1)*p):(i*p)],FUN="-")
          gp1ass=apply(gp1a^2,2,sum)
          gp1an=sweep(gp1a,MARGIN=2,sqrt(gp1ass),FUN="/")
          vi1=vi1+t(gp1an)%*%gp1an
        }
      }
      vi=vi1/q
      
      
      #### Calculating the initial estimate Ve of V. rest is the common estimate of rho in both the classes####
      rest=(sum(vi)-tr(vi))/(p*(p-1))
      ve=((1-rest^2)^(-1))*(diag(p)+rest^2*c1-rest*c2)
      #print (ve)
      #### Calculating the mle's mlv and mlsig of V and sigma respectivaly in population 2####
      mlsige=matrix(0,q,q)
      iter=0
      maxab=1
      
      while ((maxab>converge)&(iter<100))
      {
        ive=solve(ve)
        try=diag(1,p)
        try1=t(try)%*%ive%*%try
        
        
        if (p==3)
        {
          delta2=(sum(try1[1,])*d2[,1]+sum(try1[2,])*d1[,2]+sum(try1[3,])*d2[,3])/sum(try1)
        }else
        {
          delta2=(sum(try1[1,])*d2[,1]+sum(try1[2,])*d1[,2])/sum(try1)
        }
        
        
        sigma2e=matrix(0,q,q)
        
        
        #  detw=det(w)
        k2=(t(bmu)-one%x%delta2)%*%t(t(bmu)-one%x%delta2)
        s2=w2+n2*k2
        sigma2e=matrix(0,q,q)
        
        for (j in 1:n2)
        {
          subb=b[j,]
          subba=subb-bmu
          
          sub1=j1%*%t(subba)
          sub2=j2%*%t(subba)
          if (p==3){
            sub3=j3%*%t(subba)
            sigma2e=sigma2e+try1[1,1]*(sub1%*%t(sub1))+try1[2,1]*(sub1%*%t(sub2))+try1[3,1]*(sub1%*%t(sub3))+try1[1,2]*(sub2%*%t(sub1))+try1[2,2]*(sub2%*%t(sub2))+try1[3,2]*(sub2%*%t(sub3))+try1[1,3]*(sub3%*%t(sub1))+try1[2,3]*(sub3%*%t(sub2))+try1[3,3]*(sub3%*%t(sub3))
          }else{ sigma2e=sigma2e+try1[1,1]*(sub1%*%t(sub1))+try1[2,1]*(sub1%*%t(sub2))+try1[1,2]*(sub2%*%t(sub1))+try1[2,2]*(sub2%*%t(sub2))
          }
        }
        mlsig=sigma2e/(n2*p)
        absig=abs(tr(mlsige-mlsig))
        mlsige=mlsig
        imlsig=solve(mlsig)
        
        g1=tr((diag(p)%x%imlsig)%*%s2)
        
        g2=tr((c1%x%imlsig)%*%s2)
        g3=tr((c2%x%imlsig)%*%s2)
        ####solving the cubic equation ####
        s=p-1
        pp=(g3)/(-2*n2*s*q)
        qq=(2*n2*s*q-2*g1-2*g2)/(-2*n2*s*q)
        rr=g3/(-2*n2*s*q)
        
        aa=(1/3)*(3*qq-pp^2)
        bb=(1/27)*(2*pp^3 -9*pp*qq +27*rr)
        discrim=(bb^2)/4+(aa^3)/27
        if (discrim>0){
          ss1=((bb^2)/4+(aa^3)/27)^(1/2)
          ss2=Math.cbrt(-bb/2+ss1)
          ar=-bb/2-ss1
          if(ar<0){ar1=-ar
                   ar2=Math.cbrt(ar1)
                   ss3=-ar2}else{
                     ss3=Math.cbrt(ar)}
          
          ro2=ss2+ss3-pp/3
          mlv=(1-ro2^2)*solve((diag(p)+ro2^2*c1-ro2*c2)) 
          ve=mlv}
        else{
          dm=sqrt(-discrim)
          r=sqrt((bb^2)/4+dm^2)
          th=atan(-2*dm/bb)
          pplq=2*(Math.cbrt(r))*cos(th/3)
          pmiq=2*(Math.cbrt(r))*sin(th/3)
          rt1=pplq-pp/3
          rt2=-0.5*pplq-pp/3-0.5*pmiq*sqrt(3)
          rt3=-0.5*pplq-pp/3+0.5*pmiq*sqrt(3)
          if(rt1>1){rt1=0}
          if(rt2>1){rt2=0}
          if(rt3>1){rt3=0}
          ro2=max(rt1,rt2,rt3)
          mlv=(1-ro2^2)*solve((diag(p)+ro2^2*c1-ro2*c2)) 
          ve=mlv
        }
        iter=iter+1
        abr1=abs(rest-ro2)
        rest=ro2
        maxab=max(absig,abr1)
      }
      
      bro2=ro2
      bmlsig=mlsig
      bv=mlv
      ####discriminant rule ####
      dis=matrix(0, nsum,1)
      for (i in 1:n)
      {
        
        d1=log(n1sum/nsum)-0.5*log(det(amlsig%x%av))-0.5*(datasum[i,]-amu)%*%solve(amlsig%x%av)%*%t(datasum[i,]-amu)
        d2=log(n2sum/nsum)-0.5*log(det(bmlsig%x%bv))-0.5*(datasum[i,]-bmu)%*%solve(bmlsig%x%bv)%*%t(datasum[i,]-bmu)
        dis[i,]=which.max(c(d1,d2))
      }
      
      
      claserr[t,1]=1-length(which(dis[1:n1sum,]==1))/n1sum
      claserr[t,2]=1-length(which(dis[(n1sum+1):nsum,]==2))/n2sum
    }
    clserr1=mean(claserr[,1])
    sdclserr1=sd(1-claserr[,1])
    
    clserr2=mean(claserr[,2])
    sdclserr2=sd(1-claserr[,2])
    msdclserr=mean(c(sdclserr1,sdclserr2))
    predacu1=1-clserr1
    predacu2=1-clserr2
    overallacu=mean(c(predacu1,predacu2))
    overallacusd=sd(apply(claserr,1,mean))
    out=matrix(0,1,18)
    out[1]=setting[[1]][1]
    out[2]=setting[[1]][2]
    out[3]=setting[[1]][3]
    out[4]=setting[[1]][4]
    out[5]=setting[[1]][5]
    out[6]=setting[[1]][6]
    out[7]=setting[[1]][7]
    out[8]=setting[[1]][8]
    out[9]=setting[[1]][9]
    out[10]=setting[[1]][10]
    out[11]=setting[[1]][11]
     out[12]="trimmed"
     out[13]=predacu1
     out[14]=predacu2
     out[15]=sdclserr1
     out[16]=sdclserr2
     out[17]=overallacu
     out[18]=overallacusd


    file=paste ("str_AR4_trimmed_Row#",as.numeric(setting[[1]][11]), ".csv", sep="")
    write.csv(out,file=file,row.names=FALSE)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
