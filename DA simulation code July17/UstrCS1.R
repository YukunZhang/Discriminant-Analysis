if(p==5){
  for(i in 1:q)
  {
    
    data1[,(1+(i-1)*p):(i*p)]=data[,c(i,(i+q),(i+2*q),(i+3*q),(i+4*q))]   #measurement 1
    gp1bar[,(1+(i-1)*p):(i*p)]=t(matrix(apply(data1[,(1+(i-1)*p):(i*p)],2,sum)/n))
    gp1a=sweep(data1[,(1+(i-1)*p):(i*p)],MARGIN=2,gp1bar[,(1+(i-1)*p):(i*p)],FUN="-")
    gp1ass=apply(gp1a^2,2,sum)
    gp1an=sweep(gp1a,MARGIN=2,sqrt(gp1ass),FUN="/")
    vi1=vi1+t(gp1an)%*%gp1an
  }
}
if(p==3){ 
  for(i in 1:q)
  {
    
    data1[,(1+(i-1)*p):(i*p)]=data[,c(i,(i+q),(i+2*q))]   #measurement 1
    gp1bar[,(1+(i-1)*p):(i*p)]=t(matrix(apply(data1[,(1+(i-1)*p):(i*p)],2,sum)/n))
    gp1a=sweep(data1[,(1+(i-1)*p):(i*p)],MARGIN=2,gp1bar[,(1+(i-1)*p):(i*p)],FUN="-")
    gp1ass=apply(gp1a^2,2,sum)
    gp1an=sweep(gp1a,MARGIN=2,sqrt(gp1ass),FUN="/")
    vi1=vi1+t(gp1an)%*%gp1an
  }
}
vi=vi1/q

j0=matrix(0,q,q)
if (p==3)
{
  
    j1=cbind(diag(q),j0,j0)
    j2=cbind(j0,diag(q),j0)
    j3=cbind(j0,j0,diag(q))
  }else{
    j1=cbind(diag(q),j0,j0,j0,j0)
    j2=cbind(j0,diag(q),j0,j0,j0)
    j3=cbind(j0,j0,diag(q),j0,j0)
    j4=cbind(j0,j0,j0,diag(q),j0)
    j5=cbind(j0,j0,j0,j0,diag(q))
  }


w1=matrix(0,p*q,p*q)
for (i in 1:n1)
{
  subb=a[i,]                                                                                                                                                                                                                                                                       
  w1=w1+t(subb-amu)%*%(subb-amu)
}
w2=matrix(0,p*q,p*q)
for (i in 1:n2)
{
  subb=b[i,]
  w2=w2+t(subb-bmu)%*%(subb-bmu)
}
w=w1+w2
#### Calculating the initial estimate Ve of V. rest is the common estimate of rho in both the classes####
rest=(sum(vi)-tr(vi))/(p*(p-1))
ve=(1-rest)*diag(p)+rest*one%*%t(one)
#### Calculating the mle's mlv and mlsig of V and sigma respectivaly in population 2####
mlsige=matrix(0,q,q)
iter=0
maxab=1
while ((maxab>converge)&(iter<100))
{
  ive=solve(ve)
  try=diag(1,p)
  try1=t(try)%*%ive%*%try
  sigma2e=matrix(0,q,q)
  for (j in 1:n1)
  {
    subb=a[j,]
    subba=subb-amu
    
    sub1=j1%*%t(subba)
    sub2=j2%*%t(subba)
    sub3=j3%*%t(subba)
    if (p==3){
     
      sigma2e=sigma2e+try1[1,1]*(sub1%*%t(sub1))+try1[2,1]*(sub1%*%t(sub2))+try1[3,1]*(sub1%*%t(sub3))+try1[1,2]*(sub2%*%t(sub1))+try1[2,2]*(sub2%*%t(sub2))+try1[3,2]*(sub2%*%t(sub3))+try1[1,3]*(sub3%*%t(sub1))+try1[2,3]*(sub3%*%t(sub2))+try1[3,3]*(sub3%*%t(sub3))
    }else
    {     sub4=j4%*%t(subba)  
          sub5=j5%*%t(subba)
          sigma2e=sigma2e+try1[1,1]*(sub1%*%t(sub1))+try1[2,1]*(sub1%*%t(sub2))+try1[3,1]*(sub1%*%t(sub3))+try1[4,1]*(sub1%*%t(sub4))+try1[5,1]*(sub1%*%t(sub5))+try1[1,2]*(sub2%*%t(sub1))+try1[2,2]*(sub2%*%t(sub2))+try1[3,2]*(sub2%*%t(sub3))+try1[4,2]*(sub2%*%t(sub4))+try1[5,2]*(sub2%*%t(sub5))+try1[1,3]*(sub3%*%t(sub1))+try1[2,3]*(sub3%*%t(sub2))+try1[3,3]*(sub3%*%t(sub3))+try1[4,3]*(sub3%*%t(sub4))+try1[5,3]*(sub3%*%t(sub5))+try1[1,4]*(sub4%*%t(sub1))+try1[2,4]*(sub4%*%t(sub2))+try1[3,4]*(sub4%*%t(sub3))+try1[4,4]*(sub4%*%t(sub4))+try1[5,4]*(sub4%*%t(sub5))+try1[1,5]*(sub5%*%t(sub1))+try1[2,5]*(sub5%*%t(sub2))+try1[3,5]*(sub5%*%t(sub3))+try1[4,5]*(sub5%*%t(sub4))+try1[5,5]*(sub5%*%t(sub5))
    }
  }
  for (j in 1:n2)
  {
    subb=b[j,]
    subba=subb-bmu
    
    sub1=j1%*%t(subba)
    sub2=j2%*%t(subba)
    sub3=j3%*%t(subba)
    
    if (p==3){
      #   sub3=j3%*%t(subba)
      #   sigma2e=sigma2e+try1[1,1]*(sub1%*%t(sub1))+try1[2,1]*(sub1%*%t(sub2))+try1[3,1]*(sub1%*%t(sub3))+try1[1,2]*(sub2%*%t(sub1))+try1[2,2]*(sub2%*%t(sub2))+try1[3,2]*(sub2%*%t(sub3))+try1[1,3]*(sub3%*%t(sub1))+try1[2,3]*(sub3%*%t(sub2))+try1[3,3]*(sub3%*%t(sub3))
    }else{ sub4=j4%*%t(subba)  
           sub5=j5%*%t(subba)
           sigma2e=sigma2e+try1[1,1]*(sub1%*%t(sub1))+try1[2,1]*(sub1%*%t(sub2))+try1[3,1]*(sub1%*%t(sub3))+try1[4,1]*(sub1%*%t(sub4))+try1[5,1]*(sub1%*%t(sub5))+try1[1,2]*(sub2%*%t(sub1))+try1[2,2]*(sub2%*%t(sub2))+try1[3,2]*(sub2%*%t(sub3))+try1[4,2]*(sub2%*%t(sub4))+try1[5,2]*(sub2%*%t(sub5))+try1[1,3]*(sub3%*%t(sub1))+try1[2,3]*(sub3%*%t(sub2))+try1[3,3]*(sub3%*%t(sub3))+try1[4,3]*(sub3%*%t(sub4))+try1[5,3]*(sub3%*%t(sub5))+try1[1,4]*(sub4%*%t(sub1))+try1[2,4]*(sub4%*%t(sub2))+try1[3,4]*(sub4%*%t(sub3))+try1[4,4]*(sub4%*%t(sub4))+try1[5,4]*(sub4%*%t(sub5))+try1[1,5]*(sub5%*%t(sub1))+try1[2,5]*(sub5%*%t(sub2))+try1[3,5]*(sub5%*%t(sub3))+try1[4,5]*(sub5%*%t(sub4))+try1[5,5]*(sub5%*%t(sub5))
    }
  }
  mlsig=sigma2e/(n*p)
  absig=abs(tr(mlsige-mlsig))
  mlsige=mlsig
  imlsig=solve(mlsig)
  k3=tr((diag(p)%x%imlsig)%*%w)
  k4=tr(((one%*%t(one))%x%imlsig)%*%w)
  #kk3=tr((diag(p)%x%imlsig)%*%w)
  #mk3=tr((c1%x%imlsig)%*%w)
  #nk3=tr((c2%x%imlsig)%*%w)
  ####solving the cubic equation ####
  s=p-1
  # pp=(nk3)/(-2*n*s*q)
  # qq=(2*n*s*q-2*kk3-2*mk3)/(-2*n*s*q)
  # rr=nk3/(-2*n*s*q)
  ko=n*q*s*p
  pp=(ko-s*ko+k3*s^2-s*k4)/(s*ko)
  qq=(2*s*k3-ko)/(s*ko)
  rr=(k3-k4)/(s*ko)
  aa=(1/3)*(3*qq-pp^2)
  bb=(1/27)*(2*pp^3 -9*pp*qq +27*rr)
  discrim=(bb^2)/4+(aa^3)/27
  if (discrim>0){
    s1=((bb^2)/4+(aa^3)/27)^(1/2)
    s2=Math.cbrt(-bb/2+s1)
    ar=-bb/2-s1
    if(ar<0){ar1=-ar
             ar2=Math.cbrt(ar1)
             s3=-ar2}else{
               s3=Math.cbrt(ar)}
    
    ro2=s2+s3-pp/3
    mlv=(1-ro2)*diag(p)+ro2*one%*%t(one)
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
    mlv=(1-ro2)*diag(p)+ro2*one%*%t(one)
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
