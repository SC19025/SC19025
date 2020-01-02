## ------------------------------------------------------------------------
library(SC19025)
d1<-faithful$eruptions
d2<-faithful$waiting
Test<-independence(d1,d2,0.05)
Test

## ------------------------------------------------------------------------
d1<-rnorm(100)
d2<-rnorm(100,mean=0.5)
Test1<-symmetry(d1,0.05)
Test1
Test2<-symmetry(d2,0.05)
Test2

## ------------------------------------------------------------------------
n <- 1e4
y <- numeric(n)
k<-0
i<-0
c<-4*exp(1)/3
while(k<n){
  u<-runif(1)
  V<-runif(1)
  g<-2*log(3/(2*V))/3
  if(c*u*3/2*exp(-3/2*g)<g*exp(-g*g/2)){
    k<-k+1
    y[k]<-g
  }
  i<-i+1
}

i
mean(y)
hist(y, prob = TRUE, main = expression(f(x)==y*exp(-y^2/2)))
x<-seq(0,4,.01)
lines(x,x*exp(-x*x/2))

## ------------------------------------------------------------------------
x<-rnorm(1000)
y<-rnorm(1000,mean=3)
p<-0.75
z<-p*x+(1-p)*y
hist(z, prob = TRUE)
p<-seq(0.1,0.9,.1)
z<-numeric(11)
for(i in 1:9){
  z<-p[i]*x+(1-p[i])*y
  hist(z, prob = TRUE)
}
p<-2
z<-p*x+(1-p)*y
hist(z, prob = TRUE)
p<--1
z<-p*x+(1-p)*y
hist(z, prob = TRUE)

## ------------------------------------------------------------------------
d<-5
n<-10
sigma<-matrix(nrow=d,ncol=d)
for(i in 1:d){
  for(j in 1:d)

      sigma[i,j]<-min(i/j,j/i)
}
sigma
c<-eigen(sigma, symmetric = TRUE)$vectors
c
S<-matrix(0,nrow=d,ncol=d)
for(i in 1:n){
  x<-rnorm(d)
  y<-c%*%x
  S<-S+y%*%t(y)
}
S

## ------------------------------------------------------------------------
m <- 1e3
x<-runif(m,min=0,max=pi/3)
pihat <-mean(sin(x))*pi/3
pihat
abs(0.5-pihat)
m <- 1e4
x<-runif(m,min=0,max=pi/3)
pihat <-mean(sin(x))*pi/3
pihat
abs(0.5-pihat)
m <- 1e5
x<-runif(m,min=0,max=pi/3)
pihat <-mean(sin(x))*pi/3
pihat
abs(0.5-pihat)
m <- 1e6
x<-runif(m,min=0,max=pi/3)
pihat <-mean(sin(x))*pi/3
pihat
abs(0.5-pihat)

## ------------------------------------------------------------------------
n<-1e3
integrate(function(aa) exp(-aa)/(1+aa^2),0,1)
x<-runif(n)
y<-exp(-x)/(1+x^2)
z<-(exp(-x)/(1+x^2)+exp(x-1)/(1+(1-x)^2))*0.5
y2<-exp(x-1)/(1+(1-x)^2)
mean(y)
mean(y2)
mean(z)
cov(y,y2)
var(y)
var(y2)
var(z)
var(z)/var(y)

## ------------------------------------------------------------------------
M <- 10000
k <- 5
N<-M/k
n<-50
Y<-numeric(M)
for(j in 1:k){
  for(i in 1:N){
    u<-runif(1)
    x <- - log(1-(1-exp(-1))*(j-1-u)/5)
    Y[i+(j-1)*N]<-(1-exp(-1))/(1+x^2)
  }
} 
mean(Y)
sd(Y)
Y1<-numeric(100)
Y2<-numeric(100)
for(l in 1:100){
for(j in 1:k){
  for(i in 1:N){
    u<-runif(1,min=(j-1)*0.2,max=j*0.2)
    x <- - log(1 - u * (1 - exp(-1)))
    Y[i+(j-1)*N]<-(1-exp(-1))/(1+x^2)
  }
} 
Y1[l]<-mean(Y)
u<-runif(M)
x <- - log(1 - u * (1 - exp(-1)))
Y2[l]<-mean((1-exp(-1))/(1+x^2))
}
sd(Y1)
sd(Y2)

## ------------------------------------------------------------------------
n <- 20   # example6.4的算法
alpha <- .05
p<-numeric(10000)
for( i in 1:10000){
x <- rnorm(n, mean=0, sd=2)
UCL <- (n-1) * var(x) / qchisq(alpha, df=n-1)
if (UCL>4){
  p[i]=1
}
else{
  p[i]=0
}
}
mean(p)
n<-20
alpha<-0.05
p<-numeric(10000)  # t区间的覆盖概率
q<-numeric(10000)  #对同样样本计算方差置信区间的覆盖概率
for(i in 1:10000){
  x <- rnorm(n, mean=0, sd=1)
  y <- rnorm(n, mean=0, sd=1)
  z<-x^2+y^2     #n=2的卡方分布                           
  min<-mean(z)+sd(z)*qt(alpha/2,n-1)
  max<-mean(z)+sd(z)*qt(1-alpha/2,n-1)    
  if(min<2&max>2){
    p[i]=1
  }
  else{
    p[i]=0
  }
  UCL <- (n-1) * var(z) / qchisq(alpha, df=n-1)
  if(UCL>4){
    q[i]=1
  }
  else{
    q[i]=0
  }
}
mean(p)
mean(q)

## ------------------------------------------------------------------------
n<-100             #样本数量
m<-100             #每个样本的数据个数
alpha1<-0.025
alpha2<-0.05
alpha3<-0.95
alpha4<-0.975
b<-numeric(n)
for(i in 1:n){
  x<-rnorm(m)
  b[i]=sum((x-mean(x))^3)/(var(x))^(3/2)/m
}
mean(b)             #1000个样本偏度均值
var(b)              #1000个样本偏度方差
q1<-quantile(b,alpha1)
q2<-quantile(b,alpha2)
q3<-quantile(b,alpha3)
q4<-quantile(b,alpha4)
c(q1,q2,q3,q4)          #蒙特卡洛分位数
f<-function(x){exp(-(x-mean(b))^2/(2*var(b)))/((2*pi*var(b))^0.5)}      #近似的正态密度函数
sd1<-(alpha1*(1-alpha1)/n/f(q1)^2)^0.5
sd2<-(alpha2*(1-alpha2)/n/f(q2)^2)^0.5
sd3<-(alpha3*(1-alpha3)/n/f(q3)^2)^0.5
sd4<-(alpha4*(1-alpha4)/n/f(q4)^2)^0.5
c(sd1,sd2,sd3,sd4)                          #偏度分位数的方差
p1<-qnorm(alpha1,mean=0,sd=(6/m)^0.5)
p2<-qnorm(alpha2,mean=0,sd=(6/m)^0.5)
p3<-qnorm(alpha3,mean=0,sd=(6/m)^0.5)
p4<-qnorm(alpha4,mean=0,sd=(6/m)^0.5)
c(p1,p2,p3,p4)                              #大样本偏度的分位数
c(abs((q1-p1)/p1),abs((q2-p2)/p2),abs((q3-p3)/p3),abs((q4-p4)/p4))         #蒙特卡洛偏度分位数雨大样本分位数的误差比

## ------------------------------------------------------------------------
alpha <- .1                      #显著性
alpha1<-1                        #beta分布的alpha参数值
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))     #不同epslion值
N <- length(epsilon)
pwr <- numeric(N)
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
e <- epsilon[j]
sktests <- numeric(m)             
for (i in 1:m) { 
x<-numeric(n)
for(k in 1:n){                           #生成样本
  u<-runif(1)
  if(u>epsilon[j]){
    x[k]<-rbeta(1,shape1 =alpha1,shape2 = alpha1)
  }
  else{
    x[k]<-rnorm(1,mean=0.5,sd=10)              #beta分布的均值为0.5
  }
}
sk<-sum((x-mean(x))^3)/n/(var(x))^1.5
sktests[i] <- as.integer(abs(sk) >= cv)         #峰度检验
}
pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",
xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## ------------------------------------------------------------------------
alpha<-0.05        #两样本t检验
p1<-0.651
p2<-0.676
p<-0
n<-10000
ttest<-(p1-p2)/((p1*(1-p1)+p2*(1-p2))/n)^0.5
if(ttest<qt(alpha/2,df=2*n-2)|ttest>qt(1-alpha/2,df=2*n-2)){
    p<-1
}
p
sigma<-(p1-p2)*sqrt(n)/qnorm(alpha/2)          #z检验给定方差的临界值
sigma
p<-sqrt((p1-p2)^2+n*(p1-p2)^2/(qt(alpha/2,df=n-1))^2)*n      #配对t检验差值不为0个数的临界值
p

## ------------------------------------------------------------------------
library(bootstrap)
mec<-scor$mec
vec<-scor$vec
alg<-scor$alg
ana<-scor$ana
sta<-scor$sta
plot(mec,vec,main=cov(mec,vec)/sqrt(var(mec))/sqrt(var(vec)))               #作每一对成绩的散点图，标题为二者相关系数
plot(mec,alg,main=cov(mec,alg)/sqrt(var(mec))/sqrt(var(alg)))
plot(mec,ana,main=cov(mec,ana)/sqrt(var(mec))/sqrt(var(ana)))
plot(mec,sta,main=cov(mec,sta)/sqrt(var(mec))/sqrt(var(sta)))
plot(vec,alg,main=cov(vec,alg)/sqrt(var(vec))/sqrt(var(alg)))
plot(ana,alg,main=cov(ana,alg)/sqrt(var(ana))/sqrt(var(alg)))
plot(sta,alg,main=cov(sta,alg)/sqrt(var(sta))/sqrt(var(alg)))
plot(ana,vec,main=cov(ana,vec)/sqrt(var(ana))/sqrt(var(vec)))
plot(sta,vec,main=cov(sta,vec)/sqrt(var(sta))/sqrt(var(vec)))
plot(ana,sta,main=cov(ana,sta)/sqrt(var(ana))/sqrt(var(sta)))

B<-200
n<-length(mec)
rho45<-rho35<-rho34<-rho12<-numeric(B)
mecb<-vecb<-algb<-anab<-stab<-numeric(n)
for(i in 1:B){                                                     #bootstrap估计
  rd1<-sample(1:n,n,replace=T)
  for(j in 1:n){
    mecb[j]<-mec[rd1[j]]
    vecb[j]<-vec[rd1[j]]
  }
  rho12[i]<-cov(mecb,vecb)/sqrt(var(mecb))/sqrt(var(vecb))
  rd2<-sample(1:n,n,replace=T)
  for(j in 1:n){
    algb[j]<-alg[rd2[j]]
    anab[j]<-ana[rd2[j]]
  }
  rho34[i]<-cov(algb,anab)/sqrt(var(algb))/sqrt(var(anab))
  rd3<-sample(1:n,n,replace=T)
  for(j in 1:n){
    algb[j]<-alg[rd3[j]]
    stab[j]<-sta[rd3[j]]
  }
  rho35[i]<-cov(algb,stab)/sqrt(var(algb))/sqrt(var(stab))
  rd4<-sample(1:n,n,replace=T)
  for(j in 1:n){
    stab[j]<-sta[rd4[j]]
    anab[j]<-ana[rd4[j]]
  }
  rho45[i]<-cov(stab,anab)/sqrt(var(stab))/sqrt(var(anab))
}
m1<-mean(rho12)
m2<-mean(rho34)
m3<-mean(rho35)
m4<-mean(rho45)
c("rho12","rho34","rho35","rho45")
c(mean(rho12),mean(rho34),mean(rho35),mean(rho45))                #bootstrap相关系数估计值
c(cov(mec,vec)/sqrt(var(mec))/sqrt(var(vec)),cov(ana,alg)/sqrt(var(ana))/sqrt(var(alg)),cov(sta,alg)/sqrt(var(sta))/sqrt(var(alg)),cov(ana,sta)/sqrt(var(ana))/sqrt(var(sta)))         #实际数据相关系数
var1<-var2<-var3<-var4<-0
for(i in 1:B){
  var1<-var1+(rho12[i]-m1)^2/(B-1)
  var2<-var2+(rho34[i]-m1)^2/(B-1)
  var3<-var3+(rho35[i]-m1)^2/(B-1)
  var4<-var4+(rho45[i]-m1)^2/(B-1)
}
c(sqrt(var1),sqrt(var2),sqrt(var3),sqrt(var4))                    #bootstrap估计标准差

## ------------------------------------------------------------------------
alpha<-0.1                                         #置信区间的水平
m<-1000                                            #模拟次数
n<-50                                              #每次选取样本数
B<-200                                             #bootstrap循环次数
p1ln<-p1rn<-p2ln<-p2rn<-p1n<-p2n<-0
p1lb<-p1rb<-p2lb<-p2rb<-p1b<-p2b<-0
p1lp<-p1rp<-p2lp<-p2rp<-p1p<-p2p<-0                #分别记录两种分布，三种置信区间左右和总偏离概率
for(i in 1:m){
  skc<-skn<-numeric(B)
  no<-rnorm(n)
  ch<-rchisq(n,df=5)
  nob<-chb<-numeric(n)
  for(j in 1:B){                                  #bootstrap偏度估计
    x<-sample(1:n,n,replace=T)
    for(k in 1:n){
      nob[k]<-no[x[k]]
      chb[k]<-ch[x[k]]
    }
    skn[j]<-sum((nob-mean(nob))^3)/n/var(nob)^1.5
    skc[j]<-sum((chb-mean(chb))^3)/n/var(chb)^1.5
  }
  sk1<-0
  sk2<-1.6^0.5                                      #实际偏度
  norm1l<-mean(skn)+qnorm(alpha/2)*sqrt(var(skn))   #各种方法置信区间左右端点值       
  norm1r<-mean(skn)-qnorm(alpha/2)*sqrt(var(skn))
  perc1l<-quantile(skn,probs=alpha/2)
  perc1r<-quantile(skn,probs=1-alpha/2)
  basic1l<-2*mean(skn)-quantile(skn,probs=1-alpha/2)
  basic1r<-2*mean(skn)-quantile(skn,probs=alpha/2)
  norm2l<-mean(skc)+qnorm(alpha/2)*sqrt(var(skc))
  norm2r<-mean(skc)-qnorm(alpha/2)*sqrt(var(skc))
  perc2l<-quantile(skc,probs=alpha/2)
  perc2r<-quantile(skc,probs=1-alpha/2)
  basic2l<-2*mean(skc)-quantile(skc,probs=1-alpha/2)
  basic2r<-2*mean(skc)-quantile(skc,probs=alpha/2)
  if(sk1<norm1l){                                   #对每次求出的置信区间分类
    p1ln=p1ln+1/m
    p1n=p1n+1/m
  }
  else if(sk1>norm1r){
    p1rn=p1rn+1/m
    p1n=p1n+1/m
  }
  if(sk1<perc1l){
    p1lp=p1lp+1/m
    p1p=p1p+1/m
  }
  else if(sk1>perc1r){
    p1rp=p1rp+1/m
    p1p=p1p+1/m
  }
  if(sk1<basic1l){
    p1lb=p1lb+1/m
    p1b=p1b+1/m
  }
  else if(sk1>basic1r){
    p1rb=p1rb+1/m
    p1b=p1b+1/m
  }
  if(sk2<norm2l){
    p2ln=p2ln+1/m
    p2n=p2n+1/m
  }
  else if(sk2>norm2r){
    p2rn=p2rn+1/m
    p2n=p2n+1/m
  }
  if(sk2<perc2l){
    p2lp=p2lp+1/m
    p2p=p2p+1/m
  }
  else if(sk2>perc2r){
    p2rp=p2rp+1/m
    p2p=p2p+1/m
  }
  if(sk2<basic2l){
    p2lb=p2lb+1/m
    p2b=p2b+1/m
  }
  else if(sk2>basic2r){
    p2rb=p2rb+1/m
    p2b=p2b+1/m
  }
}
c(p1n,p1b,p1p)                               #正态样本三种90%置信区间不包含真实偏度的比例
c(p1ln,p1lb,p1lp)                            #正态样本三种90%置信区间右偏于真实偏度的比例
c(p1rn,p1rb,p1rp)                            #正态样本三种90%置信区间左偏于真实偏度的比例
c(p2n,p2b,p2p)                               #卡方样本三种90%置信区间不包含真实偏度的比例
c(p2ln,p2lb,p2lp)                            #卡方样本三种90%置信区间右偏于真实偏度的比例
c(p2rn,p2rb,p2rp)                            #卡方样本三种90%置信区间左偏于真实偏度的比例
alpha<-0.05                                        #置信区间的水平
m<-1000                                            #模拟次数
n<-50                                              #每次选取样本数
B<-200                                             #bootstrap循环次数
p1ln<-p1rn<-p2ln<-p2rn<-p1n<-p2n<-0
p1lb<-p1rb<-p2lb<-p2rb<-p1b<-p2b<-0
p1lp<-p1rp<-p2lp<-p2rp<-p1p<-p2p<-0                #分别记录两种分布，三种置信区间左右和总偏离概率
for(i in 1:m){
  skc<-skn<-numeric(B)
  no<-rnorm(n)
  ch<-rchisq(n,df=5)
  nob<-chb<-numeric(n)
  for(j in 1:B){                                  #bootstrap偏度估计
    x<-sample(1:n,n,replace=T)
    for(k in 1:n){
      nob[k]<-no[x[k]]
      chb[k]<-ch[x[k]]
    }
    skn[j]<-sum((nob-mean(nob))^3)/n/var(nob)^1.5
    skc[j]<-sum((chb-mean(chb))^3)/n/var(chb)^1.5
  }
  sk1<-0
  sk2<-1.6^0.5                                      #实际偏度
  norm1l<-mean(skn)+qnorm(alpha/2)*sqrt(var(skn))   #各种方法置信区间左右端点值       
  norm1r<-mean(skn)-qnorm(alpha/2)*sqrt(var(skn))
  perc1l<-quantile(skn,probs=alpha/2)
  perc1r<-quantile(skn,probs=1-alpha/2)
  basic1l<-2*mean(skn)-quantile(skn,probs=1-alpha/2)
  basic1r<-2*mean(skn)-quantile(skn,probs=alpha/2)
  norm2l<-mean(skc)+qnorm(alpha/2)*sqrt(var(skc))
  norm2r<-mean(skc)-qnorm(alpha/2)*sqrt(var(skc))
  perc2l<-quantile(skc,probs=alpha/2)
  perc2r<-quantile(skc,probs=1-alpha/2)
  basic2l<-2*mean(skc)-quantile(skc,probs=1-alpha/2)
  basic2r<-2*mean(skc)-quantile(skc,probs=alpha/2)
  if(sk1<norm1l){                                   #对每次求出的置信区间分类
    p1ln=p1ln+1/m
    p1n=p1n+1/m
  }
  else if(sk1>norm1r){
    p1rn=p1rn+1/m
    p1n=p1n+1/m
  }
  if(sk1<perc1l){
    p1lp=p1lp+1/m
    p1p=p1p+1/m
  }
  else if(sk1>perc1r){
    p1rp=p1rp+1/m
    p1p=p1p+1/m
  }
  if(sk1<basic1l){
    p1lb=p1lb+1/m
    p1b=p1b+1/m
  }
  else if(sk1>basic1r){
    p1rb=p1rb+1/m
    p1b=p1b+1/m
  }
  if(sk2<norm2l){
    p2ln=p2ln+1/m
    p2n=p2n+1/m
  }
  else if(sk2>norm2r){
    p2rn=p2rn+1/m
    p2n=p2n+1/m
  }
  if(sk2<perc2l){
    p2lp=p2lp+1/m
    p2p=p2p+1/m
  }
  else if(sk2>perc2r){
    p2rp=p2rp+1/m
    p2p=p2p+1/m
  }
  if(sk2<basic2l){
    p2lb=p2lb+1/m
    p2b=p2b+1/m
  }
  else if(sk2>basic2r){
    p2rb=p2rb+1/m
    p2b=p2b+1/m
  }
}
c(p1n,p1b,p1p)                               #正态样本三种90%置信区间不包含真实偏度的比例
c(p1ln,p1lb,p1lp)                            #正态样本三种90%置信区间右偏于真实偏度的比例
c(p1rn,p1rb,p1rp)                            #正态样本三种90%置信区间左偏于真实偏度的比例
c(p2n,p2b,p2p)                               #卡方样本三种90%置信区间不包含真实偏度的比例
c(p2ln,p2lb,p2lp)                            #卡方样本三种90%置信区间右偏于真实偏度的比例
c(p2rn,p2rb,p2rp)                            #卡方样本三种90%置信区间左偏于真实偏度的比例

## ------------------------------------------------------------------------
library(bootstrap)
n<-88                                             
lambda<-eigen(cov(scor))$values
theta.hat<-lambda[1]/sum(lambda)
theta.hat                                          #原样本参数估计值
theta<-numeric(n)                                  #n个去一参数估计值
for(i in 1:n){
  lambda<-eigen(cov(scor[-i,]))$values
  theta[i]<-lambda[1]/sum(lambda)
}
bias<-(n-1)*(mean(theta)-theta.hat)              
theta.jack<-theta.hat-bias                  
se<-sqrt((n-1)*sum((theta-mean(theta))^2)/n)
theta.jack                                        #jackknife估计值
bias                                              #偏差
se                                                #标准差

## ------------------------------------------------------------------------
library(DAAG)
magnetic<-ironslag$magnetic
chemical<-ironslag$chemical
n<-length(magnetic)
e1<-e2<-e3<-e4<-numeric(n)      
for (i in 1:n) {
y <- magnetic[-i]
x <- chemical[-i]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[i]
e1[i] <- magnetic[i] - yhat1
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[i] +
J2$coef[3] * chemical[i]^2
e2[i] <- magnetic[i] - yhat2
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[i]
yhat3 <- exp(logyhat3)
e3[i] <- magnetic[i] - yhat3
J4 <- lm(y ~ x+I(x^2)+I(x^3))
yhat4 <- J4$coef[1] + J4$coef[2] * chemical[i] +
J4$coef[3] * chemical[i]^2 + J4$coef[4] * chemical[i]^3
e4[i] <- magnetic[i] - yhat4
}
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))              #留一验证残差和

y <- magnetic
x <- chemical
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical
e1 <- magnetic - yhat1
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical +
J2$coef[3] * chemical^2
e2 <- magnetic - yhat2
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical
yhat3 <- exp(logyhat3)
e3 <- magnetic - yhat3
J4 <- lm(y ~ x+I(x^2)+I(x^3))
yhat4 <- J4$coef[1] + J4$coef[2] * chemical +
J4$coef[3] * chemical^2 + J4$coef[4] * chemical^3
e4 <- magnetic - yhat4
c(1-mean(e1^2)/var(y), 1-mean(e2^2)/var(y), 1-mean(e3^2)/var(y), 1-mean(e4^2)/var(y))     #回归判定系数

## ------------------------------------------------------------------------
R<-999
n<-100
num<-5                     #极端值个数阈值为5
T<-numeric(n)
T1<-numeric(1000)
for(i in 1:n){
  x<-rnorm(20)
  y<-rnorm(30)
  test<-numeric(R)
  for(j in 1:R){
    z<-numeric(50)
    for(k in 1:20){
      z[k]<-x[k]
    }
    for(k in 1:30){
      z[k+20]<-y[k]
    }
    z<-sample(z,50,replace=FALSE)
    for(k in 1:20){
      x[k]<-z[k]
    }
    for(k in 1:30){
      y[k]<-z[k+20]
    }
    X <- x - mean(x)
    Y <- y - mean(y)
    outx <- sum(X > max(Y)) + sum(X < min(Y))
    outy <- sum(Y > max(X)) + sum(Y < min(X))
    if(max(outx,outy)>num){
      test[j]<-1
    }
  }
  T[i]<-(1+sum(test))/(R+1)
}
mean(T)                       #使用了排列方法
for(i in 1:100){
  x<-rnorm(20)
  y<-rnorm(30)
  test<-numeric(R)
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  if(max(outx,outy)>num){
      T1[i]<-1
    }
  
}
mean(T1)                      #未使用排列方法
R<-999
n<-100
num<-8                 #极端值个数阈值为8
T<-numeric(n)
T1<-numeric(100)
for(i in 1:n){
  x<-rnorm(20)
  y<-rnorm(30)
  test<-numeric(R)
  for(j in 1:R){
    z<-numeric(50)
    for(k in 1:20){
      z[k]<-x[k]
    }
    for(k in 1:30){
      z[k+20]<-y[k]
    }
    z<-sample(z,50,replace=FALSE)
    for(k in 1:20){
      x[k]<-z[k]
    }
    for(k in 1:30){
      y[k]<-z[k+20]
    }
    X <- x - mean(x)
    Y <- y - mean(y)
    outx <- sum(X > max(Y)) + sum(X < min(Y))
    outy <- sum(Y > max(X)) + sum(Y < min(X))
    if(max(outx,outy)>num){
      test[j]<-1
    }
  }
  T[i]<-(1+sum(test))/(R+1)
}
mean(T)                   #使用了排列方法
for(i in 1:100){
  x<-rnorm(20)
  y<-rnorm(30)
  test<-numeric(R)
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  if(max(outx,outy)>num){
      T1[i]<-1
    }
  
}
mean(T1)                 #未使用排列方法

## ------------------------------------------------------------------------
library(Ball)
library(MASS)
alpha<-0.1
sigma<-matrix(c(1,0,0,1),2,2)
n<-10*1:8                                                 #样本量
powb2<-powb1<-powd1<-powd2<-numeric(8)
DCOR <- function(x, y) {                                  #距离函数
x <- as.matrix(x)
y <- as.matrix(y)
n <- nrow(x)
m <- nrow(y)
if (n != m || n < 2) stop("Sample sizes must agree")
if (! (all(is.finite(c(x, y)))))
stop("Data contains missing or infinite values")
Akl <- function(x) {
d <- as.matrix(dist(x))
m <- rowMeans(d)
M <- mean(d)
a <- sweep(d, 1, m)
b <- sweep(a, 2, m)
return(b + M)
}
A <- Akl(x)
B <- Akl(y)
dCov <- sqrt(mean(A * B))
dVarX <- sqrt(mean(A * A))
dVarY <- sqrt(mean(B * B))
dCor <- sqrt(dCov / sqrt(dVarX * dVarY))
list(dCov=dCov, dCor=dCor, dVarX=dVarX, dVarY=dVarY)
}
for(i in 1:8){
  bt1<-numeric(100)
  dt1<-dt2<-bt2<-numeric(100)
  for(j in 1:100){                                #循环一百次
    
    x<-mvrnorm(n=n[i],rep(0,2),sigma)
    e<-mvrnorm(n=n[i],rep(0,2),sigma)
    y1<-x/4+e
    y2<-x/4*e
    p1<-p2<-0 
    for(k in 1:100){                                #距离方法
      z1<-rbind(x,y1)
      z2<-rbind(x,y2)
      z3<-z1
      z4<-z2
      ran<-sample(1:(2*n[i]),2*n[i],replace=FALSE)
      for(l in 1:(2*n[i])){
        z3[l,]<-z1[ran[l],]
        z4[l,]<-z2[ran[l],]
      }
      X1<-z3[1:n[i],]
      Y1<-z3[1:n[i],]
      X2<-z4[1:n[i],]
      Y2<-z4[1:n[i],]
      p1<-p1+DCOR(X1,Y1)$dCor/100
      p2<-p2+DCOR(X2,Y2)$dCor/100
    }
    dt1[j]<-p1>alpha
    dt2[j]<-p2>alpha
    bt1[j]<-bd.test(x,y1,R=999,seed=j*12345)$p.value>alpha            #ball方法
    bt2[j]<-bd.test(x,y2,R=999,seed=j*12345)$p.value>alpha
  }
  powd1[i]<-mean(dt1)
  powd2[i]<-mean(dt2)
  powb1[i]<-mean(bt1)
  powb2[i]<-mean(bt2)
}                                             #四种功效
powd1
powd2
powb1
powb2                                         
plot(10*1:8,powd1,type="l")
lines(10*1:8,powb1,col="red")
plot(10*1:8,powd2,type="l")
lines(10*1:8,powb2,col="red")

## ------------------------------------------------------------------------
library(GeneralizedHyperbolic)
sigma=c(0.125,0.25,0.5,1,2,4,8)              #the variances we choose
N<-2000                                  #sample size              
x0<-25
rw.Metropolis <- function( sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= exp(-abs(y)) / exp(-abs(x[i-1]))  ){
x[i] <- y
}
else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
}
rw1 <- rw.Metropolis( sigma[1], x0, N)
rw2 <- rw.Metropolis( sigma[2], x0, N)
rw3 <- rw.Metropolis( sigma[3], x0, N)
rw4 <- rw.Metropolis( sigma[4], x0, N)
rw5 <- rw.Metropolis( sigma[5], x0, N)
rw6 <- rw.Metropolis( sigma[6], x0, N)
rw7 <- rw.Metropolis( sigma[7], x0, N)
print(c(rw1$k, rw2$k, rw3$k, rw4$k,rw5$k,rw6$k,rw7$k))      #the acceptance rates

plot(1:N,rw1$x,type="l",xlab="sample",ylab="X",main="sigma=0.125")    #the chains generated
plot(1:N,rw2$x,type="l",xlab="sample",ylab="X",main="sigma=0.25")
plot(1:N,rw3$x,type="l",xlab="sample",ylab="X",main="sigma=0.5")
plot(1:N,rw4$x,type="l",xlab="sample",ylab="X",main="sigma=1")
plot(1:N,rw5$x,type="l",xlab="sample",ylab="X",main="sigma=2")
plot(1:N,rw6$x,type="l",xlab="sample",ylab="X",main="sigma=4")
plot(1:N,rw7$x,type="l",xlab="sample",ylab="X",main="sigma=8")

a <- c(.05, seq(.1, .9, .1), .95)
Q <- qskewlap(a)
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x, rw5$x, rw6$x, rw7$x)
mc <- rw[501:N, ]
Qrw <- apply(mc, 2, function(x) quantile(x, a))
print(round(cbind(Q, Qrw), 3)) 

a <- ppoints(100)
qqplot(qskewlap(a) ,quantile(rw3$x[501:N], a),xlab="laplace quantile",ylab="sample quantile",main="sigma=0.5")  #q-q plot
lines(c(qskewlap(0.005),qskewlap(0.995)),c(qskewlap(0.005),qskewlap(0.995)))
qqplot(qskewlap(a) ,quantile(rw4$x[501:N], a),xlab="laplace quantile",ylab="sample quantile",main="sigma=1")
lines(c(qskewlap(0.005),qskewlap(0.995)),c(qskewlap(0.005),qskewlap(0.995)))
qqplot(qskewlap(a) ,quantile(rw5$x[501:N], a),xlab="laplace quantile",ylab="sample quantile",main="sigma=2")
lines(c(qskewlap(0.005),qskewlap(0.995)),c(qskewlap(0.005),qskewlap(0.995)))
qqplot(qskewlap(a) ,quantile(rw6$x[501:N], a),xlab="laplace quantile",ylab="sample quantile",main="sigma=4")
lines(c(qskewlap(0.005),qskewlap(0.995)),c(qskewlap(0.005),qskewlap(0.995)))
qqplot(qskewlap(a) ,quantile(rw7$x[501:N], a),xlab="laplace quantile",ylab="sample quantile",main="sigma=8")
lines(c(qskewlap(0.005),qskewlap(0.995)),c(qskewlap(0.005),qskewlap(0.995)))


## ------------------------------------------------------------------------
x<-c(0.01,0.1,0.25,0.5,1,2,4,10,100)
for(i in 1:9){                                  # 比较x和exp(log(x))              
  print(isTRUE(x[i]==exp(log(x[i]))))
  print(isTRUE(all.equal(x[i],exp(log(x[i])))))
}
for(i in 1:9){                                  # 比较x和log(exp(x))                           
  print(isTRUE(x[i]==log(exp(x[i]))))
  print(isTRUE(all.equal(x[i],log(exp(x[i])))))
}
for(i in 1:9){                                  # 比较log(exp(x))和exp(log(x))  
  print(isTRUE(exp(log(x[i]))==log(exp(x[i]))))
  print(isTRUE(all.equal(exp(log(x[i])),log(exp(x[i])))))
}

## ------------------------------------------------------------------------
kk<-c(4:25,100,500,1000)
solution5<-solution4<-numeric(25)
for(i in 1:25){                                 #第4题的方程的解
f<-function(x){
  pt(sqrt(x^2*(kk[i]-1)/(kk[i]-x^2)),df=kk[i]-1)-pt(sqrt(x^2*kk[i]/(kk[i]+1-x^2)),df=kk[i])
}
solution4[i] <- uniroot(f,c(0.0001,min(sqrt(kk[i]-0.001),3)))$root
}

g<-function(x,k){                             #两个积分差的函数（系数做了处理）
  if(x>=sqrt(k)-0.001){
    return(NA)
  }
  if(k<=3){
    return(0)
  }
  c1<-sqrt(k/(k-1))
  c2<-exp(2*lgamma(k/2)-lgamma((k+1)/2)-lgamma((k-1)/2))
  f<-function(u,K){
    (1+u*u/(K-1))^(-K/2)
  }
  in1<-integrate(f,lower=0,upper=sqrt(x*x*(k-1)/(k-x*x)),K=k)$value
  in2<-integrate(f,lower=0,upper=sqrt(x*x*k/(k+1-x*x)),K=k+1)$value
  return(c1*c2*in1-in2)
}
for(i in 1:25){                             #二分法求解                 
  x1<-1
  x2<-1.99
  for(j in 1:24){                           #二分法循环24次，误差<2^-24<10^-7
    x<-(x1+x2)/2
    if(g(x,kk[i])*g(x1,kk[i])<0){
      x2<-x
    }
    else{
      x1<-x
    }
  }
  solution5[i]<-x
}
print(cbind(kk,solution4,solution5))

## ------------------------------------------------------------------------
na<-28
nb<-24
no<-41
nab<-70
library(rootSolve)
library(stats4)
f<-function(p ){                                       #初始最大似然方程
  f1<-na/p[1]-na/(2-p[1]-2*p[2])-2*nb/(2-2*p[1]-p[2])-2*no/(1-p[1]-p[2])+nab/p[1]
  f2<-nb/p[2]-nb/(2-p[2]-2*p[1])-2*na/(2-2*p[2]-p[1])-2*no/(1-p[1]-p[2])+nab/p[2]
  c(F1=f1,F2=f2)
}
q1 <-multiroot(f = f, start = c(0.4,0.4))$root
q2<-numeric(2)
logl<-numeric(20)
n=0
while(sum(abs(q1-q2))>10^(-8)){                       #迭代
  n=n+1
  q2[1]<-q1[1]
  q2[2]<-q1[2]
  naa<-na*q2[1]/(2-q2[1]-2*q2[2])
  nbb<-na*q2[2]/(2-q2[2]-2*q2[1])

  f<-function(p){                                  #最大似然估计方程
    f1<-na/p[1]-na/(1-p[1]-p[2])-nb/(1-p[1]-p[2])+nab/p[1]+naa/p[1]-naa/(1-p[1]-p[2])-2*no/(1-p[1]-p[2])-nbb/(1-p[1]-p[2])
    f2<-nb/p[2]-nb/(1-p[1]-p[2])-na/(1-p[1]-p[2])+nab/p[2]+nbb/p[2]-nbb/(1-p[1]-p[2])-2*no/(1-p[1]-p[2])-naa/(1-p[1]-p[2])
    c(F1=f1,F2=f2)
  }
  q1 <-multiroot(f = f, start = c(q2[1],q2[2]))$root                          
  logl[n]<-2*naa*log(q1[1])+2*nbb*log(q1[2])+2*no*log(1-q1[1]-q1[2])+(na-naa)*log(q1[1]*(1-q1[1]-q1[2]))+(nb-nbb)*log(q1[2]*(1-q1[1]-q1[2]))+nab*log(q1[1]*q1[2])                                #对数最大似然函数值   
}
q1                                              #最终结果p，q
logl[1:n]
plot(1:n,logl[1:n],type="l")

## ------------------------------------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
la1<-lapply(formulas, lm, data = mtcars)
la1
fl1<-list(length(formulas))
for(i in seq_along(formulas)){
  fl1[[i]]<-lm(formulas[[i]],data=mtcars)
}
fl1
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
la2<-lapply(bootstraps,lm,formula=mpg ~ disp)
la2
fl2<-list(length(bootstraps))
for (i in seq_along(bootstraps)){
  fl2[[i]]<- lm(mpg ~ disp,data=bootstraps[[i]])
}
fl2

## ------------------------------------------------------------------------
library(Rcpp)
library(SC19025)
library(microbenchmark)
sigma=c(0.05,0.5,2,16)              #the variances we choose
N<-2000                                  #sample size              
x0<-25
rw.Metropolis <- function( sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= exp(-abs(y)) / exp(-abs(x[i-1]))  ){
x[i] <- y
}
else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
}
cppFunction('NumericVector MetropolisC(double sigma, double x0, int N){
double y;
NumericVector x(N+1), u(N);
int k = 0;
x[0] = x0;
u = runif(N,0,1);
for(int i = 0; i < N; i++){
  y = rnorm(1, x[i], sigma)[0];
  if(u[i+1] <= (exp(-abs(y))/exp(-abs(x[i])))){
    x[i+1] = y;
    k++;
  } else {
    x[i+1] = x[i];
  }
}
  x[N] = k;
  return(x);
}')

ts <- microbenchmark(rw1 <- rw.Metropolis( sigma[1], x0, N),
rw2 <- rw.Metropolis( sigma[2], x0, N),
rw3 <- rw.Metropolis( sigma[3], x0, N),
rw4 <- rw.Metropolis( sigma[4], x0, N),
rwc1<-MetropolisC(sigma[1],x0,N),
rwc2<-MetropolisC(sigma[2],x0,N),
rwc3<-MetropolisC(sigma[3],x0,N),
rwc4<-MetropolisC(sigma[4],x0,N))
summary(ts)[,c(1,3,5,6)]

print(c(1-rwc1[N+1]/N,1- rwc2[N+1]/N,1- rwc3[N+1]/N, 1-rwc4[N+1]/N)) 

plot(1:N,rwc1[1:N],type="l",xlab="sample",ylab="X",main="sigma=0.05")    #the chains generated
plot(1:N,rwc2[1:N],type="l",xlab="sample",ylab="X",main="sigma=0.5")
plot(1:N,rwc3[1:N],type="l",xlab="sample",ylab="X",main="sigma=2")
plot(1:N,rwc4[1:N],type="l",xlab="sample",ylab="X",main="sigma=16")

a <- ppoints(100)
qqplot(quantile(rwc1[501:N], a) ,quantile(rw1$x[501:N], a),xlab="C program",ylab="R program",main="sigma=0.05")
lines(quantile(rwc1[501:N], a),quantile(rwc1[501:N], a))
qqplot(quantile(rwc2[501:N], a) ,quantile(rw2$x[501:N], a),xlab="C program",ylab="R program",main="sigma=0.5")
lines(quantile(rwc2[501:N], a),quantile(rwc2[501:N], a))
qqplot(quantile(rwc3[501:N], a) ,quantile(rw3$x[501:N], a),xlab="C program",ylab="R program",main="sigma=2")
lines(quantile(rwc3[501:N], a),quantile(rwc3[501:N], a))
qqplot(quantile(rwc4[501:N], a) ,quantile(rw4$x[501:N], a),xlab="C program",ylab="R program",main="sigma=16")
lines(quantile(rwc4[501:N], a),quantile(rwc4[501:N], a))

