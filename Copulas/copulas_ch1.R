install.packages('copula')
library(copula)

d=2
ic=indepCopula(dim=d)

set.seed(2008)
u=runif(d)
(Pi=pCopula(u,copula=ic)) # the value of the independence copula at u

stopifnot(all.equal(Pi, prod(u))) #check nemerical equality of the samples

wireframe2(ic,FUN=pCopula,  #surface plot of the independence copula
           col.4=adjustcolor('black', alpha.f=0.25))

contourplot2(ic, FUN=pCopula) #contour plot of the independence copula

a=c(1/4,1/2) #lower left end point
b=c(1/3,1) #upper right end point
stopifnot(0<=a,a<=1,0<=b,b<=1,a<=b) #check
p=(b[1]-a[1])*(b[2]-a[2]) #maual computation
stopifnot(all.equal(prob(ic,l=a,u=b),p)) #check

n=1000 #sample size
set.seed(271) # set a seed
U=rCopula(n, copula=ic) # generate a sample of the independence copula
plot(U, xlab=quote(U[1]), ylab=quote(U[2]))

set.seed(271)
stopifnot(all.equal(U,matrix(runif(n*d), nrow=n)))

set.seed(314)
U=rCopula(1e6,copula=ic) #large sample size for good approximation
## Approximate the Pi-volume by the aforementioned proportion
p.sim=mean(a[1]<U[,1] & U[,1]<=b[1] & a[2]<U[,2] & U[,2]<=b[2])
stopifnot(all.equal(p.sim ,p, tol=1e-2))


d=2  #dimension
theta=-9 #copula parameter
fc=frankCopula(theta,dim=d) #define a Frank copula

set.seed(2010)
n=5 #number of evaluation points
u=matrix(runif(n*d), nrow=n)  # n randompoints in [0,1]^d
pCopula(u,copula=fc) #copula values at u

dCopula(u, copula=fc) #density values at u

wireframe2(fc, FUN=pCopula, #wireframe plot(copula)
           draw.4.pCoplines=FALSE)
wireframe2(fc,FUN=dCopula, delta=0.001) #wireframe plot(density)
contourplot2(fc, FUN=pCopula)  #contour plot(copula)
contourplot2(fc,FUN=dCopula, n.grid=72, #contour plot(density)
             lwd=1/2)


library(lattice)
set.seed(1946)
n=1000
U=rCopula(n,copula=fc)
U0=rCopula(n, copula=setTheta(fc, value=0))
U9=rCopula(n, copula=setTheta(fc, value=9))
plot(U, xlab=quote(U[1]), ylab=quote(U[2]))
plot(U0, xlab=quote(U[1]), ylab=quote(U[2]))
plot(U9, xlab=quote(U[1]), ylab=quote(U[2]))


d=3
cc=claytonCopula(4,dim=d) #theta=4

set.seed(2013)
n=5
u=matrix(runif(n*d),nrow=n) #random points in the unit hypercube
pCopula(u,copula=cc) #copula values at u


dCopula(u,copula=cc)


set.seed(271)
U=rCopula(1000, copula=cc)
splom2(U,cex=0.3, col.mat='black')


gc=gumbelCopula(3) #theta=3(note the defalt dim=2)

set.seed(1993)
U=rCopula(1000, copula=gc)
plot(U,xlab=quote(U[1]), ylab=quote(U[2]))
wireframe2(gc, dCopula, delta=0.025)


# the frechet-hoeffdings bounds

set.seed(1980)
U=runif(100)
plot(cbind(U,1-U), xlab=quote(U[1]), ylab=quote(U[2]))
plot(cbind(U,U), xlab=quote(U[1]), ylab=quote(U[2]))


u=seq(0,1,length.out=40)# subdivision points in each dimension
u12=expand.grid('u[1]'=u,'u[2]'=u) #build a grid
W=pmax(u12[,1]+u12[,2]-1,0) #values of W on grid
M=pmin(u12[,1], u12[,2]) # values of M on grid
val.W=cbind(u12, 'W(u[1], u[2])'=W) #append grid
val.M=cbind(u12, 'M(u[1] ,u[2])'=M)  #append grid

wireframe2(val.W)
wireframe2(val.M)
contourplot2(val.W, xlim=0:1, ylim=0:1)
contourplot2(val.M, xlim=0:1, ylim=0:1)

# marshall-olkin copula
C=function(u, alpha){
  pmin(u[,1]*u[,2]^(1-alpha[2]), u[,1]^(1-alpha[1])*u[,2])}
alpha=c(0.2,0.8)
val=cbind(u12,'C(u[1],u[2])'=C(u12,alpha=alpha)) #append C values

#generate data

set.seed(712)
V=matrix(runif(1000*3), ncol=3)
U=cbind(pmax(V[,1]^(1/(1-alpha[1])), V[,3]^(1/alpha[1])),
        pmax(V[,2]^(1/(1-alpha[2])), V[,3]^(1/alpha[2])))

#plots
wireframe2(val)
plot(U,xlab=quote(U[1]), ylab=quote(U[2]))

install.packages('mvtnorm')
library(mvtnorm)

d=2 #dimension
rho= 0.7  #off-diagonal entry of the correlation matrix P
P=matrix(rho, nrow=d, ncol=d) #build the correlation matrix P
diag(P)=1
set.seed(64)
u=runif(d) 
x=qnorm(u)
pmvnorm(upper=x, corr=P) #evaluate the copula C at u


nc=normalCopula(rho)
pCopula(u, copula=nc)

nu=3  #degrees of freedom
x.=qt(u,df=nu)
pmvt(upper=x., corr=P, df=nu)  #evaluate the t copula at u

try(pmvt(upper=x., corr=P, df=3.5))


tc=tCopula(rho, dim=d, df=nu)
pCopula(u, copula=tc)


H.obj=mvdc(claytonCopula(1), margins=c('norm', 'exp'),
           paramMargins=list(list(mean=1, sd=2), list(rate=3)))

set.seed(1979)
z=cbind(rnorm(5,mean=1, sd=2), rexp(5,rate=3)) # evaluation points
pMvdc(z,mvdc=H.obj)

dMvdc(z,mvdc=H.obj)  #values of the corresponding density at z

set.seed(1975)
X=rMvdc(1000, mvdc=H.obj)


plot(X, cex=0.5, xlab=quote(X[1]), ylab=quote(X[2]))
contourplot2(H.obj, FUN=dMvdc, xlim=range(X[,1]), ylim=range(X[,2]),
             n.grid=257)

# install.packages('nor1mix')
library(nor1mix)

#normal mixtures
nm1=norMix(c(1,-1), sigma=c(.5,1), w=c(.2,.8))
plot(nm1, p.comp=TRUE  )
nm2=norMix(c(0,2), sigma=c(1.5,1), w=c(.3,.7))
plot(nm2, p.comp=TRUE)

H.obj.m=mvdc(claytonCopula(1), margins=c('norMix', 'norMix'),
             paramMargins=list(nm1,nm2))

set.seed(271)

X=rMvdc(1000,mvdc=H.obj.m)

plot(X, cex=0.5, xlab=quote(X[1]), ylab=quote(X[2]))
contourplot2(H.obj.m, FUN=dMvdc, xlim=range(X[,1]), ylim=range(X[,2]),
             n.grid=129)


## Define parameters of the three margins

th=2.5 #pareto parameter
m=10 #mean of the lognormal
v=20 #variance of the lognormal
s=4 #shape of the gamma underlying the loggamma
r=5 #rate of the gamma underying the loggamma

## Define list of marginal dfs

qF=list(qPar=function(p) (1-p)^(-1/th)-1,
        qLN=function(p) qlnorm(p, meanlog=log(m)-log(1+v/m^2)/2,
                               sdlog=sqrt(log(1+v/m^2))),
        qLg=function(p) exp(qgamma(p,shape=s, rate=r)))

# generate the data

set.seed(271)
X=sapply(qF, function(mqf) mqf(runif(2500)))


##title nonparametric VaR estimate uder a t copula
## param X loss matrix
## param alpha confidence level(s)
## param rho correlation parameter of the t copula
## param df degrees of freedom parameter of the t copula
## return nonparametric Var estimate under the t copula(nemeric)

VaR=function(X, alpha, rho, df=3.5){
  stopifnot(is.matrix(X), 0<=rho, rho<=1, length(rho)==1,
            0<alpha, alpha<1, length(alpha)>=1)
  n=nrow(X) #sample size
  d=ncol(X) #dimension
  ## simulate from a t copula with d.o.f. parameter 3.5 and exchangeable
  # correlation matrix with off-diagonal entry rho. also compute the
  #componentwise ranks.
  #Note: we can set the seed here as we can estimate VaR for all
  #confidence levels based on the same  copula sample.
  # we even should set seed here to minimize the variance 
  # of the estimator and make the results more comparable
  
  set.seed(271)
  U=rCopula(n, copula=tCopula(rho, dim=d, df=df))
  rk=apply(U,2,rank)
  #componentwise reorder the data according to these ranks to
  #mimic the corresponding t copula dependence among the losses
  Y=sapply(1:d, function(j) sort(X[,j])[rk[,j]])
  #build row sums to mimic a sample from the distribution of the
  #sum under the corresponding t copula
  S=rowSums(Y)
  #nonparametrically estimate Var for all confidence levels alpha
  # quantile function here
  quantile(S, probs=alpha, type=1, names=FALSE)
  
}

alpha=c(
  0.001, 0.01, 0.05,0.1, 0.2, 0.3,0.4,0.5,
  0.6, 0.7, 0.8, 0.9,0.95,0.99,0.999
)  # confidence levels
rho=seq(0,1,by=0.1)  #parameter of the homogeneous t copula
grid=expand.grid('alpha'=alpha, 'rho'=rho)[,2:1] # build a grid
VaR.fit=sapply(rho, function(r)
  VaR(X, alpha=alpha, rho=r)) #(alpha, rho)
res=cbind(grid, "VaR[alpha](L^'+')"=as.vector(VaR.fit))

wireframe2(res)
# install.packages('qrmtools')
library(TTR)
library(quantmod)
library(qrmtools)

worst.VaR=sapply(alpha, function(a) mean(ARA(a,qF=qF)$bounds))
plot(alpha, worst.VaR, type='b', col=2,
     xlab=quote(alpha), ylab=quote(VaR[alpha](L^'+')),
     ylim=range(VaR.fit, worst.VaR)) #computed with the ARA
lines(alpha, apply(VaR.fit, 1, max), type='b', col=1) #sumulated
legend('topleft', bty='n', lty=rep(1,2), col=2:1,
       legend=c(expression('Worst'~VaR[alpha]~'according to ARA()'),
                expression('Worst'~VaR[alpha]~'under'~t[3.5]~'copulas')))


## computing worst VaR in the three-dimensional case
wVaR=ARA(0.99, qF=qF) #compute worst Var(bounds)
X=wVaR[['X.rearranged']]$up #extract rearranged matrix(upper bound)
U=pobs(X) #compute psudo-observations
pairs2(U) # approx. sample of a copula leading to worst VaR for our marg. dfs

##computing worst VaR in the bivariate case
wVar.=ARA(0.99, qF=qF[1:2]) #compute worst Var(bounds)
X.=wVar.[['X.rearranged']]$up #extract rearranged matrix(upper bound)
U.=pobs(X.) #compute psudo-observations
plot(U., xlab=quote(U[1]), ylab=quote(U[2]))



n=1000
d=2 #dimension
rho=0.7 # off-diagonal entry in the correlation matrix P
P=matrix(rho, nrow=d, ncol=d) #build the correlation matrix P
diag(P)=1 
nu=3.5 #degrees of freedom
set.seed(271) 
X=rmvt(n, sigma=P, df=nu) # n ind. multivariate t observations
U=pt(X,df=nu) #n ind. realizations from the corresponding copula


set.seed(271)
U.=rCopula(n, tCopula(rho, dim=d, df=nu))
stopifnot(all.equal(U,U.)) #test of equality

plot(U., xlab=quote(U[1]), ylab=quote(U[2]))
plot(U, xlab=quote(U[1]), ylab=quote(U[2]))


##plot function highlighting three point

plotABC=function(x, ind3, col=adjustcolor('black', 1/2), pch=19,...){
  cols=adjustcolor(c('red', 'blue', 'magenta'), offset=-c(1,1,1,1.5)/4)
  par(pty='s')
  plot(x,col=col, asp=1,...)
  xy=x[ind3, , drop=FALSE]
  points(xy,pch=pch, col=cols)
  text(xy,label=names(ind3), adj=c(0.5, -0.6), col=cols, font=2)
}

ind3=c(A=725, B=351, C=734) #found via plot(X):identify(X)
##Scatter plot of observations from the multivariate t distribution
plotABC(X, ind3=ind3, xlab=quote(X[1]), ylab=quote(X[2]))
##scatter plot of observations from the corresponding t copula
plotABC(U, ind3=ind3, xlab=quote(U[1]), ylab=quote(U[2]))
## scatterplot of observation from the meta-t distribution

Y=qnorm(U)  # transform U(t copula) to normal margins
plotABC(Y, ind3=ind3, xlab=quote(Y[1]), ylab=quote(Y[2]))



rho=0.6
P=matrix(c(1,rho, rho,1),ncol=2) #the correlation matrix
C=function(u) pCopula(u, copula=normalCopula(rho)) #normal copula
Htilde=function(x){
  apply(cbind(log(x[,1]), -log((1-x[,2])/x[,2])),1,function(x.) 
    pmvnorm(upper=x., corr=P))}
qF1tilde=function(u) exp(qnorm(u))
qF2tilde=function(u) 1/(1+exp(-qnorm(u)))
Ctilde=function(u) Htilde(cbind(qF1tilde(u[,1]), qF2tilde(u[,2])))
set.seed(31)
u=matrix(runif(5*2), ncol=2) # 5 random evaluation points
stopifnot(all.equal(Ctilde(u), C(u)))

set.seed(721)
X=rmvnorm(1000, mean=c(0,0), sigma=P) #sample from N(0,P)
## 'sample' the copula of X directly
U=pnorm(X)
## transform the sample X componentwise
TX=cbind(exp(X[,1]), plogis(X[,2])) #note:plogis(x)=1/(1+exp(-x))
##apply the marginal dfs to get a sample from the copula of TX
##note: qlogis(p)== logit(p)== log(p/(1-p))
V=cbind(pnorm(log(TX[,1])), pnorm(qlogis(TX[,2])))
stopifnot(all.equal(V,U))  #the samples of the two copulas are the same


cc=claytonCopula(2)
set.seed(271)
U=rCopula(1000, copula=cc) #sample from the clayton copula
V=1-U  #sample from the survival clayton copulas
plot(U, xlab=quote(U[1]), ylab=quote(U[2])) # scatter plot
plot(V, xlab=quote(V[1]), ylab=quote(V[2])) #for the survival copula

wireframe2(cc , FUN=dCopula, delta=0.025)
wireframe2(rotCopula(cc), FUN=dCopula, delta=0.025)

contourplot2(tCopula(0.7, df=3.5), FUN=dCopula, n.grid=64, lwd=1/2)
contourplot2(gumbelCopula(2), FUN=dCopula, n.grid=64, lwd=1/4,
             pretty=FALSE, cuts=42,
             col.regions=gray(seq(0.5,1, length.out=128)))


##Evaluate the density of C for h_1(u)=2*u*(u-1/2)*(u-1),
##h_2(u)=theta*u*(1-u) and two different thetas
u=seq(0,1,length.out=20) # subdivision points in each dimension
u12=expand.grid('u[1]'=u, 'u[2]'=u) #build a grid
dC=function(u,th) 1+th*(6*u[,1]*(u[,1]-1)+1)*(1-2*u[,2])
wireframe2(cbind(u12, 'c(u[1],u[2])'=dC(u12,th=-1)))
wireframe2(cbind(u12, 'c(u[1],u[2])'=dC(u12,th=1)))



n=1000
set.seed(314)
Z=rnorm(n)
U=runif(n)
V=rep(1,n)
V=rep(1,n)
V[U<1/2]=-1  #V in{-1,1}, each with probability 1/2
X=cbind(Z,Z*V) #(X_1,X_2)
stopifnot(cor.test(X[,1], X[,2])$p.value>=0.05) #H0: cor=0 not rejected
Y=matrix(rnorm(n*2), ncol=2) #independent N(0,1)
## plots

plot(X, xlab=quote(X[1]), ylab=quote(X[2]))
plot(Y, xlab=quote(Y[1]), ylab=quote(Y[2]))

corBoundLN=function(s, bound=c('max', 'min')){
  ## s=(sigma_1, sigma_2)
  if(!is.matrix(s)) s=rbind(s)
  bound=match.arg(bound)
  if(bound=='min') s[,2]=-s[,2]
  (exp((s[,1]+s[,2])^2/2)-exp((s[,1]^2+s[,2]^2)/2))/
    sqrt(expm1(s[,1]^2)*exp(s[,1]^2)*expm1(s[,2]^2)*exp(s[,2]^2))}

##evaluate correlation bounds on a grid
s=seq(0.01, 5, length.out=20) #subdivision points in each dimension
s12=expand.grid('sigma[1]'=s, 'sigma[2]'=s)  #build a grid

##plots

wireframe2(cbind(s12, 'underline(Cor)(sigma[1],sigma[2])'=
                   corBoundLN(s12, bound='min')))
wireframe2(cbind(s12, 'bar(Cor)(sigma[1],sigma[2])'=corBoundLN(s12)))


theta=-0.7
stopifnot(all.equal(rho(normalCopula(theta)), 6/pi*asin(theta/2)))
stopifnot(all.equal(tau(normalCopula(theta)), 2/pi*asin(theta)))
theta=2
stopifnot(all.equal(tau(claytonCopula(theta)), theta/(theta+2)))
stopifnot(all.equal(tau(gumbelCopula(theta)),1-1/theta))

theta=(0:8)/16
stopifnot(all.equal(iRho(normalCopula(), rho=6/pi*asin(theta/2)), theta))
stopifnot(all.equal(iTau(normalCopula(), tau=2/pi*asin(theta)), theta))
theta=1:20
stopifnot(all.equal(iTau(claytonCopula(), theta/(theta+2)), theta))
stopifnot(all.equal(iTau(gumbelCopula(), 1-1/theta), theta))


theta=3
iRho(claytonCopula(), rho=rho(claytonCopula(theta)))


#sparman's rho can be estimated with cor(, method='spearman)
theta=iRho(claytonCopula(), rho=0.6)  #true spearman's rho=0.6
set.seed(974)
U=rCopula(1000, copula=claytonCopula(theta)) 
rho.def=cor(apply(U,2,rank))[1,2] #sparman's rho manually
rho.R=cor(U,method = 'spearman')[1,2] #spaerman's rho from R
stopifnot(all.equal(rho.def, rho.R)) #the same
rho.R #close to -0.6


#kendall's tau can be estimated with cor(, method='kendall'):
theta=iTau(normalCopula(), tau=-0.5) #true kendall's tau=-0.5
set.seed(974)
U=rCopula(1000, copula=normalCopula(theta))
p.n=0
for(i in 1:(n-1)) #number of concordant pairs(obviously inefficient)
  for(j in (i+1):n)
    if(prod(apply(U[c(i,j),], 2, diff)) >0) p.n=p.n+1
tau.def=4*p.n/(n*(n-1))-1 #kendall's tau manually
tau.R=cor(U, method='kendall')[1,2] #kendall's tau from R
stopifnot(all.equal(tau.def, tau.R)) #the same
tau.R #close to -0.5

set.seed(75)
X=rnorm(100)
Y=-X^3 #perfect negative dependence
rho.counter=cor(X,Y, method='spearman')
tau.counter=cor(X,Y, method='kendal')
stopifnot(rho.counter==-1, tau.counter==-1)
Z=exp(X)
rho.co=cor(X,Z, method='spearman')
tau.co=cor(X,Z,method='kendall')
stopifnot(rho.co==1, tau.co==1)


rho=seq(-1,1,by=0.01) #correlation parameters of normal copulas
rho.s=(6/pi)*asin(rho/2) #corresponding Spearman's rho
tau=(2/pi)*asin(rho) #corresponding kendall's tau
plot(rho, rho.s, type='l', col=2, lwd=2,
     xlab=expression('Correlation parameter'~rho~'of'~C[rho]^n),
     ylab=expression('Corresponding'~rho[s]~'and'~tau))
abline(a=0, b=1, col=1, lty=2, lwd=2)
lines(rho, tau, col=3, lwd=2)
legend('bottomright', bty='n', col=1:3, lty=c(2,1,1), lwd=2,
       legend=c('Diagonal', expression(rho[s]), expression(tau)))
plot(rho, rho.s-rho, type='l', yaxt='n',lwd=2,
     xlab=expression(rho), ylab=expression(rho[s]-rho))
mdiff=max(rho.s-rho)
abline(h=c(-1,1)*mdiff,lty=2, lwd=2)
rmdiff=round(mdiff,4)
axis(2,at=c(-mdiff, -0.01, 0, 0.01, mdiff),
     labels=as.character(c(-rmdiff, -0.01, 0, 0.01, rmdiff)))

## kendall's tau and corresponding copula parameters
tau=0.7
th.n=iTau(normalCopula(),tau=tau)
th.t=iTau(tCopula(df=3), tau=tau)
th.c=iTau(claytonCopula(), tau=tau)
th.g=iTau(gumbelCopula(), tau=tau)

#samples from the corresponding 'mvdc' objects
set.seed(271)
n=10000
N01m=list(list(mean=0, sd=1), list(mean=0, sd=1))
X.n=rMvdc(n, mvdc=mvdc(normalCopula(th.n), c('norm','norm'),N01m))
X.t=rMvdc(n,mvdc=mvdc(tCopula(th.t,df=3), c('norm','norm'), N01m))
X.c=rMvdc(n,mvdc=mvdc(claytonCopula(th.c), c('norm', 'norm'), N01m))
X.g=rMvdc(n, mvdc=mvdc(gumbelCopula(th.g), c('norm', 'norm'), N01m))


plotCorners=function(X,qu, lim, smooth=False, ...){
  plot(X, xlim=lim, ylim=lim, xlab=quote(X[1]), ylab=quote(X[2]),
       col=adjustcolor('black', 0.5), ...)#or pch=16
  abline(h=qu, v=qu, lty=2, col=adjustcolor('black', 0.6))
  ll=sum(apply(X<=qu[1], 1, all))*100/n
  ur=sum(apply(X>=qu[2],1,all))*100/n
  mtext(sprintf('Lower left:%.2f%%, upper right:%.2f%%', ll,ur),
        cex=0.9, side=1, line=-1.5)
  invisible()}

a.=0.005
q=qnorm(c(a.,1-a.))
lim=range(q,X.n, X.t, X.g)
lim=c(floor(lim[1]), ceiling(lim[2]))
plotCorners(X.n,qu=q,lim=lim, cex=0.4)
plotCorners(X.t,qu=q,lim=lim, cex=0.4)
plotCorners(X.c,qu=q,lim=lim, cex=0.4)
plotCorners(X.g,qu=q,lim=lim, cex=0.4)

theta=3
lam.c=lambda(claytonCopula(theta))
stopifnot(all.equal(lam.c[['lower']],2^(-1/theta)),
          all.equal(lam.c[['upper']],0))

lam.g=lambda(gumbelCopula(theta))
stopifnot(all.equal(lam.g[['lower']],0),
          all.equal(lam.g[['upper']], 2-2^(1/theta)))

rho=0.7
nu=3

lam.n=lambda(normalCopula(rho))
stopifnot(all.equal(lam.n[['lower']],0),
          all.equal(lam.n[['lower']],lam.n[['upper']]))

lam.t=lambda(tCopula(rho,df=nu))
stopifnot(
  all.equal(lam.t[['lower']],
            2*pt( -sqrt((nu+1)*(1-rho)/(1+rho)) , df=nu+1 )),
  all.equal(lam.t[['lower']], lam.t[['upper']]))


#coefficient of tail dependence as a function of rho
rho=seq(-1,1,by=0.01)
nu=c(3,4,8,Inf)
n.nu=length(nu)
lam.rho=sapply(nu,function(nu.) # (rho,nu)matrix
  sapply(rho, function(rho.) lambda(tCopula(rho., df=nu.))[['lower']]))
expr.rho=as.expression(lapply(1:n.nu,function(j)
  bquote(nu==.(if(nu[j]==Inf) quote(infinity) else nu[j]))))
matplot(rho, lam.rho,type='l', lty=1, lwd=2,col=1:n.nu,
        xlab=quote(rho), ylab=quote(lambda))
legend('topleft', legend=expr.rho, bty='n', lwd=2, col=1:n.nu)
##coefficient of tail dependence as a function of nu
nu.=c(seq(3,12,by=0.2), Inf)
rho.=c(-1,-0.5, 0, 0.5, 1)
n.rho=length(rho.)
lam.nu=sapply(rho., function(rh) #(nu,rho)matrix 
  sapply(nu., function(nu) lambda(tCopula(rh,df=nu))[['lower']]))
expr=as.expression(lapply(1:n.rho, function(j) bquote(rho==.(rho.[j]))))
matplot(nu., lam.nu, type='l', lty=1, lwd=2, col=1:n.rho,
        xlab=quote(nu), ylab=quote(lambda))
legend('right', expr, bty='n', lwd=2, col=1:n.rho)



##note: all calculations here are deterministic
u=seq(0.95, to=0.9999, length.out=128) # levels u of P(U_1>u, U_2>u)
rho=c(0.75,0.5) #correlation parameter rho
nu=c(3,4,8,Inf) #degrees of freedom
len=length(rho)*length(nu)
tail.prob=matrix(u,nrow=length(u), ncol=1+len)# tail probabilities
expr=vector('expression', length=len) #vector of expressions
ltys=cols=numeric(len)

for(i in seq_along(rho)){ #rho
  for(j in seq_along(nu)){ #degrees of freedom
    k=length(nu)*(i-1)+j
    ##create the copula
    cop=ellipCopula('t', param=rho[i], df=nu[j])
    ##Evaluate P(u_1>u,U_2>u)=P(U_1<=1-u, U_2<=1-u)
    tail.prob[,k+1]=pCopula(cbind(1-u, 1-u), copula=cop)
    ## Create plot informaion
    expr[k]=as.expression(
      substitute(group('(',list(rho,nu),')')==
                   group('(',list(RR,NN),')'),
                 list(RR=rho[i],
                      NN=if(is.infinite(nu[j]))
                        quote(infinity) else nu[j])))
    ltys[k]=length(rho)-i+1
    cols[k]=j}}

##standardize w,r,t gauss case
tail.prob.fact=tail.prob #for comparison to Gauss case
tail.prob.fact[,2:5]=tail.prob[,2:5]/tail.prob[,5]
tail.prob.fact[,6:9]=tail.prob[,6:9]/tail.prob[,9]

#plot tail probabilities
matplot(tail.prob[,1], tail.prob[,-1], type='l', lwd=2, lty=ltys,
        col=cols, xlab=quote(P(U[1]>u,U[2]>u)~~"as a function of u"),
        ylab='')
legend('topright', expr, bty='n', lwd=2, lty=ltys, col=cols)
##plot standardized tail probabilities
matplot(tail.prob.fact[,1], tail.prob.fact[,-1], log='y', type='l',
        lty=ltys, col=cols, lwd=(wd=2*c(1,1,1,1.6,1,1,1,1)),
        xlab=quote(P(U[1]>u,U[2]>u)~~
                     'as a function of u standardized by Gauss case'),
        ylab='')
legend('topleft', expr, bty='n', lwd=wd, lty=ltys, col=cols)