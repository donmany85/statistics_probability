#mixed-typed distributions 
#mixed-typed distribution은 연속적이지도 이산적이지도 않은 변수들의 분포라고 한다


# zero-modified beta distribution 

dbeta0M <- function(x, prob, a, b){
  dbeta(x,a,b)*(1-prob)*(x !=1)+prob*(x==1)
}
pbeta0M <- function(q, prob, a, b){
  pbeta(q,a,b)*(1-prob)+prob*(q>=1)
}



#The Maxwell-Boltzmann Bore-Einstein Fermi-Dirac(MBBEFD) distribution
#맥스웰-볼츠만 보어-아인슈타인 페르미-디랙분포
#물리통계학에 쓰이기만 할것같지만 재보험 손실계산에도 쓰인다고 한다(Bernegger 1997)

dMBBEFD <- function(x,a,b){
  -a*(a+1)*b^x*log(b)/(a+b^x)^2+(a+1)*b/(a+b)*(x==1)
}

pMBBEFD <- function(x,a,b){
  a*((a+1)/(a+b^x)-1)*(x<1)+1*(x>=1)
}


#mixed-gamma-pareto distribution

library(actuar)
dmixgampar <- function(x, prob, nu, lambda, alpha, theta){
   prob*dgamma(x, nu, lambda) + (1-prob)*dpareto(x, alpha, theta)}
pmixgampar <- function(q, prob, nu, lambda, alpha, theta){
   prob*pgamma(q, nu, lambda) + (1-prob)*ppareto(q, alpha, theta)}