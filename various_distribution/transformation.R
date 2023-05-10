# 모든 데이터가 classical distribution에 적합하진 않는다.
# 그래서 필요한것이 transformation이다.
# transformation의 방법

# X+c
# scaling lambda*X
# Power X^alpha(멱변환)
# Inverse (->1/X)
# 로그화
# 지수화 
# X/(1-X)로 odds비 산출 

# 감마분포 transform 예시

f<- function(x) dgamma(x,2)
f1 <- function(x) f(x-1)
f2 <- function(x) f(x/2)/2
f3 <- function(x) 2*x*f(x^2)
f4 <- function(x) f(1/x)/x^2
f5 <- function(x) f(exp(x))*exp(x)
f6 <- function(x) f(log(x))/x
x=seq(0,10,by=.025)
plot(x,f(x), ylim=c(0, 1.3), xlim=c(0, 10), main="Theoretial densities",
        lwd=2,
       type="l", xlab="x", ylab="")
lines(x,f1(x), lty=2, lwd=2)
lines(x,f2(x), lty=3, lwd=2)
lines(x,f3(x), lty=4, lwd=2)
lines(x,f4(x), lty=1, col="grey", lwd=2)
lines(x,f5(x), lty=2, col="grey", lwd=2)
lines(x,f6(x), lty=3, col="grey", lwd=2)
legend("topright", lty=1:4, col=c(rep("black", 4), rep("grey", 3)),
         leg=c("X","X+1","2X", "sqrt(X)", "1/X", "log(X)", "exp(X)"))
       

#kernel-based densities

set.seed(123)
x <- rgamma(100, 2)
x1 <- x+1
x2 <- 2*x
x3 <- sqrt(x)
x4 <- 1/x
x5 <- log(x)
x6 <- exp(x)
plot(density(x), ylim=c(0, 1), xlim=c(0, 10), main="Empirical densities",
lwd=2, xlab="x", ylab="f_X(x)")
lines(density(x1), lty=2, lwd=2)
lines(density(x2), lty=3, lwd=2)
lines(density(x3), lty=4, lwd=2)
lines(density(x4), lty=1, col="grey", lwd=2)
lines(density(x5), lty=2, col="grey", lwd=2)
lines(density(x6), lty=3, col="grey", lwd=2)
     