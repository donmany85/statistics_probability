# 피어슨 체계
# 어떤 pdf fX를 가정할경우

# (1/fX(x))*(dfX(x)/dx)=(a+x)/(c0+c1x+c2x^2)
#
# 0: c1=c2=0 정규분포
# 1: c0+c1x+c2x^2  , fX(x)=K((x-a1)^m1)*(a2-x)^m2 ,부호가 다른두 실근 베타분포
# 2: 1에서 m1=m2=m인 경우
# 3: c2=0 1차 다항함수 c0+c1x 로서 감마분포 fX(x)=K((c0+c1x)^m)*e^(x+c1)
# 4: 다항함수 p(x)=c0+c1x+c2x^2 로 허근,  p(x)=C0+c2(x+C1)^2로
#fX(x)=K(C0+c2(x+C1)^2)e^(k*1/tan((x+c1)/sqrt(c0/c2))) -> Bardoff-Nielsen의 역가우시안분포
# 5: p는 완전제곱으로 p(x)=(x+C1)^2 fX(x)=K((x+C1)^(-1/c2))*e^(k/(x+C1))  
# 5의 특수한 경우로 k=0, c2>0인경우 8,  c2<0 인경우 9
# 6: p는 두 실근 a1, a2를 갖고있는 경우 fX(x)=K((x-a1)^m1)*(x-a2)^m2로 일반화 베타분포
# 7: a=c1=0  fX(x)=K(c0+c2x^2)^(-1/2c2) 인 경우

library(PearsonDS)
x <- seq(-1, 6, by=1e-3)
y0 <- dpearson0(x, 2, 1/2)
y1 <- dpearsonI(x, 1.5, 2, 0, 2)
y2 <- dpearsonII(x, 2, 0, 1)
y3 <- dpearsonIII(x, 3, 0, 1/2)
y4 <- dpearsonIV(x, 2.5, 1/3, 1, 2/3)
y5 <- dpearsonV(x, 2.5, -1, 1)
y6 <- dpearsonVI(x, 1/2, 2/3, 2, 1)
y7 <- dpearsonVII(x, 3, 4, 1/2)
plot(x, y0, type="l", ylim=range(y0, y1, y2, y3, y4, y5, y7),
       ylab="f(x)", main="The Pearson distribution system")
lines(x[y1 != 0], y1[y1 != 0], lty=2)
lines(x[y2 != 0], y2[y2 != 0], lty=3)
lines(x[y3 != 0], y3[y3 != 0], lty=4)
lines(x, y4, col="grey")
lines(x, y5, col="grey", lty=2)
lines(x[y6 != 0], y6[y6 != 0], col="grey", lty=3)
lines(x[y7 != 0], y7[y7 != 0], col="grey", lty=4)
legend("topright", leg=paste("Pearson", 0:7), lty=1:4,
         col=c(rep("black", 4), rep("grey", 4)))
