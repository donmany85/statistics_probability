gen.ruin = function(n, x.cnt, y.cnt, x.p){
  x.cnt.c = x.cnt
  y.cnt.c = y.cnt
  x.rnd = rbinom(n, 1, p=x.p)
  x.rnd[x.rnd==0] = -1
  y.rnd = x.rnd*-1
  x.cum.sum = cumsum(x.rnd)+x.cnt
  y.cum.sum = cumsum(y.rnd)+y.cnt
  
  ruin.data = cumsum(x.rnd)+x.cnt
  
  if( any( which(ruin.data>=x.cnt+y.cnt) ) | any( which(ruin.data<=0) ) ){ cut.data = 1+min( which(ruin.data>=x.cnt+y.cnt), which(ruin.data<=0) )
  
  ruin.data[cut.data:length(ruin.data)] = 0
  
  }
  
  return(ruin.data)
  
}
n.reps = 10000
ruin.sim = replicate(n.reps, gen.ruin(n=1000, x.cnt=5, y.cnt=10, x.p=.6))
ruin.sim[ruin.sim==0] = NA
hist( apply(ruin.sim==15 | is.na(ruin.sim), 2, which.max) , nclass=100, col='8', main="Distribution of Number of Turns",
      xlab="Turn Number")
abline(v=mean(apply(ruin.sim==15 | is.na(ruin.sim), 2, which.max)), lwd=3, col='red')
abline(v=median(apply(ruin.sim==15 | is.na(ruin.sim), 2, which.max)), lwd=3, col='green')
x.annihilation = apply(ruin.sim==15, 2, which.max)
( prob.x.annilate = length(x.annihilation[x.annihilation!=1]) / n.reps )
state.cnt = ruin.sim
state.cnt[state.cnt!=15] = 0
state.cnt[state.cnt==15] = 1
mean.state = apply(ruin.sim, 1, mean, na.rm=T)
plot(mean.state, xlim=c(0,which.max(mean.state)), ylim=c(0,20), ylab="Points", xlab="Number of Plays", pch=16, cex=.5, col='green')
lines(mean.state, col='green')
points(15-mean.state, pch=16, cex=.5, col='blue')
lines(15-mean.state, col='blue')