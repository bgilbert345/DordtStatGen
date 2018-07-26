
library(Matrix)
source("continuous_simulator_Liu1")


#define Qstat
top=function(x){tail(sort(x),causal)}

Qstat = function(gen, phen.perm){
  mu.hat = rep(mean(phen), length(phen))
  S = crossprod(gen, phen.perm-mu.hat)
  t=apply(abs(S),2,top)
  apply(t^2,2,sum)
}




non.settings=c(8,48,100,500)
pow.trials = 10000
permutes = 1000
causal = 8 #number of causal variants
n = 2000 #number of individuals
beta = c(rep(0.75, causal/2),rep(-0.75, causal/2)) #effect size
maf1 = .01
maf2=.001
s = 1 #model standard deviation
mu = 0 #model mean

pow.mat = matrix(nrow=length(non.settings),ncol=pow.trials)

for(k in 1:length(non.settings))
  
{
  
  noncausal= non.settings[k] #number of noncausal variants
  print(noncausal)
  ps = vector("numeric",length=pow.trials)
  
  for(j in 1:pow.trials){
    #simulate data
    sim = cont_sim(n, maf1,maf2, beta, mu, s, causal,causal+noncausal)
    gen = Matrix(sim$gen,sparse=T)
    phen = sim$phen
    phen.perm = matrix(nrow=length(phen), ncol=permutes, replicate(permutes,sample(phen)))
    
        test.statistic = Qstat(gen, phen)
        #permutation test
        reps=Qstat(gen,phen.perm)
        ps[j]=mean(reps>=test.statistic)
    }
  pow.mat[k,]=ps
  
}
pow=vector("numeric") 
for (h in 1:length(pow.mat[,1])){
   pow=c(pow,mean(pow.mat[h,]<0.05))}


params=list(non.settings=non.settings, pow.trials=pow.trials, permutes=permutes, causal=causal,
n=n, beta=beta, maf1=maf1,maf2=maf2,s=s, mu=mu)

result=list(pow.mat=pow, params=params)

save(result, file="top10liuother.Rdata")

EOF