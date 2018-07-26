
library(Matrix)
source("continuous_simulator")


#define Qstat
Qstat = function(gen, phen.perm, p){
  mu.hat = rep(mean(phen), length(phen))
  S = crossprod(gen, phen.perm-mu.hat)
  Sp = abs(S)^p
  colSums(Sp)^(1/p)
}




non.settings=c(10,50,100,500)
pow.trials = 10000
p=c(3,5)
permutes = 1000
causal = 10 #number of causal variants
n = 1000 #number of individuals
beta = rep(0.75, causal) #effect size
maf=0.01
s = 1 #model standard deviation
mu = 0 #model mean

pow.mat = matrix(nrow=length(non.settings), ncol=length(p))

for(k in 1:length(non.settings))
  
{
  
  noncausal= non.settings[k] #number of noncausal variants
  print(noncausal)
  ps = matrix(ncol=length(p),nrow=pow.trials)
  
  for(j in 1:pow.trials){
    #simulate data
    sim = cont_sim(n, maf, beta, mu, s, causal,causal+noncausal)
    gen = Matrix(sim$gen,sparse=T)
    phen = sim$phen
    phen.perm = matrix(nrow=length(phen), ncol=permutes, replicate(permutes,sample(phen)))
    for(m in 1:length(p)){
        test.statistic = Qstat(gen, phen,p[m])
        #permutation test
        reps=Qstat(gen,phen.perm,p[m])
        ps[j,m]=mean(reps>=test.statistic)
    }
  }
  for(r in 1:length(p)){
    pow.mat[k,r] = mean(ps[,r]<.05)
  }
}

colnames(pow.mat) = sapply(p, toString)
row.names(pow.mat) = sapply(non.settings, toString)

params=list(non.settings=non.settings, pow.trials=pow.trials, p=p, permutes=permutes, causal=causal,
n=n, beta=beta, maf=maf, s=s, mu=mu)

result=list(pow.mat=pow.mat, params=params)

save(result, file="oddother.Rdata")

EOF