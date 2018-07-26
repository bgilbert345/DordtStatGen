## cont_sim
## Purpose: Simulates a group of individuals with some number of alleles of interest
## and a continuous phenotype

## Output
## Returns a matrix of minor allele counts and a continuous phenotype for each individual
## Each row contains represents 1 individual
## Each column contains minor allele counts, except for the last column, which contains the phenotype

## Inputs
## n: a real number representing the number of simulations desired
## maf: a vector (length l) of minor allele frequencies
## betas: a vector (length l) of betas, which is the effect slope of minor alleles
## mu: a real number representing base phenotype level
## v: a real number representing the variance of the error distribution for the betas

cont_sim <- function(n, maf, betas, mu, s, causal,l) {
  ## G: a matrix (size n*l) that contains the number of minor alleles for each allele
  ## for each individual and their phenotype
  G <- matrix(nrow=n,ncol=l)
  ## phen: a vector (length n) that contains the simulated phenotypes
  phen = vector(length=n)
  ## each iteration of the loop simulates 1 individual, the minor allele counts are simulated
  ## by sampling from a binomial distribution, phenotype variation is sampled from a normal distribution
  for(i in 1:n) {
    G[i,] = rbinom(l, 2, maf)
    
    phen[i] = mu+betas%*%G[i,1:causal]+rnorm(1,0, s)
  }
  ## Returns the genotype matrix and the phenotype vector
  list(gen= G, phen = phen)
}
