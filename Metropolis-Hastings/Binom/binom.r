trueP = 0.4

# Making test data
experiments <- function(n, prob){
  return(rbinom(n, 10, prob))
}
x = experiments(45, trueP)

# Prior distribution 
prior <- function(param){
  return(dbeta(param, 1, 1, log = T))  
}

# Likelihood distribution
likelihood <- function(param){
  singleLikelihoods = dbinom(x, 10, param, log = T)
  return(sum(singleLikelihoods))
}
  
# Posterior distribution
posterior <- function(param){
  return(prior(param) + likelihood(param))
}

# Proposal distribution - (E(X) is param)
proposalDist <- function(param){
  return(rbeta(1, 1, (1 - param)/param))
}

# E(X) is param
q <- function(param, paramP){
  return(log(dbeta(param, 1, (1 - paramP)/paramP)))
}
  
# Metropolis-Hastings algorithm
MHA <- function(startValue, iterations){
  chain = array(dim = iterations + 1)
  chain[1] = startValue
  
  for(i in 1:iterations){
    
    proposal = proposalDist(chain[i])
    probab = exp(posterior(proposal) + q(chain[i], proposal) - posterior(chain[i]) - q(proposal, chain[i]))
    
    if(runif(1) < probab){
      chain[i + 1] = proposal
    }
    else{
      chain[i + 1] = chain[i]
    }
    
  }  
  return(chain)
}

# Starting values
startValue = 0.5
iterations = 50000

# Run MHA
chain = MHA(startValue, iterations)

# Burn first 1000 iterations for better results
burnIn = 1000

par(mfrow = c(1, 2))

# Histogram (with thinning)
hist(chain[seq(from = burnIn, to = iterations, by = 10)], main = "Histogram", breaks = seq(from = 0, to = 1, by = 0.01), xlab = "(true value, mean) = (red, blue)")
abline(v = mean(chain), col = "blue")
abline(v = trueP, col = "red")

# Plot chain
plot(chain[seq(from = burnIn, to = iterations, by = 10)], type = "l", xlab="(true value, mean) = (red, blue)" , main = "Chain values of p" )
abline(h = mean(chain), col = "blue")
abline(h = trueP, col = "red")