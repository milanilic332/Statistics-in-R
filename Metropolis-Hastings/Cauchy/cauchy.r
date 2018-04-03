# Multivariate Cauchy distribution 
Cauchy <- function(x, x0 = c(50, 50), g = 1){
  return(log((1/(2*pi))*(g/((x[1] - x0[1])^2 + (x[2] - x0[2])^2 + g^2)^1.5)))
}

# Proposal distribution - asymmetric
proposalDist <- function(param){
  return(rnorm(2, mean = param + 1, sd = 5))
}

q <- function(x, xP){
  x1 = dnorm(x[1], mean = xP[1] + 1, sd = 5)
  x2 = dnorm(x[2], mean = xP[2] + 1, sd = 5)
  return(log(x1*x2))
}

# Metropolis-Hastings algorithm
MA <- function(startValue, iterations){
  chain = matrix(data = NA, ncol = 2, nrow = iterations + 1)
  chain[1, ] = startValue

  for(i in 1:iterations){
    
    proposal = proposalDist(chain[i, ])
    probab = exp(Cauchy(proposal) + q(chain[i, ], proposal) - Cauchy(chain[i, ]) - q(proposal, chain[i, ]))
    
    if(runif(1) < probab){
      chain[i + 1, ] = proposal
    }
    else{
      chain[i + 1, ] = chain[i, ]
    }
    
  }
  return(chain)
}

# Start values
startValue = c(10, 10)
iterations = 50000

# Run MA
chain = MA(startValue, iterations)

# Burn first 5000 values
burnIn = 5000

# Histograms (with thinning)
par(mfrow = c(1, 2))

n1 = length(chain[seq(from = burnIn, to = iterations, by = 10), 1])
k1 = floor(log(n1, base = 2)) + 1
d1 = diff(range(chain[seq(from = burnIn, to = iterations, by = 10), 1]))/k1

n2 = length(chain[seq(from = burnIn, to = iterations, by = 10), 2])
k2 = floor(log(n2, base = 2)) + 1
d2 = diff(range(chain[seq(from = burnIn, to = iterations, by = 10), 2]))/k2

sortedChain = matrix(data = NA, ncol = 2, nrow = n1)

sortedChain[, 1] = sort(chain[seq(from = burnIn, to = iterations, by = 10), 1])
sortedChain[, 2] = sort(chain[seq(from = burnIn, to = iterations, by = 10), 2])

br1 = sortedChain[, 1][1] + 0:k1*d1
br2 = sortedChain[, 2][1] + 0:k2*d2

hist(sortedChain[, 1], breaks = br1, main = "Distribution of x[1]")
hist(sortedChain[, 2], breaks = br2, main = "Distribution of x[2]")

# Plot chains
par(mfrow = c(1, 2))
plot(chain[seq(from = burnIn, to = iterations, by = 10), 1], type = 'l')
plot(chain[seq(from = burnIn, to = iterations, by = 10), 2], type = 'l')

# Plot chains together
par(mfrow = c(1, 1))
plot(chain[seq(from = burnIn, to = iterations, by = 10), 1], chain[seq(from = burnIn, to = iterations, by = 10), 2], type = 'l', xlab = "x[2] values", ylab = "x[1] values", col = "coral1")