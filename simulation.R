nrow = 2500  # Number of sites
n.spp = 100L # Number of species
dim = 5      # Dimension of x (number of environmental variables)

logistic = binomial()$linkinv

# Bernoulli random number generator that outputs matrices instead of vectors
rbern = function(p){
  out = rbinom(length(p), size = 1, prob = p)
  dim(out) = dim(p)
  out
}

# x is filled with independent standard Gaussians
x = matrix(rnorm(nrow * dim, sd = 1), ncol = dim)

# true response to environment is also independent Gaussians, but SD is 1/2
w = matrix(rnorm(dim*n.spp, mean = 0, sd = 0.5), nrow = dim, ncol = n.spp)


# Initialize the simulated landscape with all zeros
y = matrix(
  0, 
  nrow = nrow, 
  ncol = n.spp
)

# All species' intercepts are zero
b = rep(0, n.spp)

# l, for "lateral" contains all negative values. Magnitude is exponentially
# distributed.  By convention, the diagonal is zero and the matrix is symmetric
l = matrix(0, nrow = n.spp, ncol = n.spp)
l[upper.tri(l)] = -rexp(choose(n.spp, 2), rate = 5)
l = l + t(l)


# 1000 rounds of Gibbs sampling (takes a minute or two)
for(i in 1:1E3){
  for(j in sample(1:n.spp)){
    # Update species j:
    y[ , j] = rbern(
      logistic(
        # intercept + environment (x * w) + species interactions(y * l)
        b[j] + x %*% w[,j] + y %*% l[, j]
      )
    )
  }
}
