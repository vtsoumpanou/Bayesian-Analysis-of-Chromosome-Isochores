
#===============================================================================
#Further Analysis in the 3rd dataset from the Isochore Analysis project
#Model 2 from Isochore Analysis is assumed, thus the data come from 2 isochores
#===============================================================================

chr3 <-
  c(1776, 1929, 1885, 2058, 1907, 2143, 2206, 2112, 2131, 1948, 
    2291, 2112, 2122, 2985, 2690, 2634, 2339, 2924, 2805, 2755, 2901, 
    2527, 2159, 2337, 2237, 2224, 2251, 2384, 2340, 2234, 2642, 2422, 
    2381, 2527, 2428, 2167, 2275, 2589, 2364, 2407, 2736, 2452, 2600, 
    2554, 2640, 2396, 2488, 2660, 2512, 2746, 2559, 2399, 2770, 2870, 
    2300, 2391, 1936, 2308, 2128, 1940, 1925, 2108, 2343, 2223, 2174, 
    2277, 2480, 2749, 2794, 2846, 2792, 2668, 2712, 2793, 2209, 2381, 
    2218, 2368, 2503, 2619, 2378, 2262, 2145, 2008, 2186, 1994, 2075, 
    1933, 2214, 2408, 2049, 2110, 2391, 2218, 2337, 2808, 2244, 2213, 
    2569, 2455)

#===============================================================================
# 1. Priors speciffication for Model 2:
# x_1, x_2, ..., x_l, l=length(chr3) are the observations in each bade "window"
# Let t be the Change Point index, t in {1, 2, ..., l-1}. Then:
# x_1, ..., x_t ~ Binom(n=5000, p1), x_(t+1), ..., x_l ~ Binom(n, p2)
# where p1, p2: the prob. of C base. The parameters for this problem are:
# p1, p2, t. We assume: p1, p2 ~ Beta(0.5, 0,5) - Jeffrey's prior
# and t ~ Discrete_Uniform(1, l1), l1:=l-1
# The priors are chosen to be uninformative due to: 1. lack of expert knowledge
# 2. large volume of data ---> from asymptotics we know that the priors do not
# heavily influence the inference.
#===============================================================================
# 2. The joint posterior is (up to a constant):
# Beta(0.5+sum_{i=1}^{i=t}(x_i), 0.5+n-sum_{i=1}^{i=t}(x_i)) X
# Beta(0.5+sum_{i=t+1}^{i=l1}(x_i), 0.5+n-sum_{i=t+1}^{i=l1}(x_i)) X
# I[t in {1,...., l1}].
# The conditional posteriors are easily extracted from the joint
#===============================================================================
# 3. Gibbs sampler                          ====================================
#===============================================================================

gibbs <- function(data = chr3,
                  n = 5000,
                  a1 = 0.5,
                  b1 = 0.5,
                  a2 = 0.5,
                  b2 = 0.5,
                  nburn = 0,
                  ndraw = 5000) {
  # Gibbs sampling for Binomial change point problem
  l  <- length(data)
  l1 <- l - 1
  sum_total <- sum(data)
  cum_sumx  <- cumsum(data)
  
  #initial values:drawn from prior
  p1 <- rbeta(1, a1, b1)
  p2 <- rbeta(1, a2, b2)
  t  <- sample(1:l1, 1)
  
  # create matrix to record draws:
  draws <- matrix(NA, ndraw, 3)
  
  # MCMC LOOP FOLLOWS:
  
  iter <- -nburn
  while (iter < ndraw) {
    iter <- iter + 1
    
    sum_curr <- cum_sumx[t]
    
    # draw p1
    p1 <- rbeta(1, a1 + sum_curr, b1 + n * t - sum_curr)
    
    # draw p2
    p2 <- rbeta(1,
                a2 + (sum_total - sum_curr),
                b2 + n * (l - t) - (sum_total - sum_curr))
    
    
    logp <- numeric(l1)
    for (j in 1:l1) {
      sum_first <- cum_sumx[j]
      sum_second <- sum_total - sum_first
      
      logp[j] <- sum_first * log(p1) + (n * j - sum_first) * log(1 - p1) +
        sum_second * log(p2) + (n * (l - j) - sum_second) * log(1 - p2)
    }
    
    
    # log trick for numerical stability
    p <- exp(logp - max(logp))
    # normalize p to be a probability vector
    p <- p / sum(p)
    # sample t according to p
    t <- sample(1:l1, 1, prob = p)
    
    # record draws
    if (iter > 0)
      draws[iter, ] <- c(p1, p2, t)
  }
  
  draws
}

#===============================================================================
# 4. Plots and Inference                       =================================
#===============================================================================

# Convergence of the algorithm
draws <- gibbs(data = chr3, ndraw = 1000)
p1 <- draws[, 1]
p2 <- draws[, 2]
t <- draws[, 3]

par(mfrow = c(1, 3))

plot(p1, type = "l", main = "Draws p1")
plot(p2, type = "l", main = "Draws p2")
plot(t, type = "l", main = "Draws t")

# Histograms
par(mfrow = c(2, 2))

hist(t, probab = T)
hist(p1, probab = T)
hist(p2, probab = T)
hist(p1 - p2, probab = T)
