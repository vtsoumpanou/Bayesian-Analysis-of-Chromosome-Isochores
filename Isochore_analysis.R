#===============================================================================
# INTRODUCTION            ======================================================
#===============================================================================
# Analysis will be conducted on 5 sets of data from 5 chromosomes.
# The data are in the form of 100 counts of CG bases 
# out of 5000 base windows in chromosomes.
# The objective of this analysis is to make inferences to conclude whether 
# the data in each set come from one or two different isochores. 
# A Bayesian approach is followed, taking into account expert judgement.
#===============================================================================

# Get Data

chr1 <- c(2260, 2249, 2177, 1694, 2289, 2304, 2472, 2076, 2201, 1932, 
        2290, 2201, 2271, 2243, 2209, 2243, 2270, 2167, 2273, 2172, 2044, 
        2285, 2114, 1999, 2224, 2208, 2260, 2544, 2482, 2487, 2255, 2505, 
        2278, 2014, 2223, 2028, 2480, 2354, 2359, 2179, 2504, 1887, 2123, 
        2162, 2638, 2333, 2348, 2229, 2052, 2159, 2323, 2209, 2052, 1732, 
        2284, 2350, 1908, 2144, 2225, 2259, 2182, 2341, 2107, 2423, 2475, 
        2180, 2337, 2332, 2441, 2242, 2318, 2239, 2391, 2198, 1814, 2278, 
        2430, 2439, 2412, 2519, 2166, 2141, 2115, 2048, 2006, 1979, 1984, 
        2501, 1995, 2233, 2384, 1933, 2142, 2138, 2040, 2280, 2093, 2213, 
        2342, 2235)

chr2 <- c(2249, 2362, 2263, 2469, 2352, 2320, 2290, 2267, 2235, 2269, 
        2272, 2273, 2185, 2496, 2195, 2155, 1992, 2323, 2427, 2189, 2377, 
        2282, 2317, 2174, 2400, 2300, 2120, 2259, 2262, 2271, 2334, 2833, 
        2375, 2426, 2489, 2446, 2693, 2088, 2096, 2121, 2150, 2129, 2111, 
        2145, 2185, 1915, 2274, 2193, 2171, 2154, 2143, 2066, 2100, 2188, 
        2058, 2327, 2034, 2197, 2169, 2171, 2084, 2081, 2475, 2187, 2147, 
        1931, 2163, 1988, 1997, 2261, 1917, 1920, 2095, 2064, 2082, 1954, 
        2134, 1980, 1917, 1954, 1925, 2026, 2014, 2178, 2208, 2055, 2117, 
        2030, 1982, 2056, 2118, 2072, 2185, 2112, 2197, 2027, 2104, 2637, 
        2363, 2019)

chr3 <- c(1776, 1929, 1885, 2058, 1907, 2143, 2206, 2112, 2131, 1948, 
        2291, 2112, 2122, 2985, 2690, 2634, 2339, 2924, 2805, 2755, 2901, 
        2527, 2159, 2337, 2237, 2224, 2251, 2384, 2340, 2234, 2642, 2422, 
        2381, 2527, 2428, 2167, 2275, 2589, 2364, 2407, 2736, 2452, 2600, 
        2554, 2640, 2396, 2488, 2660, 2512, 2746, 2559, 2399, 2770, 2870, 
        2300, 2391, 1936, 2308, 2128, 1940, 1925, 2108, 2343, 2223, 2174, 
        2277, 2480, 2749, 2794, 2846, 2792, 2668, 2712, 2793, 2209, 2381, 
        2218, 2368, 2503, 2619, 2378, 2262, 2145, 2008, 2186, 1994, 2075, 
        1933, 2214, 2408, 2049, 2110, 2391, 2218, 2337, 2808, 2244, 2213, 
        2569, 2455)

chr4 <- c(2557, 2441, 2488, 2435, 2297, 2325, 2769, 2642, 2620, 2569, 
        2248, 2316, 2185, 2236, 2305, 2526, 1970, 2185, 1902, 2136, 2073, 
        2193, 2191, 2517, 2148, 2107, 2154, 1955, 2210, 2143, 2187, 2285, 
        2280, 2045, 2236, 1993, 1895, 2087, 1919, 2084, 2077, 2165, 2177, 
        2047, 2099, 2127, 1948, 2283, 1868, 2057, 2076, 1924, 2053, 1798, 
        1907, 1912, 2081, 2042, 1969, 1849, 1976, 2121, 1970, 2022, 2262, 
        2104, 2177, 2211, 1969, 1946, 2105, 1945, 2116, 2126, 2113, 2623, 
        2908, 2749, 3020, 2182, 2363, 1989, 2099, 2333, 2094, 2195, 2144, 
        2241, 2059, 2009, 2334, 2027, 2135, 2315, 2015, 2015, 1981, 1869, 
        2218, 1976)
  
chr5 <- c(2253, 2041, 2230, 2082, 1904, 2254, 2084, 1999, 2014, 2080, 
        2141, 2092, 2081, 1974, 1966, 2227, 2179, 2248, 2116, 2140, 2084, 
        2029, 2009, 2066, 1943, 2034, 2065, 2256, 2042, 2048, 2032, 2146, 
        2197, 2240, 2341, 2041, 2246, 2267, 1845, 2281, 2681, 2466, 2513, 
        2512, 2498, 2568, 2596, 2542, 2776, 2360, 2326, 2230, 2163, 2360, 
        2026, 2103, 2202, 1923, 1967, 2027, 2139, 2224, 2107, 2053, 1991, 
        1987, 2208, 2278, 2319, 2411, 2258, 2683, 3156, 2783, 2630, 2802, 
        2671, 2598, 2640, 2982, 2850, 2759, 2607, 2690, 2240, 2426, 2416, 
        2910, 2411, 2228, 2268, 2494, 2534, 2375, 2287, 2394, 2175, 2267, 
        2428, 2511)

#===============================================================================
# 1. Model Specification for the data             ==============================
#===============================================================================
# A Binomial Model is chosen for the data as the #of CG bases can be interpreted
# as #of successes in 5000 trials.
#===============================================================================
#===============================================================================
# 2. Model Specification
# - chosen prior is jeffreys prior for Beta - Beta(0.5, 0.5)
#===============================================================================
# M_1 - In the one isochore case the most fitting model is Binomial(n,theta)
# The term sum n choose xi is common for both models, thus will be omitted


Model_1 <- function(chromosome,
                    a = 0.5,
                    b = 0.5,
                    n = 5000) {
  l <- length(chromosome)
  S <- sum(chromosome)
  LogEv1 <- lbeta(a + S, b + (n * l) - S) - lbeta(a, b)
  
}

# M_2 - in the two isochores case the most fitting model is the mixture of two
# Binomials with n=5000 and theta_1 for i<=t,
# theta_2 for i>t the probabilities of CG in each isochore (t the change point)
# For the change point t we assume a Uniform(1,l)

Model_2 <- function(chromosome,
                    a = 0.5,
                    b = 0.5,
                    n = 5000) {
  l <- length(chromosome)
  l_1 <- l - 1
  
  S_total <- sum(chromosome)
  cum_chr <- cumsum(chromosome)
  
  logpr_t <- numeric(l_1)   # log p(data | t)
  
  for (t in 1:l_1) {
    # segment 1
    S1 <- cum_chr[t]
    T1 <- n * t
    
    # segment 2
    S2 <- S_total - S1
    T2 <- n * (l - t)
    
    
    logpr_t[t] <-
      (lbeta(a + S1, b + T1 - S1) - lbeta(a, b)) +
      (lbeta(a + S2, b + T2 - S2) - lbeta(a, b))
  }
  
  # Uniform prior for change point
  lprior_t <- log(1 / l_1)
  
  # log-sum-exp for numeric stability
  m <- max(logpr_t + lprior_t)  #the prior is uniform, does not affect the max
  LogEv2 <- m + log(sum(exp(logpr_t + lprior_t - m)))
  
}



#===============================================================================
# 3. Model Selection - Prior specification        ==============================
#===============================================================================
# The model prior chosen is Discrete_Uniform (1,2)


Model_Selection <- function(chromosome,
                            prior_M1 = 0.5,
                            prior_M2 = 0.5) {
  #Evidence for each Model
  LogEv1 <- Model_1(chromosome)
  LogEv2 <- Model_2(chromosome)
  
  #Calculation of Posterior Model probabilities
  Post_prob_M1 <- 1 / (1 + exp(LogEv2 - LogEv1))
  Post_prob_M2 <- 1 - Post_prob_M1
  
  inference <- sprintf(
    "Model 1 has posterior probability %.3f, while Model 2 has %.3f",
    Post_prob_M1,
    Post_prob_M2
  )
  inference
}


#===============================================================================
# 4. Results - Inferences              =========================================
#===============================================================================

data <- list(chr1, chr2, chr3, chr4, chr5)

lapply(data, Model_Selection)


#===============================================================================
# Further Analysis                     =========================================
#===============================================================================

par(mfcol = c(1, 5))
plots <- function(chromosome) {
  CG_content <- chromosome / 5000
  
  plot(CG_content,
       type = "l",
       xlab = "Window",
       ylab = "CG Content")
}



lapply(data, plots)




# Bayes factor - calculation of the logarithm for numeric stability

log_BF <- function(chromosome) {
  log_BF <- Model_2(chromosome) - Model_1(chromosome)
  
}

lapply(data, log_BF)




# Data histogram

data_hist <- function(chromosome) {
  hist(
    chromosome / 5000,
    main = "CG Content Distribution",
    xlab = "CG Proportion",
    col = "lightgrey",
    breaks = 20
  )
  lines(density(chromosome / 5000, adjust = 0.5),
        col = 'red',
        lwd = 2)
}

lapply(data, data_hist)