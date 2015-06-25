# Power Calculations for Survey 4 

library(pwr)
library(dplyr)
library(reshape2)
library(ggplot2)

# Calculations
n_size <- seq(from=800, to=1200, by=20)
p_0.2 <- ((n_size/2) %>% pwr.t.test(d=0.2, sig.level=.05, type="two.sample"))$power
p_0.25 <- ((n_size/2) %>% pwr.t.test(d=0.25, sig.level=.05, type="two.sample"))$power
p_0.3 <- ((n_size/2) %>% pwr.t.test(d=0.3, sig.level=.05, type="two.sample"))$power
power_calculator <- function(mu_t, mu_c, sigma, alpha=0.05, N){
  lowertail <- (abs(mu_t - mu_c)*sqrt(N))/(2*sigma)
  uppertail <- -1*lowertail
  beta <- pnorm(lowertail- qnorm(1-alpha/2), lower.tail=TRUE) + 1-  pnorm(uppertail- qnorm(1-alpha/2), lower.tail=FALSE)
  return(beta)
}

# Simulations
power_sim <- function(possible.ns, alpha, sims, tau){
  powers <- rep(NA, length(possible.ns))           # Empty object to collect simulation estimates
  #### Outer loop to vary the number of subjects ####
  for (j in 1:length(possible.ns)){
    N <- possible.ns[j]                              # Pick the jth value for N
    
    significant.experiments <- rep(NA, sims)         # Empty object to count significant experiments
    
    #### Inner loop to conduct experiments "sims" times over for each N ####
    for (i in 1:sims){
      Y0 <-  rnorm(n=N, mean=0, sd=1)              # control potential outcome
      Y1 <- Y0 + tau                                 # treatment potential outcome
      Z.sim <- rbinom(n=N, size=1, prob=.5)          # Do a random assignment
      Y.sim <- Y1*Z.sim + Y0*(1-Z.sim)               # Reveal outcomes according to assignment
      fit.sim <- lm(Y.sim ~ Z.sim)                   # Do analysis (Simple regression)
      p.value <- summary(fit.sim)$coefficients[2,4]  # Extract p-values
      significant.experiments[i] <- (p.value <= alpha) # Determine significance according to p <= 0.05
    }
    
    powers[j] <- mean(significant.experiments)       # store average success rate (power) for each N
  }
  return(powers)
}
possible.ns <- seq(from=800, to=1200, by=20)     # The sample sizes we'll be considering
alpha <- 0.05                                    # Standard significance level
sims <- 1000                                 # Number of simulations to conduct for each N
set.seed(20349)
sim_0.2 <- power_sim(possible.ns = possible.ns, alpha = alpha, sims = 1000, tau = 0.2)
set.seed(50303)
sim_0.25 <- power_sim(possible.ns = possible.ns, alpha = alpha, sims = 1000, tau = 0.25)
set.seed(99202)
sim_0.3 <- power_sim(possible.ns = possible.ns, alpha = alpha, sims = 1000, tau = 0.3)

# plot a graph
pp_1 <- data.frame(n_size, es_0.2 = p_0.2, es_0.25 = p_0.25, 
                   es_0.3 = p_0.3, type = "Calculated")
pp_2 <- data.frame(n_size, es_0.2 = sim_0.2, es_0.25 = sim_0.25, 
                   es_0.3 = sim_0.3, type = "Simulated")
pp <- melt(data = rbind(pp_1, pp_2), id = c("n_size", "type"))
levels(pp$variable) <- c("0.2", "0.25", "0.3")
ggplot(data = pp, mapping = aes(x = n_size, y = value,
                                color = variable,
                                shape = variable)) + 
  geom_point() + xlab("N") + ylab("Power: Significance Level = 0.05") + 
  theme_bw() + 
  scale_shape_discrete(name=expression(tau)) + 
  scale_color_discrete(name=expression(tau)) + facet_grid(. ~ type)
ggsave("~/Dropbox/confounding/preanalysis/images/power_graph.pdf", 
       width = 8, height = 3.5)

