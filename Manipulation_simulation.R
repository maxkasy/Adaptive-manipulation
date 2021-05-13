library(tidyverse)
library(furrr)
source("Bayesian_optimization_routines.R")
source("DGP_feature_manipulation.R")

FUN = lossFreshIndiv # function used to sample loss - based on the feature manipulation DGP
K = 4 # dimension of argument of loss

replications = 32 # number of simulation replications

nEval = 300 # number of simulation draws
nStart = 20 # size of warm-up sample (number of randomly drawn observations)
Beta = matrix(runif(nStart * 4), nStart, 4) # rows corresponds to the arguments X that we are optimizing over
Loss = matrix(apply(Beta, 1, FUN), nStart, 1) # this vector corresponds to Y - the loss that we are minimizing

theta = matrix(rep(.5, K + 2), K + 2, 1) # initial values for hyper-parameters
theta = maximize_marginal_LLH_update(X = Beta, Y = Loss, theta) # tune hyper-parameters using the warm-up sample
print(theta)

betamin = rep(0, K) # lower and upper bound for the coefficients
betamax = rep(1, K)
nFeatures = 1000 # number of features used for spectral-sampling approximation to the GP posterior
M = 5 # number of restarts for optimization in Thompson sampling step


# run simulations repatedly
plan(multisession)
sim_data = future_map(1:replications,
                        ~ BayesianOptimization(FUN, nEval, 
                                X = Beta, Y = Loss, 
                                theta, betamin, betamax, nFeatures, M))

# combine loss from simulation replications into a matrix
loss_matrix = 
    do.call(cbind,
    sim_data %>% 
        map("Y") %>% 
        map(as.vector)) 

loss_matrix %>% 
    as_tibble() %>% 
    write_csv("Simulation_output/loss_matrix_simulated.csv")

# data-frame of average loss across simulations, for each time period
average_loss = tibble(
    Period = 1:nrow(loss_matrix),
    Loss = rowMeans(loss_matrix) 
)
   
# plot results
average_loss %>% 
    filter(Period > nStart) %>%
    ggplot(aes(x= Period, y = Loss)) +
    geom_line() +
    geom_smooth() +
    labs(title = "Average simulated loss",
         caption = paste0("Average loss across ", replications, 
                          " simulations. Dimension of regressors: ", K, 
                          ". Blue line shows smoothed average loss.")) +
    theme_light()

ggsave("Simulation_output/average_loss_simulated.png",
       width =7, height = 4)
