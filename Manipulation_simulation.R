library(tidyverse)
library(furrr)
source("Bayesian_optimization_routines.R")
source("DGP_feature_manipulation.R")

FUN = lossFreshIndiv
K = 4 # dimension of argument of loss

replications = 32

nEval = 300
nStart = 20
Beta = matrix(runif(nStart * 4), nStart, 4) # rows corresponds to X we are optimizing over
Loss = matrix(apply(Beta, 1, FUN), nStart, 1) # corresponds to Y - loss we are minimizing


theta = matrix(rep(.5, K + 2), K + 2, 1) # hyper-parameters
theta = maximize_marginal_LLH_update(X = Beta, Y = Loss, theta)

print(theta)
xmin = rep(0, K)
xmax = rep(1, K)
nFeatures = 1000
M = 5


# run simulations repatedly
plan(multisession)
sim_data = future_map(1:replications,
                        ~ BayesianOptimization(FUN, nEval, 
                                X = Beta, Y = Loss, 
                                theta, xmin, xmax, nFeatures, M))

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

ggsave("Simulation_output/average_lss_simulated.png",
       width =7, height = 4)
