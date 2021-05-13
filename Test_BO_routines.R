library(tidyverse)
# library(furrr)
source("Bayesian_optimization_routines.R")

FUN = function(x) {
    2 * (x - .5) ^ 2 + cos(10 * x) + .5 * rnorm(1)
}


# testing GP regression
nStart = 100
X = matrix(runif(nStart), nStart, 1)
Y = matrix(map_dbl(X, FUN), nrow(X), 1)
theta = matrix(rep(.5,3),3,1)
# theta = c(.1, 2, .05) # eyeballing solution
theta = maximize_marginal_LLH_update(X, Y, theta)


# plan(multisession)
plot_llh = expand.grid(
    theta1 = seq(0.01, 1, length.out = 20),
    theta2 = seq(0.01, 3, length.out = 20)
) %>%
    mutate(llh = map2_dbl(theta1, theta2, ~ marginal_llh(X, Y, c(.x, .y, .09))))

plot_llh %>%
    ggplot(aes(x = theta1, y = theta2, z = llh)) +
    geom_contour_filled() +
    theme_minimal() ->
    p1


data = tibble(
    y = as.vector(Y),
    x = as.vector(X),
    yhat = PosteriorMean(x, x, y, theta),
    e = y - yhat
)

p2 = data %>%
    ggplot(aes(x = x, y = y)) +
    geom_line(aes(y = yhat)) +
    geom_point() +
    theme_minimal()


profile_llh = tibble(theta2 = seq(0.01, 4, length.out = 40),
                     llh = map_dbl(theta2, ~ marginal_llh(X, Y, c(.13, .x, .09))))

profile_llh %>%
    ggplot(aes(x = theta2, y = llh)) +
    geom_line() +
    theme_minimal()



# testing bayesian optimization
nEval = 149
nStart = 20
X = matrix(runif(nStart), nStart, 1)
Y = matrix(map_dbl(X, FUN), nrow(X), 1)
theta = matrix(rep(.5, 3), 3, 1)
theta = maximize_marginal_LLH_update(X, Y, theta)
print(theta)
xmin = 0
xmax = 1
nFeatures = 1000
M = 5


tmp = BayesianOptimization(FUN, nEval, X, Y, theta, xmin, xmax, nFeatures, M)

data = tibble(
    y = as.vector(tmp$Y),
    x = as.vector(tmp$X),
    yhat = PosteriorMean(x, x, y, tmp$theta_opt),
    e = y - yhat
)


print(tmp[3:5])

print(mean(data$y))
print(var(data$y))
print(var(data$e))

p = data %>%
    ggplot(aes(x = x, y = y)) +
    geom_line(aes(y = yhat)) +
    geom_point() +
    theme_minimal()

hist(data$x, nclass = 50)

p
