# Load required libraries
library(MASS)   # for generating multivariate normal distribution
library(ggplot2) # for plotting

# Function to generate bivariate lognormal data
generate_bivariate_lognormal <- function(n, logmu1, logmu2, logsd1, logsd2, rho) {
  # Transform log parameters to original scale
  mu1 <- exp(logmu1)
  mu2 <- exp(logmu2)
  sd1 <- exp(logsd1)
  sd2 <- exp(logsd2)
  
  # Generate correlated multivariate normal data in log-space
  Sigma <- matrix(c(logsd1^2, rho * logsd1 * logsd2, rho * logsd1 * logsd2, logsd2^2), nrow = 2, ncol = 2)
  Mu <- c(log(mu1), log(mu2))
  X <- mvrnorm(n, Mu, Sigma)
  
  # Transform to lognormal space
  X1 <- exp(X[,1])
  X2 <- exp(X[,2])
  
  return(data.frame(X1 = X1, X2 = X2))
}

# Parameters for the distribution (using logged parameters)
n <- 10000000
logmu1 <- 4.70915535794278  # log(mu1)
logmu2 <- 4.875076242644642  # log(mu2)
logsd1 <- 6.234717665387251  # log(sd1)
logsd2 <- 6.957650566833429  # log(sd2)
rho <- 0.9735646160233167

# Generate data
df <- generate_bivariate_lognormal(n, logmu1, logmu2, logsd1, logsd2, rho)



ggplot(df, aes(x = X1, y = X2)) + geom_bin2d(bins = 20) +  scale_fill_gradientn(colours = c("lightgrey", "blue")) +  labs(x = "X1", y = "X2", title = "IRA_FRA_plots.R") + theme_minimal() + scale_x_continuous(trans = "log10", limits = c(1e-9, 1e14), breaks = c(1e0, 1e14)) + scale_y_continuous(trans = "log10", limits = c(1e-9, 1e14), breaks = c( 1e0, 1e14))

