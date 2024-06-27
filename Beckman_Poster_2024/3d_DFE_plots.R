install.packages('mvtnorm')
install.packages('rgl')
library(mvtnorm)
library(rgl)

#establish parameter values
log_mu1 <- 4.70915535794278
log_sigma1 <- 6.234717665387251
log_mu2 <- 4.875076242644642
log_sigma2 <- 6.957650566833429
rho <- 0.9735646160233167

#convert parameter values to not log scale
mu1 <- exp(log_mu1)
sigma1 <- exp(log_sigma1)
mu2 <- exp(log_mu2)
sigma2 <- exp(log_sigma2)

#make grid
n <- 50
x <- seq(0, 3, length.out = n)
y <- seq(0, 3, length.out = n)
z <- matrix(0, n, n)

#calculate values
for (i in 1:n){
  for (j in 1:n){
    z[i, j] <- dmvnorm(c(x[i], y[j]), mean = c(mu1, mu2), sigma = matrix(c(sigma1^2, rho * sigma1 * sigma2, rho * sigma1 * sigma2, sigma2^2), nrow = 2))
  }
}

#plot the distribution
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue", xlab = "Mmd_IRA", ylab = "Mmd_FRA", zlab = "Density", main = "Asymmetric Bivariate Lognormal Distribution")
contour(x, y, z, add = TRUE, color = 'black')
