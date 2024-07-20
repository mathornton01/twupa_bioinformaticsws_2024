# Load necessary libraries
library(fitdistrplus)
library(ggplot2)

# Parameters
n <- 10000  # number of trials
k <- 4  # number of categories
p <- rep(1/k, k)  # uniform probabilities for each outcome
num_samples <- 10000  # number of samples for the simulation

# Function to generate samples and compute the minimum values
simulate_minimum_multinomial <- function(n, p, num_samples) {
  min_values <- numeric(num_samples)
  for (i in 1:num_samples) {
    sample <- rmultinom(1, n, p)
    min_values[i] <- min(sample)
  }
  return(min_values)
}

# Generate the samples
sim_min_values <- simulate_minimum_multinomial(n, p, num_samples)

# Scale the minimum values by the number of trials
scaled_min_values <- sim_min_values / n

# Ensure values are within the interval [0, 1]
scaled_min_values <- pmin(pmax(scaled_min_values, 0), 1)

# Fit candidate distributions to the scaled sample data
fit_gamma <- fitdist(scaled_min_values, "gamma")
fit_normal <- fitdist(scaled_min_values, "norm")
fit_beta <- fitdist(scaled_min_values, "beta")
fit_weibull <- fitdist(scaled_min_values, "weibull")

# Compare the fits using AIC
fits <- list(gamma = fit_gamma, normal = fit_normal, beta = fit_beta, weibull = fit_weibull)
aic_values <- sapply(fits, function(f) f$aic)
best_fit <- names(which.min(aic_values))

# Print the AIC values and the best-fit distribution
print(aic_values)
print(paste("Best fit based on AIC:", best_fit))

# Create a data frame for the empirical density function
empirical_density <- density(scaled_min_values)
empirical_df <- data.frame(x = empirical_density$x, y = empirical_density$y)

# Function to calculate the PDF of the gamma distribution with estimated parameters
pdf_gamma <- function(x) {
  dgamma(x, shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"])
}

# Function to calculate the PDF of the normal distribution with estimated parameters
pdf_normal <- function(x) {
  dnorm(x, mean = fit_normal$estimate["mean"], sd = fit_normal$estimate["sd"])
}

# Function to calculate the PDF of the beta distribution with estimated parameters
pdf_beta <- function(x) {
  dbeta(x, shape1 = fit_beta$estimate["shape1"], shape2 = fit_beta$estimate["shape2"])
}

# Function to calculate the PDF of the weibull distribution with estimated parameters
pdf_weibull <- function(x) {
  dweibull(x, shape = fit_weibull$estimate["shape"], scale = fit_weibull$estimate["scale"])
}

# Function to calculate the PDF of the best-fit distribution
pdf_best_fit <- function(x) {
  if (best_fit == "gamma") {
    pdf_gamma(x)
  } else if (best_fit == "normal") {
    pdf_normal(x)
  } else if (best_fit == "beta") {
    pdf_beta(x)
  } else if (best_fit == "weibull") {
    pdf_weibull(x)
  }
}

# Create data frames for all PDFs within a focused range around the distributions
x_vals <- seq(0.2, 0.3, length.out = 100)
gamma_pdf_df <- data.frame(x = x_vals, y = pdf_gamma(x_vals), distribution = "Gamma")
normal_pdf_df <- data.frame(x = x_vals, y = pdf_normal(x_vals), distribution = "Normal")
beta_pdf_df <- data.frame(x = x_vals, y = pdf_beta(x_vals), distribution = "Beta")
weibull_pdf_df <- data.frame(x = x_vals, y = pdf_weibull(x_vals), distribution = "Weibull")

# Combine all PDF data frames
all_pdfs_df <- rbind(gamma_pdf_df, normal_pdf_df, beta_pdf_df, weibull_pdf_df)

# Create a data frame for the best-fit PDF
best_fit_pdf_df <- data.frame(x = x_vals, y = sapply(x_vals, pdf_best_fit))

# Plot the empirical density, the PDFs of all fitted distributions, and the best-fit PDF
ggplot() +
  geom_line(data = empirical_df, aes(x = x, y = y, color = "Empirical Density"), size = 1) +
  geom_line(data = best_fit_pdf_df, aes(x = x, y = y, color = "Best-fit Distribution PDF"), size = 1) +
  geom_line(data = all_pdfs_df, aes(x = x, y = y, color = distribution), linetype = "dashed") +
  scale_color_manual(values = c("Empirical Density" = "blue", "Best-fit Distribution PDF" = "red",
                                "Gamma" = "green", "Normal" = "purple", "Beta" = "orange", "Weibull" = "brown")) +
  labs(x = "Scaled Minimum Value", y = "Density",
       title = paste("Empirical Density and Fits for Scaled Minimum Value\n(Best fit based on AIC:", best_fit, ")")) +
  theme_minimal() +
  ggtitle("Empirical Density and Fitted Distributions") +
  theme(legend.position = "right") +
  xlim(0.225, 0.250) + ylim(0,200)

