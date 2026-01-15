#Load necessary libraries
library(lme4)
library(ggplot2)

# read in data
data <- read.delim('/Users/Natalie/Desktop/Nat_spz_data3.txt',header=T)
library(minpack.lm)

# Step 1: Define a better logistic model function
logistic_model <- function(log10conc, L, k, x0) {
  L / (1 + exp(-k * (log10conc - x0)))
}

# Step 2: Define a function to fit the logistic model using nlsLM (from minpack.lm package)
logistic_fit <- function(df) {
  # Calculate good starting values based on the data
  L_start <- max(df$motility, na.rm = TRUE)  # Maximum motility observed
  x0_start <- median(df$log10conc, na.rm = TRUE)  # Median log10conc as x0
  k_start <- 0.1  # Initial guess for the steepness
  
  # Fit the model with increased maxiter for more iterations
  fit <- tryCatch({
    nlsLM(motility ~ logistic_model(log10conc, L, k, x0), 
          data = df,
          start = list(L = L_start, k = k_start, x0 = x0_start),
          control = nls.lm.control(maxiter = 1024,maxfev=10000))  # Increase maxiter to 200
  }, error = function(e) {
    return(e)  # Return the error if the model fitting fails
  })
  return(fit)
}

time_points <- unique(data$time)
fit_results <- list()

for (t in time_points) {
  df_time <- subset(data, time == t)
  fit_results[[as.character(t)]] <- logistic_fit(df_time)
}

# Step 4: Check the results
fit_results


# Set specific colors for the time points
time_colors <- c("0" = "#1A85FF", "15" = "#D41159", "60" = "#5D3A9B")

# Compute mean and standard deviation for motility at each time point and concentration
summary_data <- data %>%
  group_by(time, log10conc) %>%
  summarize(
    mean_motility = mean(motility, na.rm = TRUE),
    sd_motility = sd(motility, na.rm = TRUE),
    n = n()
  ) %>%
  ungroup()

# Plot with fitted logistic lines and error bars
ggplot(data, aes(x = log10conc, y = motility, color = factor(time))) +
  stat_smooth(
    method = "nlsLM",
    formula = y ~ logistic_model(x, L, k, x0),
    method.args = list(start = list(L = max(data$motility), k = 0.1, x0 = median(data$log10conc))),
    se = FALSE,
    size = 1
  ) +
  geom_errorbar(
    data = summary_data,
    aes(
      x = log10conc,
      ymin = mean_motility - sd_motility,
      ymax = mean_motility + sd_motility,
      y = mean_motility
    ),
    width = 0.1,
    inherit.aes = FALSE,  # Prevent inheriting aesthetics from the main dataset
    color = "black",
    alpha = 0.8
  ) +
  geom_point(
    data = summary_data,
    aes(x = log10conc, y = mean_motility),
    inherit.aes = FALSE,
    size = 3,
    color = "black",
    shape = 21,
    fill = "white"
  ) +
  facet_wrap(~time, scales = "free_y") +
  labs(
    title = "Fitted Logistic Curves with Error Bars (Mean ± SD)",
    x = "Log10(Concentration)",
    y = "Motility",
    color = "Time Point"
  ) +
  theme_minimal() +
  scale_y_continuous(limits=c(0,2))+
  theme(legend.position = "bottom")+
  scale_color_manual(values = time_colors) +  # Apply custom colors for time
  labs(
    title = "Fitted Logistic Curves with Error Bars (Mean ± SD)",
    x = "Log10(Concentration)",
    y = "Motility",
    color = "Time Point"
  ) 


# Extract and display the coefficients (parameters L, k, x0) for each time point
fit_parameters <- lapply(fit_results, function(fit) {
  if (inherits(fit, "nls")) {
    coef(fit)
  } else {
    # If the model fitting fails, return the error message
    fit$message
  }
})

# Print the estimated parameters for each time point
fit_parameters

###Determine are if the lines are signficantly different from each other
set.seed(123)  # Set seed for reproducibility
n_bootstrap <- 50  # Number of bootstrap iterations
bootstrap_results <- data.frame(time = integer(), L = numeric(), k = numeric(), x0 = numeric())

for (t in unique(data$time)) {
  df_time <- subset(data, time == t)
  
  for (i in 1:n_bootstrap) {
    df_bootstrap <- df_time[sample(nrow(df_time), replace = TRUE), ]  # Resample with replacement
    fit <- logistic_fit(df_bootstrap)
    
    # Check if the model fit was successful
    if (!is.null(fit)) {
      params <- coef(fit)
      
      # Only add results if parameters are valid
      if (length(params) == 3) {
        bootstrap_results <- rbind(bootstrap_results, data.frame(time = t, L = params['L'], k = params['k'], x0 = params['x0']))
      }
    }
  }
}



# Now calculate confidence intervals
library(dplyr)

bootstrap_results = bootstrap_results %>% filter(L <=2)

# Calculate 95% confidence intervals for each parameter across time points
ci_results <- bootstrap_results %>%
  group_by(time) %>%
  summarize(
    L_lower = quantile(L, 0.025),
    L_upper = quantile(L, 0.975),
    k_lower = quantile(k, 0.025),
    k_upper = quantile(k, 0.975),
    x0_lower = quantile(x0, 0.025),
    x0_upper = quantile(x0, 0.975)
  )

print(ci_results)

# Kruskal-Wallis test for each parameter across time points
kruskal_L <- kruskal.test(L ~ time, data = bootstrap_results)
kruskal_k <- kruskal.test(k ~ time, data = bootstrap_results)
kruskal_x0 <- kruskal.test(x0 ~ time, data = bootstrap_results)

# Print the test results - x0 is the important once, it tells us that the concentration at which the slope is steepest differs, ie time dependency
kruskal_L
kruskal_k
kruskal_x0

# Perform Dunn's test for pairwise comparisons
library(FSA)

# Dunn's test for parameter L
dunn_L <- dunnTest(L ~ time, data = bootstrap_results, method = "bonferroni")

# Dunn's test for parameter k
dunn_k <- dunnTest(k ~ time, data = bootstrap_results, method = "bonferroni")

# Dunn's test for parameter x0
dunn_x0 <- dunnTest(x0 ~ time, data = bootstrap_results, method = "bonferroni")

# Print the Dunn's test results
dunn_L
dunn_k
dunn_x0

###Here, L is the upper value of the response (motility), k is growth rate (steepness of the curve) and x0 is the concentration at which
#the slope is steepest

# Extract parameters (L, k, x0) from each fit
parameter_estimates <- data.frame(time = integer(), L = numeric(), k = numeric(), x0 = numeric())

for (t in time_points) {
  fit <- fit_results[[as.character(t)]]
  
  # Check if the model was successfully fitted
  if (!inherits(fit, "error")) {
    params <- coef(fit)
    parameter_estimates <- rbind(parameter_estimates, data.frame(
      time = t,
      L = params["L"],
      k = params["k"],
      x0 = params["x0"]
    ))
  }
}





# Add predictions to the original data using the fitted parameters
data_with_predictions <- data.frame()

for (t in time_points) {
  df_time <- subset(data, time == t)
  fit <- fit_results[[as.character(t)]]
  
  if (inherits(fit, "nls")) {
    # Generate predicted motility based on the logistic model for each log10conc
    df_time$predicted_motility <- predict(fit, newdata = df_time)
    data_with_predictions <- rbind(data_with_predictions, df_time)
  }
}

# Ensure concentration and time are factors
data$concentration <- factor(data$concentration, levels = c(0, sort(unique(data$concentration[data$concentration != 0]))))
data$time <- factor(data$time)

# Create an empty list to store ANOVA results for each time point
anova_results_list <- list()
post_hoc_results_list <- list()

# Step 1: Loop over each time point
for(t in unique(data$time)) {
  # Subset the data for each time point
  data_time_point <- subset(data, time == t)
  
  # Step 2: Run ANOVA to test the effect of concentration within this time point
  anova_results <- aov(motility ~ concentration, data = data_time_point)
  
  # Store ANOVA results for the current time point
  anova_results_list[[as.character(t)]] <- summary(anova_results)
  
  # Step 3: Run Tukey's HSD test to compare concentrations within this time point
  post_hoc_results <- TukeyHSD(anova_results, "concentration")
  
  # Store post-hoc results for the current time point
  post_hoc_results_list[[as.character(t)]] <- post_hoc_results
}

# Step 4: Print the ANOVA and post-hoc results
anova_results_list
post_hoc_results_list



###This approach tests whether log10(concentration) has a significant effect on motility across time point - ie concentration as a numeric variable
# Loop over each time point
regression_results_list <- list()

for(t in unique(data$time)) {
  # Subset the data for each time point
  data_time_point <- subset(data, time == t)
  
  # Linear regression with concentration as a continuous variable
  regression_model <- lm(motility ~ log10conc, data = data_time_point)
  
  # Store regression results
  regression_results_list[[as.character(t)]] <- summary(regression_model)
}

# Print regression results
regression_results_list
