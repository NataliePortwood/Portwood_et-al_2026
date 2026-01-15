
# Load necessary libraries

library(MASS)  # For negative binomial regression

library(pscl)  # For zero-inflated models

#you would have an excel file with ‘oocysts’, ‘day’,’treatment’ as column headers
#Call your treatments treatment and control because then alphabetically everything is compared to control. 
#You want day to be categorical, so call it day 9, day 12 etc.

# Example dataset: replace with your actual data

#data <- read.delim("C:/Users/Natalie Portwood/Desktop",header=T)
data <- Oocyst_berghei_data_Binomial



# Negative Binomial Regression (NB)

nb_model <- glm.nb(oocysts ~ treatment + day, data = data)

summary(nb_model)



library(emmeans)

###Multiple comparisons

emm = emmeans(nb_model, pairwise ~ treatment | day, adjust = "Tukey")

emm

### to do this per day 

results <- list()

for (d in unique(data$day)) {
  data_sub <- subset(data, day == d)
  nb_model <- glm.nb(oocysts ~ treatment, data = data_sub)
  emm <- emmeans(nb_model, pairwise ~ treatment, adjust = "Tukey")
  results[[as.character(d)]] <- list(model = nb_model, emmeans = emm)
  
  cat("Day:", d, "\n")
  print(emm$contrasts)
}

