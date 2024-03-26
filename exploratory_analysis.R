### Exploratory analyses for the Biology 306 flour beetle population growth project

# Packages needed:
library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Dataset needed:
beetle.full <- read.csv("beetle_data_FULL.csv", header = TRUE)

#### STEP 1: Data wrangling ####

# Adding column for temperature
beetle.full$temperature <- sapply(str_split(beetle.full$code, "-"), function(x) tail(x, n = 1))

# Creating a new data frame for total count
beetle <- beetle.full %>%
  group_by(code, week, temperature) %>%
  summarise(total.count = sum(count))

# Calculating mean final density
beetle.max.week <- beetle %>%
  group_by(code) %>%
  filter(week == max(week)) %>%
  group_by(temperature) %>%
  summarise(mean.total.count = mean(total.count))

#### STEP 2: Visualizations ####

#X <- 35 # Swap out for different temperatures
#beetle.subset <- beetle %>%
#  filter(temperature == X)

palette <- rev(brewer.pal(n = 6, name = "RdBu"))

(gg.all <- ggplot(aes(x = week, y = total.count, color = temperature, group = code), 
                  data = beetle) +
  scale_color_manual(values = palette) +
  facet_wrap(~ temperature) +
  theme_classic() +
  geom_line() +
  geom_point() +
  labs(x = "Week", y = "Total count", color = "Temperature (°C)"))
ggsave("all_temps.png", plot = gg.all, width = 10, height = 5)
  

#### STEP 3: Fitting models ####

# Define logistic function (see https://eligurarie.github.io/EFB370/labs/lab6/Lab6_FittingLogisticCurves.html)
N.logistic <- function(x, N0, K, r0) K/(1 + ((K - N0)/N0) * exp(-r0*x))

beetle.subset <- subset(beetle, code == levels(factor(beetle$code))[14]) # specify temp and population
plot(total.count ~ week, data = beetle.subset) #visualize

logistic.fit <- nls(total.count ~ N.logistic(week, N0, K, r0), 
                    data = beetle.subset,
                    start = list(N0 = 5, K = 50, r0 = 1))
# Datasets throwing error: 1, 2, 5, 7, 8, 11, 12, 14, 16, 17, 18, 19, 22, 23, 24, 25, 26, 27
summary(logistic.fit)







#### STOP!!! ####


# The following is a subset of code from Sean, to be tweaked and added in. 
# See original file for appropriate context (Lab_7_Analysis)



#-- Do nls() fits of logistic to time series data (mean_count vs. lab_week)
# Define logistic function (see https://eligurarie.github.io/EFB370/labs/lab6/Lab6_FittingLogisticCurves.html)
N.logistic <- function(x, N0, K, r0) K/(1 + ((K - N0)/N0)*exp(-r0*x))
#--Low temperature [fitted estimates are r0 = 1.33, K = 79.09]
fit_low <- nls(mean_count ~ N.logistic(lab_week, N0, K, r0), data = subset(pooled_average, pooled_average$treatment=='low'), 
               start = list(N0 = 5, K = 60, r0 = 1.5))
summary(fit_low)
#--Med temperature [fitted estimates are r0 = 1.63, K = 57.63]
fit_med <- nls(mean_count ~ N.logistic(lab_week, N0, K, r0), data = subset(pooled_average, pooled_average$treatment=='med'), 
               start = list(N0 = 5, K = 60, r0 = 1.5))
summary(fit_med)
#--High temperature [fitted estimates are r0 = 1.65, K = 35.94]
fit_high <- nls(mean_count ~ N.logistic(lab_week, N0, K, r0), data = subset(pooled_average, pooled_average$treatment=='high'), 
                start = list(N0 = 5, K = 60, r0 = 1.5))
summary(fit_high)

# Make loop to fit model and store parameters in a dataframe (https://stackoverflow.com/questions/71131348/extract-model-parameters-from-nls-model-for-loop)
library(tidyverse)
temp_treatment <- unique(pooled_average$treatment)

# define start to obtain the number and name of fit parameters
start <- list(N0 = 5, K = 60, r0 = 1.5)
# create empty data.frame to store IDs and parameters
params <- data.frame(matrix(nrow = length(temp_treatment), ncol = 1+length(start)))
names(params) <- c("temp_treatment",names(start))

#for loop to run through all the treatments
# (changed to 'seq_along()' to get consecutive 'i's
#  for easier indexing of the 'param'-data.frame)
for (i in seq_along(temp_treatment)){
  #Creating our dataframe for treatment "i"
  Individual_DFs <- pooled_average %>% filter (treatment %in% temp_treatment[i]) 
  
  #Fit model for treatment "i"
  nls.floop <- nls(mean_count ~ N.logistic(lab_week, N0, K, r0),
                   data = Individual_DFs, 
                   start = start)
  
  # store IDs
  params[i,1] <- temp_treatment[i]
  # store fit parameters
  params[i,2:ncol(params)] <- nls.floop$m$getPars()
}

#--Make a quick plot of averaged data, with fitted logistic models (https://stackoverflow.com/questions/46800516/how-to-plot-the-output-from-an-nls-model-fit-in-ggplot2)
plot(pooled_average$lab_week, pooled_average$mean_count,
     pch = 19,
     col = factor(pooled_average$treatment))
library(ggplot2)
ggplot(pooled_average, aes(x=lab_week, y=mean_count, color=treatment)) + 
  geom_point(size=3) + 
  geom_smooth(method = "nls", method.args = list(formula = y ~ K/(1 + ((K - N0)/N0)*exp(-r0*x)), start = start), 
              data = pooled_average,
              se = FALSE,
              aes(color = factor(treatment))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

#--Make some Arrhenius plots
# First, calculate Boltzmann temperature
temp$invBT_eV <- 1/(0.00008617*(temp$temp_C + 273.15))
# Calculate mean 1/kT
library(plyr)
temp_means <- ddply(temp, .(treatment), summarize, invBT_eV_mean=mean(invBT_eV, na.rm=T))
# Add temp data to dataframe with fitted r and K
params$invBT_eV <- temp_means[match(params$temp_treatment,temp_means$treatment),"invBT_eV_mean"]

# Plot K vs. 1/kT
library(ggplot2)
library(scales)
ggplot(params, aes(x=invBT_eV, y=K, color=temp_treatment)) + 
  geom_smooth(method = "lm", se=FALSE, color="black") + 
  geom_point(size=3) + 
  xlab(expression(paste('Inverse temperature ', '<1/',italic('kT'),'>', ' (',  eV^{-1}, ')'))) +
  ylab(expression(paste('Carrying capacity ',italic('K'),' (ind)'))) +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
# Estimate E for K
K_Arrhenius <- lm(log(K) ~ invBT_eV, params)
summary(K_Arrhenius)
confint(K_Arrhenius)

# Plot r vs. 1/kT
ggplot(params, aes(x=invBT_eV, y=r0, color=temp_treatment)) + 
  geom_smooth(method = "lm", se=FALSE, color="black") + 
  geom_point(size=3) + 
  xlab(expression(paste('Inverse temperature ', '<1/',italic('kT'),'>', ' (',  eV^{-1}, ')'))) +
  ylab(expression(paste('Intrinsic rate of increase ',italic('r'),' (ind/ind/t)'))) +
  #ylab(expression("Intrinsic rate of increase, r")) +
  scale_y_continuous(trans="log", breaks = trans_breaks("log", function(x) exp(x), n=3), 
                     labels = trans_format("log", math_format(e^.x))) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
# Estimate E for K
r_Arrhenius <- lm(log(r0) ~ invBT_eV, params)
summary(r_Arrhenius)
confint(r_Arrhenius)


########################
# Fit each group's curve separately

### NOTE: having issues with convergence that I don't have time to solve now, see https://stackoverflow.com/questions/27547548/solving-error-message-step-halving-factor-reduced-below-minimum-in-nls-step-a

# Make loop to fit model and store parameters in a dataframe (https://stackoverflow.com/questions/71131348/extract-model-parameters-from-nls-model-for-loop)
library(tidyverse)
group <- unique(pooled_sum$section_group)

# define start to obtain the number and name of fit parameters
start <- list(N0 = 5, K = 60, r0 = 1.5)
# create empty data.frame to store IDs and parameters
params_group <- data.frame(matrix(nrow = length(group), ncol = 1+length(start)))
names(params_group) <- c("group",names(start))

#for loop to run through all the groups
# (changed to 'seq_along()' to get consecutive 'i's
#  for easier indexing of the 'param'-data.frame)
for (i in seq_along(group)){
  #Creating our dataframe for group "i"
  Individual_DFs2 <- pooled %>% filter (section_group %in% group[i]) 
  
  #Fit model for treatment "i"
  nls.floop <- nls(count ~ N.logistic(lab_week, N0, K, r0),
                   data = Individual_DFs2, 
                   start = start)
  
  # store IDs
  params_group[i,1] <- group[i]
  # store fit parameters
  params_group[i,2:ncol(params_group)] <- nls.floop$m$getPars()
}



