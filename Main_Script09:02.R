# Uploading the Ferrier & Navvarro dataset
data <- read.csv("Dataset/DATA_MACRO_FN.csv")

GDP           <- ts(data$GDP, start = c(1913, 2), frequency = 4)
GOV           <- ts(data$GOV, start = c(1913, 2), frequency = 4)
RZ            <- ts(data$NEWS, start = c(1913, 2), frequency = 4)
MTR           <- ts(data$MTR, start = c(1913, 2), frequency = 4)
ATR           <- ts(data$ATR, start = c(1913, 2), frequency = 4)

ts_data = ts(data, start = c(1913, 2), frequency = 4)


# In order to build my Smooth Transition Local Projection estimator I need to build the
# "State" variable, the gamma (referring to the Navarro's paper names)

gamma_raw = (MTR - ATR) / (1-ATR)
plot(gamma_raw)

# after having created my state variable I need to standardize it, in order to obtain Z scores
# will take values between -3 +3, this is good for our build in logistic function in the 
# smooth transition LP
gamma = (gamma_raw - mean(gamma_raw, na.rm = T)) / sd(gamma_raw, na.rm = T)
plot(gamma)

####### Creating the regressors for my specification

# What I need
# Delta Yt+h (dependent variable)
# Delta Gt+h (endogenous variable that need to be instumented)
# Instruments BP and RZ
# State Variable: Progressivity variable (gamma)
# Controls: lagged of growth rates of log GDP, lagged growth rates of gov, 
# lagged Tax Rates (MTR), and non linear and non linear time trends (t,t^2,t^3,t^4)

# In the first model specification I will use log GDP
# and log GOV, hence I will obtain elasticities of the multipliers.
# This means that I will need to re convert the multipliers in dollars
# I will need to scale them


# Building the controls 
# GDP
plot(GDP)
# Transforming GDP in logs
gdp <- log(GDP) 
plot(gdp)
# Taking the First difference of Log GDP
gdp_d <- diff(log(GDP))
plot(gdp_d)

## GOV

plot(GOV)
# transforming the GOV in logs
gov <- log(GOV)
plot(gov)

# Taking first differences of the log GOV
gov_d <- diff(log(GOV))
plot(gov_d)

library(dplyr)

# now I want to create 8 lags of the first difference of log GDP
gdp_d_1 <- stats::lag(gdp_d, -1)
gdp_d_2 <- stats::lag(gdp_d, -2)
gdp_d_3 <- stats::lag(gdp_d, -3)
gdp_d_4 <- stats::lag(gdp_d, -4)
gdp_d_5 <- stats::lag(gdp_d, -5)
gdp_d_6 <- stats::lag(gdp_d, -6)
gdp_d_7 <- stats::lag(gdp_d, -7)
gdp_d_8 <- stats::lag(gdp_d, -8)

# creating the lags for the first difference of the GOV
gov_d_1 <- stats::lag(gov_d, -1)
gov_d_2 <- stats::lag(gov_d, -2)
gov_d_3 <- stats::lag(gov_d, -3)
gov_d_4 <- stats::lag(gov_d, -4)
gov_d_5 <- stats::lag(gov_d, -5)
gov_d_6 <- stats::lag(gov_d, -6)
gov_d_7 <- stats::lag(gov_d, -7)
gov_d_8 <- stats::lag(gov_d, -8)

# creating the lags for the Tax Rate MTR
MTR_1 <- stats::lag(MTR, -1)
MTR_2 <- stats::lag(MTR, -2)
MTR_3 <- stats::lag(MTR, -3)
MTR_4 <- stats::lag(MTR, -4)
MTR_5 <- stats::lag(MTR, -5)
MTR_6 <- stats::lag(MTR, -6)
MTR_7 <- stats::lag(MTR, -7)
MTR_8 <- stats::lag(MTR, -8)


### Creating the time trends
T <- length(GDP)
time_index <- ts(1:T, start = c(1913, 1), frequency = 4)

t   <- time_index
t_2 <- time_index^2
t_3 <- time_index^3
t_4 <- time_index^4


#####################     Creating the BP Shock from a trivariate SVAR       #######################
library(vars)
Tax <- MTR * GDP
tax <- log(Tax)
tax_d <- diff((tax))
plot(tax_d)
var_data_ts <- ts.intersect(gov_d, gdp_d, tax_d)
# converting this into a data fram
var_data_df <- as.data.frame(var_data_ts)

# now let's run the Reduced From VAR with 4 lags

RF_VAR_BP <- VAR(var_data_df, p = 4, type = "const")
# Since Gov ('dg') is the FIRST variable, its orthogonalized shock is
# identical to its raw residual. We don't need SVAR() to prove this.
BP_Shock_Raw <- residuals(RF_VAR_BP)[, "gov_d"]

# Convert to Time Series with Correct Frequency
# The residuals start 4 quarters AFTER the data used in the VAR (due to 4 lags).
# We assume var_data_ts has the correct start/frequency

start_res <- tsp(var_data_ts)[1] + 4/4 # Shift start date forward by 1 year (4 quarters)
BP_Shock_TS <- ts(BP_Shock_Raw, start = start_res, frequency = 4)
plot(BP_Shock_Raw)


# now I need to "Orthogonalize" the BP shock, I need to "remove" the "Anticipated" portion
# of the shock leaving just the "Unanticipated" one
####### How? by regressing the the residuals of the SVAR on the NEWS (RZ) Shock
# collecting the residuals of this regression I'm extracting just the uncorrelated part of the BP
# shock from the RZ variable

aligned_data <- ts.intersect(BP_Shock_TS, RZ = ts_data[, "NEWS"])

BP_Orthogonalized_Model <- lm(BP_Shock_TS ~ RZ, data = aligned_data)


BP_Orthogonalides <- residuals(BP_Orthogonalized_Model)
BP <- ts(BP_Orthogonalides, start = start(aligned_data), frequency = frequency(aligned_data))
plot(BP)

#### Combining everything in a dataset with all the variable that I need for my STLP_IV
# adding the GDP and GOV in levels as they are needed to create the multipliers variables
LP_data <- ts.intersect(
  
  # In levels variables
  GDP,
  GOV,
  
  # In logs variabbles 
  gdp,
  gov,
  
  # the Instruments
  BP,
  RZ = RZ,
  
  # the state variable
  gamma_z = gamma, # standardized gamma
  
  # GDP growth rates lags
  gdp_d_1, gdp_d_2, gdp_d_3, gdp_d_4,
  gdp_d_5, gdp_d_6, gdp_d_7, gdp_d_8,
  
  # GOV growth rates lags
  gov_d_1, gov_d_2, gov_d_3, gov_d_4,
  gov_d_5, gov_d_6, gov_d_7, gov_d_8,
  
  # Tax rate lags
  MTR_1, MTR_2, MTR_3, MTR_4,
  MTR_5, MTR_6, MTR_7, MTR_8,
  
  # time trends
  t, t_2, t_3, t_4
)

##############################

# Now is the time to build the STLP-IV
# unfortunately the lpirfs package doesn't work with over identification, hence I will need to create
# the optimal instrument separatly
# Remember that I'm computing are elasticities, later I have to re convert 
# them in dollar terms


# Using the same dates that they use in the paper
LP_FN_dates <- window(LP_data, end = c(2006, 4))
# we need to transform the data from "time series object" to as data frame in
# order to run the LP
df_lp <- as.data.frame(LP_FN_dates)

# creating my list of CONTROLS (Lags + Trend)
controls_list <- c("gdp_d_1", "gdp_d_2", "gdp_d_3", "gdp_d_4",
                   "gdp_d_5", "gdp_d_6", "gdp_d_7", "gdp_d_8",
                   "gov_d_1", "gov_d_2", "gov_d_3", "gov_d_4",
                   "gov_d_5", "gov_d_6", "gov_d_7", "gov_d_8",
                   "MTR_1",   "MTR_2",   "MTR_3",   "MTR_4",
                   "MTR_5",   "MTR_6",   "MTR_7",   "MTR_8",
                   "t", "t_2", "t_3", "t_4")

################## SMOOTH TRANSITION INSTRUMENTAL VARIABLES LOCAL PROJECTIONS #######################

# FIRST SPECIFICATION, THE ONE THAT YIELDS THE RESULTS CONSISTENT WITH THE F&N PAPER
# just identified STLP-IV with the NEWS (RZ) Shock


############################ SPECIFICATION 1 : JUST RZ INSTRUMENT ##################################

library(lpirfs)

STLP_IV_RZ <- lp_nl_iv(
  
  # --- Variables ---
  endog_data      = df_lp[, "gdp", drop = FALSE],
  lags_endog_nl   = 0,                           
  
  shock           = df_lp[, "gov", drop = FALSE],
  instr           = df_lp[, "RZ",  drop = FALSE],
  
  # --- Controls ---
  contemp_data    = df_lp[, controls_list],
  
  # --- Smooth Transition ---
  switching       = df_lp[, "gamma_z", drop = FALSE],
  use_logistic    = TRUE,
  gamma           = 3,
  lag_switching   = FALSE,
  
  # --- Settings ---
  use_hp          = FALSE,
  lambda          = NaN,
  trend           = 0,
  cumul_mult      = TRUE,
  confint         = 1.96, # Confidence Interval
  hor             = 20   # Number of h Horizon of the projecting
)

# quick summary
summary(STLP_IV_RZ)

# Plotting the ELASTICITIES MULTIPLIERS --> NEED TO BE CONVERTED IN DOLLARS
plot(STLP_IV_RZ)
# good sign! both the IRFs are in line with the paper!!

# Now I want to obtain the Multipliers in dollars!
############################## CONVERTING THE MULTIPLIERS IN DOLLARS ##############################

# for doing so, I need the GDP/GOV ratio to convert Elasticity to Dollar Multiplier
scaling_factor_GDP_GOV <- mean(LP_FN_dates[, "GDP"] / LP_FN_dates[, "GOV"], na.rm = T)

############################### Plotting the Dollar Multipliers ###############################

# Now I need to collect the parameters of the IRF of the 
#                   Regime 1 (Low-Progressivity)

low_prog_df_RZ <- data.frame(
  Horizon = 0:(length(STLP_IV_RZ$irf_s1_mean) - 1),
  Mean    = as.numeric(STLP_IV_RZ$irf_s1_mean) * scaling_factor_GDP_GOV, # Mean --> the Beta coefficients
  Lower   = as.numeric(STLP_IV_RZ$irf_s1_low)  * scaling_factor_GDP_GOV, # Low  --> the lower bound of the 95% CI
  Upper   = as.numeric(STLP_IV_RZ$irf_s1_up)   * scaling_factor_GDP_GOV, # High --> the upper bound of the 95% CI
  Regime  = "Low Progressivity"
)

#                   Regime 2 (High-Progressivity)

high_prog_df_RZ <- data.frame(
  Horizon = 0:(length(STLP_IV_RZ$irf_s2_mean) - 1),
  Mean    = as.numeric(STLP_IV_RZ$irf_s2_mean) * scaling_factor_GDP_GOV,
  Lower   = as.numeric(STLP_IV_RZ$irf_s2_low)  * scaling_factor_GDP_GOV,
  Upper   = as.numeric(STLP_IV_RZ$irf_s2_up)   * scaling_factor_GDP_GOV,
  Regime  = "High Progressivity"
)

# collecting everything in one "data frame"
scaled_multipliers_2regimes_RZ <- bind_rows(low_prog_df_RZ, high_prog_df_RZ)
# Plotting the two regimes multipliers into a single one image

###### ONE SINGLE PLOT OF THE IRFs OF STLP-IV_RZ
library(ggplot2)

ggplot(scaled_multipliers_2regimes_RZ, aes(x = Horizon, y = Mean, color = Regime, fill = Regime)) + 
  # Plotting the Zero Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Plotting 95% Confidence Intervals
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) + 
  # Plotting the Betas (mean) Lines 
  geom_line(size = 1.2) + 
  # Customing the colors (Blue High prog, Red Low prog)
  scale_color_manual(values = c("High Progressivity" = "blue", "Low Progressivity" = "red")) + 
  scale_fill_manual (values = c("High Progressivity" = "royalblue2", "Low Progressivity" = "orangered")) +
  # Lables and Titles
  labs(
    title    = "Government Spending Multipliers in $",
    subtitle = "Smooth Transition LP-IV",
    y        = "Dollar Change in GDP / Dollar Change in Spending",
    x        = "Horizon, Quarters",
    color    = "Regime\n(Shaded: 95% Newey-West CI)",
    fill     = "Regime\n(Shaded: 95% Newey-West CI)",
    caption = paste0("Scaled by 1913-2006 avg GDP/GOV ratio")
  ) + 
  # Theme 
  theme_minimal() + 
  theme(
    legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )


#### High Progressivity Multiplier

ggplot(high_prog_df_RZ, aes(x = Horizon, y = Mean, color = Regime, fill = Regime)) + 
  # Plotting the Zero Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Plotting 95% Confidence Intervals
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) + 
  # Plotting the Betas (mean) Lines 
  geom_line(size = 1.2) + 
  # Customing the colors (Blue High prog, Red Low prog)
  scale_color_manual(values = c("High Progressivity" = "blue")) + 
  scale_fill_manual (values = c("High Progressivity" = "royalblue2")) +
  # Lables and Titles
  labs(
    title    = "Government Spending Multipliers in $: HIGH PROGRESSIVITY",
    subtitle = "Smooth Transition LP-IV",
    y        = "Dollar Change in GDP / Dollar Change in Spending",
    x        = "Horizon, Quarters",
    color    = "Regime\n(Shaded: 95% Newey-West CI)",
    fill     = "Regime\n(Shaded: 95% Newey-West CI)",
    caption = paste0("Scaled by 1913-2006 avg GDP/GOV ratio")
  ) + 
  # Theme 
  theme_minimal() + 
  theme(
    legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

#### Low Progressivity Multiplier

ggplot(low_prog_df_RZ, aes(x = Horizon, y = Mean, color = Regime, fill = Regime)) + 
  # Plotting the Zero Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Plotting 95% Confidence Intervals
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) + 
  # Plotting the Betas (mean) Lines 
  geom_line(size = 1.2) + 
  # Customing the colors (Blue High prog, Red Low prog)
  scale_color_manual(values = c("Low Progressivity" = "red")) + 
  scale_fill_manual (values = c("Low Progressivity" = "orangered")) +
  # Lables and Titles
  labs(
    title    = "Government Spending Multipliers in $: LOW PROGRESSIVITY",
    subtitle = "Smooth Transition LP-IV",
    y        = "Dollar Change in GDP / Dollar Change in Spending",
    x        = "Horizon, Quarters",
    color    = "Regime\n(Shaded: 95% Newey-West CI)",
    fill     = "Regime\n(Shaded: 95% Newey-West CI)",
    caption = paste0("Scaled by 1913-2006 avg GDP/GOV ratio")
  ) + 
  # Theme 
  theme_minimal() + 
  theme(
    legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

# Specification with NEWS shock is consistent with Ferriere&Navarro, different magnitute



## OTHER SPECIFICATIONS

############################### SPECIFICATION 2 : JUST BP (Orthogonalized) INSTRUMENT ################################

STLP_IV_BP <- lp_nl_iv(
  
  # --- Variables ---
  endog_data      = df_lp[, "gdp", drop = FALSE],
  lags_endog_nl   = 0,                           
  
  shock           = df_lp[, "gov", drop = FALSE],
  instr           = df_lp[, "BP",  drop = FALSE],
  
  # --- Controls ---
  contemp_data    = df_lp[, controls_list],
  
  # --- Smooth Transition ---
  switching       = df_lp[, "gamma_z", drop = FALSE],
  use_logistic    = TRUE,
  gamma           = 3,
  lag_switching   = FALSE,
  
  # --- Settings ---
  use_hp          = FALSE,
  lambda          = NaN,
  trend           = 0,
  cumul_mult      = TRUE,
  confint         = 1.96, # Confidence Interval
  hor             = 20   # Number of h Horizon of the projecting
)


# Now I need to collect the parameters of the IRF of the 
#                   Regime 1 (Low-Progressivity)

low_prog_df_BP <- data.frame(
  Horizon = 0:(length(STLP_IV_BP$irf_s1_mean) - 1),
  Mean    = as.numeric(STLP_IV_BP$irf_s1_mean) * scaling_factor_GDP_GOV, # Mean --> the Beta coefficients
  Lower   = as.numeric(STLP_IV_BP$irf_s1_low)  * scaling_factor_GDP_GOV, # Low  --> the lower bound of the 95% CI
  Upper   = as.numeric(STLP_IV_BP$irf_s1_up)   * scaling_factor_GDP_GOV, # High --> the upper bound of the 95% CI
  Regime  = "Low Progressivity"
)

#                   Regime 2 (High-Progressivity)

high_prog_df_BP <- data.frame(
  Horizon = 0:(length(STLP_IV_BP$irf_s2_mean) - 1),
  Mean    = as.numeric(STLP_IV_BP$irf_s2_mean) * scaling_factor_GDP_GOV,
  Lower   = as.numeric(STLP_IV_BP$irf_s2_low)  * scaling_factor_GDP_GOV,
  Upper   = as.numeric(STLP_IV_BP$irf_s2_up)   * scaling_factor_GDP_GOV,
  Regime  = "High Progressivity"
)

scaled_multipliers_2regimes_BP <- bind_rows(low_prog_df_BP, high_prog_df_BP)

# Plotting the two regimes multipliers into a single one image

ggplot(scaled_multipliers_2regimes_BP, aes(x = Horizon, y = Mean, color = Regime, fill = Regime)) + 
  # Plotting the Zero Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Plotting 95% Confidence Intervals
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) + 
  # Plotting the Betas (mean) Lines 
  geom_line(size = 1.2) + 
  # Customing the colors (Blue High prog, Red Low prog)
  scale_color_manual(values = c("High Progressivity" = "blue", "Low Progressivity" = "red")) + 
  scale_fill_manual (values = c("High Progressivity" = "royalblue2", "Low Progressivity" = "orangered")) +
  # Lables and Titles
  labs(
    title    = "Government Spending Multipliers in $",
    subtitle = "Smooth Transition LP-IV",
    y        = "Dollar Change in GDP / Dollar Change in Spending",
    x        = "Horizon, Quarters",
    color    = "Regime\n(Shaded: 95% Newey-West CI)",
    fill     = "Regime\n(Shaded: 95% Newey-West CI)",
    caption = paste0("Scaled by 1913-2006 avg GDP/GOV ratio")
  ) + 
  # Theme 
  theme_minimal() + 
  theme(
    legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

# this gives the exact OPPOSIT result wrt the other, and wrt the Ferriere & Navarro paper
# --> Interpretation: Low progressivity funded Goverment spending increases have positive and significant  multipliers
# and on the other hand the high progressivity ones are statistically equal to 0
# With the Surprise Blanchard&Perotti shock the Progressivity channel disappears, agents need lack of foresight



#### High Progressivity Multiplier 

ggplot(high_prog_df_BP, aes(x = Horizon, y = Mean, color = Regime, fill = Regime)) + 
  # Plotting the Zero Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Plotting 95% Confidence Intervals
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) + 
  # Plotting the Betas (mean) Lines 
  geom_line(size = 1.2) + 
  # Customing the colors (Blue High prog, Red Low prog)
  scale_color_manual(values = c("High Progressivity" = "blue")) + 
  scale_fill_manual (values = c("High Progressivity" = "royalblue2")) +
  # Lables and Titles
  labs(
    title    = "Government Spending Multipliers in $: HIGH PROGRESSIVITY",
    subtitle = "Smooth Transition LP-IV",
    y        = "Dollar Change in GDP / Dollar Change in Spending",
    x        = "Horizon, Quarters",
    color    = "Regime\n(Shaded: 95% Newey-West CI)",
    fill     = "Regime\n(Shaded: 95% Newey-West CI)",
    caption = paste0("Scaled by 1913-2006 avg GDP/GOV ratio")
  ) + 
  # Theme 
  theme_minimal() + 
  theme(
    legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

#### Low Progressivity Multiplier

ggplot(low_prog_df_BP, aes(x = Horizon, y = Mean, color = Regime, fill = Regime)) + 
  # Plotting the Zero Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Plotting 95% Confidence Intervals
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) + 
  # Plotting the Betas (mean) Lines 
  geom_line(size = 1.2) + 
  # Customing the colors (Blue High prog, Red Low prog)
  scale_color_manual(values = c("Low Progressivity" = "red")) + 
  scale_fill_manual (values = c("Low Progressivity" = "orangered")) +
  # Lables and Titles
  labs(
    title    = "Government Spending Multipliers in $: LOW PROGRESSIVITY",
    subtitle = "Smooth Transition LP-IV",
    y        = "Dollar Change in GDP / Dollar Change in Spending",
    x        = "Horizon, Quarters",
    color    = "Regime\n(Shaded: 95% Newey-West CI)",
    fill     = "Regime\n(Shaded: 95% Newey-West CI)",
    caption = paste0("Scaled by 1913-2006 avg GDP/GOV ratio")
  ) + 
  # Theme 
  theme_minimal() + 
  theme(
    legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

# this gives the exact OPPOSIT result wrt the other, and wrt the Ferriere & Navarro paper
# --> Interpretation: Low progressivity funded Goverment spending increases have positive and significant  multipliers
# and on the other hand the high progressivity ones are statistically equal to 0

# Why I have obtain this result?
# --> Econometric Puzzle!
# Lack of foresight of households





############# SPECIFICATION 3 : Over Identification RZ + BP (clean) INSTRUMENTS ###################
# If I have two instruments, why don't I use it both!!
# as Professor Chiara Dal Bianco tought me 
# I need to create the "optimal Instrument" that combine both my instruments

# gov ~ RZ + BP + [All Controls]
full_formula <- as.formula(paste("gov ~ RZ + BP +", paste(controls_list, collapse = " + ")))

# Run the First Stage Regression
first_stage_model <- lm(full_formula, data = df_lp)

# Store the 'Optimal Instrument' (Fitted Values)
df_lp$Optimal_Instrument <- fitted(first_stage_model)

STLP_IV_OVERID <- lp_nl_iv(
  
  # --- Variables ---
  endog_data      = df_lp[, "gdp", drop = FALSE],
  lags_endog_nl   = 0,                           
  
  shock           = df_lp[, "gov", drop = FALSE],
  instr           = df_lp[, "Optimal_Instrument",  drop = FALSE],
  
  # --- Controls ---
  contemp_data    = df_lp[, controls_list],
  
  # --- Smooth Transition ---
  switching       = df_lp[, "gamma_z", drop = FALSE],
  use_logistic    = TRUE,
  gamma           = 3,
  lag_switching   = FALSE,
  
  # --- Settings ---
  use_hp          = FALSE,
  lambda          = NaN,
  trend           = 0,
  cumul_mult      = TRUE,
  confint         = 1.96, # Confidence Interval
  hor             = 20   # Number of h Horizon of the projecting
)

# Summary of the STLP-IV
summary(STLP_IV_OVERID)

plot(STLP_IV_OVERID)
# the two impulse responses look equal.. let's see what happens with the dollar multiplier

#                   Regime 1 (Low-Progressivity)

low_prog_df_OVERID <- data.frame(
  Horizon = 0:(length(STLP_IV_OVERID$irf_s1_mean) - 1),
  Mean    = as.numeric(STLP_IV_OVERID$irf_s1_mean) * scaling_factor_GDP_GOV, # Mean --> the Beta coefficients
  Lower   = as.numeric(STLP_IV_OVERID$irf_s1_low)  * scaling_factor_GDP_GOV, # Low  --> the lower bound of the 95% CI
  Upper   = as.numeric(STLP_IV_OVERID$irf_s1_up)   * scaling_factor_GDP_GOV, # High --> the upper bound of the 95% CI
  Regime  = "Low Progressivity"
)

#                   Regime 2 (High-Progressivity)

high_prog_df_OVERID <- data.frame(
  Horizon = 0:(length(STLP_IV_OVERID$irf_s2_mean) - 1),
  Mean    = as.numeric(STLP_IV_OVERID$irf_s2_mean) * scaling_factor_GDP_GOV,
  Lower   = as.numeric(STLP_IV_OVERID$irf_s2_low)  * scaling_factor_GDP_GOV,
  Upper   = as.numeric(STLP_IV_OVERID$irf_s2_up)   * scaling_factor_GDP_GOV,
  Regime  = "High Progressivity"
)

scaled_multipliers_2regimes_OVERID <- bind_rows(low_prog_df_OVERID, high_prog_df_OVERID)

# Plotting the two regimes multipliers into a single one image

ggplot(scaled_multipliers_2regimes_OVERID, aes(x = Horizon, y = Mean, color = Regime, fill = Regime)) + 
  # Plotting the Zero Line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Plotting 95% Confidence Intervals
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) + 
  # Plotting the Betas (mean) Lines 
  geom_line(size = 1.2) + 
  # Customing the colors (Blue High prog, Red Low prog)
  scale_color_manual(values = c("High Progressivity" = "blue", "Low Progressivity" = "red")) + 
  scale_fill_manual (values = c("High Progressivity" = "royalblue2", "Low Progressivity" = "orangered")) +
  # Lables and Titles
  labs(
    title    = "Government Spending Multipliers in $",
    subtitle = "Smooth Transition LP-IV",
    y        = "Dollar Change in GDP / Dollar Change in Spending",
    x        = "Horizon, Quarters",
    color    = "Regime\n(Shaded: 95% Newey-West CI)",
    fill     = "Regime\n(Shaded: 95% Newey-West CI)",
    caption = paste0("Scaled by 1913-2006 avg GDP/GOV ratio")
  ) + 
  # Theme 
  theme_minimal() + 
  theme(
    legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )


# I have a specification that yield their results, but is a different specification wrt their, 
# and my other specifications are complitely different, the just identified by BP specification 
# yields the right opposit result and the over identified one, the specification with the optimal
# instrument gives that there is no heterogeneity and at the 95% level there is no multiplier

# Loss of identification of the progressivity channel


##################################################################################################

########################## Think Globally.. Project Locally !!!   #################################







