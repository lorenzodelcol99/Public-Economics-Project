# Uploading the Ferrier & Navvarro dataset
data <- read.csv("Dataset/DATA_MACRO_FN.csv")

GDP           <- ts(data$GDP, start = c(1913, 2), frequency = 4)
GOV           <- ts(data$GOV, start = c(1913, 2), frequency = 4)
RZ            <- ts(data$NEWS, start = c(1913, 2), frequency = 4)
MTR           <- ts(data$MTR, start = c(1913, 2), frequency = 4)
ATR           <- ts(data$ATR, start = c(1913, 2), frequency = 4)


ts_data = ts(data, start = c(1913, 2), frequency = 4)

### Old BP shock, the "Naive One"

# creating the BP Shock

library(dplyr)
data_reg <- data %>%
  mutate(
    # Logs
    ln_GOV = log(GOV),
    ln_GDP = log(GDP),
    
    # Growth Rates (Differences)
    d_gov = c(NA, diff(ln_GOV)),
    d_gdp = c(NA, diff(ln_GDP)),
    
    # Create Lag Columns explicitly (lags 1 to 4)
    # This makes the regression line much cleaner later
    gov_L1 = dplyr::lag(d_gov, 1),
    gov_L2 = dplyr::lag(d_gov, 2),
    gov_L3 = dplyr::lag(d_gov, 3),
    gov_L4 = dplyr::lag(d_gov, 4),
    
    gdp_L1 = dplyr::lag(d_gdp, 1),
    gdp_L2 = dplyr::lag(d_gdp, 2),
    gdp_L3 = dplyr::lag(d_gdp, 3),
    gdp_L4 = dplyr::lag(d_gdp, 4)
  )

BP_model <- lm(d_gov ~ gov_L1 + gov_L2 + gov_L3 + gov_L4 + 
                 gdp_L1 + gdp_L2 + gdp_L3 + gdp_L4, 
               data = data_reg,
               na.action = na.exclude)

# storing the Shock
BP <- residuals(BP_model)
BP           <- ts(BP, start = c(1913, 2), frequency = 4)

plot(BP)

# A shorter way to do the same stuff
library(dynlm)
# If I use the library(dynlm) package, I can directly estimate this Dinamic Linear Regression
# in a single row
# --> d() = first difference --> L( d() , 1:4) from t to t-4 using the Lag operator

#
bp_model_raw <- dynlm(d(log(GOV)) ~ L(d(log(GOV)), 1:4) + L(d(log(GDP)), 1:4), 
                      data = ts_data)

# Store the raw residuals
BP_Raw <- residuals(bp_model_raw)

plot(BP_Raw)

############################# Creating the Clean BP Shock ################################
library(dynlm)

# I order to "clean" the BP shock, I need to "regress" the BP_raw shock over the "NEWS" RZ shock
# and take the residuals --> which are going to be "all the variability of the BP shock that
# do not correlates with the NEWS shock
# those residuals are going to be my Clean_BP shock

# ok before we need to 

# 1. Create a temporary dataset aligning BP_Raw and RZ
# ts.intersect drops non-overlapping dates, ensuring 1914:Q3 matches 1914:Q3
cleaning_data <- ts.intersect(BP_Raw, RZ = ts_data[,"NEWS"])

# 2. Run the Cleaning Regression on this ALIGNED data
# We regress BP on Current News + 4 Lags of News
# Note: We use 'data = cleaning_data' to ensure it uses the aligned series
BP_Clean_Model <- dynlm(BP_Raw ~ RZ + L(RZ, 1:4), data = cleaning_data)

# 3. Extract the Clean Shock (Residuals)
BP_Clean <- residuals(BP_Clean_Model)

# 4. Visual Check
par(mfrow=c(1,1))
plot.ts(BP_Raw, col = "red", lwd = 2, ylab = "Shock", 
        main = "Orthogonalization Effect: Raw (Red) vs Clean (Blue)")
lines(BP_Clean, col = "blue", lwd = 1)
#legend("topright", legend=c("Raw BP", "Clean BP"), col=c("red", "blue"), lty=1)

# -------------------------------------------------------------------------
# PREPARING LP_DATA (Updated)
# -------------------------------------------------------------------------
# Now we use the aligned BP_Clean in the final intersection
LP_data <- ts.intersect(
  GDP, GOV, gdp, gov,
  BP_Clean,          # <--- The correctly aligned clean shock
  RZ = RZ,
  gamma_z = gamma,
  gdp_d_1, gdp_d_2, gdp_d_3, gdp_d_4, gdp_d_5, gdp_d_6, gdp_d_7, gdp_d_8,
  gov_d_1, gov_d_2, gov_d_3, gov_d_4, gov_d_5, gov_d_6, gov_d_7, gov_d_8,
  MTR_1, MTR_2, MTR_3, MTR_4, MTR_5, MTR_6, MTR_7, MTR_8,
  t, t_2, t_3, t_4
)

plot(BP_Clean)

par(mfrow=c(2,1))
plot.ts(BP_Raw, main="Raw BP Shock (Raw)", col="red", ylab="")
plot.ts(BP_Clean, main="Orthogonalized BP Shock (Clean)", col="blue", ylab="")


par(mfrow=c(1,1))
# 2. Plot the "Raw" (Dirty) Shock first in RED
# We set a generous y-axis limit (ylim) to make sure both lines fit
plot.ts(BP_Raw, 
        col = "red", 
        lwd = 3.5,
        ylab = "% of GDP / Shock Size",
        main = "The Effect of Orthogonalization: Raw vs. Clean BP Shock")

# 3. Add the "Clean" (Orthogonalized) Shock in BLUE
# lines() automatically matches the dates if they are ts objects
lines(BP_Clean, col = "blue", lwd = 1.5)

# 4. Add a grid and legend for clarity
grid()
legend("topright", 
       legend = c("Raw BP (Contains Old News)", "Clean BP (Pure Surprise)"),
       col = c("red", "blue"), 
       lty = 1, 
       lwd = 1.5,
       bg = "white")


### Ok now I want to see if the if the STLP-IV give different results with the Clean BP

# Hence I plug here the previous part of the code
# --> 

# I order to build my Smooth Transition Local Projection estimator I need to build the
# "State" variable, the gamma (referring to the Navarro's paper names)

gamma_raw = (MTR - ATR) / (1-ATR)
plot(gamma_raw)

# after having created my state variable I need to standardize it, in order to obtain Z scores
# will take values between -3 +3, this is good for our build in logistic function in the 
# smooth transition LP
gamma = (gamma_raw - mean(gamma_raw, na.rm = T)) / sd(gamma_raw, na.rm = T)
plot(gamma)

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
  BP_Clean,
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

# Now is the tiem to build the STLP-IV
# unfortunatly the lpirfs package doesn't work with over identification
# problem, this that I'm computing are elasticities, later I have to re convert 
# them in dollar terms

library(dynlm)
library(lpirfs)


LP_FN_dates <- window(LP_data, end = c(2006, 4))

library(lpirfs)

# 1. Prepare Data
df_lp <- as.data.frame(LP_FN_dates)

# 2. Controls List (Your manual lags + trends)
controls_list <- c("gdp_d_1", "gdp_d_2", "gdp_d_3", "gdp_d_4",
                   "gdp_d_5", "gdp_d_6", "gdp_d_7", "gdp_d_8",
                   "gov_d_1", "gov_d_2", "gov_d_3", "gov_d_4",
                   "gov_d_5", "gov_d_6", "gov_d_7", "gov_d_8",
                   "MTR_1",   "MTR_2",   "MTR_3",   "MTR_4",
                   "MTR_5",   "MTR_6",   "MTR_7",   "MTR_8",
                   "t", "t_2", "t_3", "t_4")

# 3. Run STLP-IV with the CLEAN BP Instrument

# 1) STLP-IV with BP Shock


STLP_IV_BP <- lp_nl_iv(
  
  # --- Variables ---
  endog_data      = df_lp[, "gdp", drop = FALSE],
  lags_endog_nl   = 0,                           
  
  shock           = df_lp[, "gov", drop = FALSE],
  instr           = df_lp[, "BP_Clean",  drop = FALSE],
  
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

scaling_factor_GDP_GOV <- mean(LP_FN_dates[, "GDP"] / LP_FN_dates[, "GOV"], na.rm = T)

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


# is not working. let's see if with the over identified setup it works

#################################### 2) OVER IDENFITICATION STLP-IV 

# I need to create the "optimal Instrument" that combine both my instruments

# gov ~ RZ + BP + [All Controls]
full_formula <- as.formula(paste("gov ~ RZ + BP_Clean +", paste(controls_list, collapse = " + ")))

# Run the First Stage Regression
first_stage_model <- lm(full_formula, data = df_lp)

# Store the 'Optimal Instrument' (Fitted Values)
df_lp$Optimal_Instrument <- fitted(first_stage_model)

####### DA FINIRE DI SMANETTARE


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

# Also here it doesn't work
