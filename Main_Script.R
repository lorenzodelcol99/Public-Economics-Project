# Uploading the Ferrier & Navvarro dataset
data <- read.csv("Dataset/DATA_MACRO_FN.csv")

GDP           <- ts(data$GDP, start = c(1913, 2), frequency = 4)
GOV           <- ts(data$GOV, start = c(1913, 2), frequency = 4)
NEWS          <- ts(data$NEWS, start = c(1913, 2), frequency = 4)
MTR           <- ts(data$MTR, start = c(1913, 2), frequency = 4)
ATR           <- ts(data$ATR, start = c(1913, 2), frequency = 4)


ts_data = ts(data, start = c(1913, 2), frequency = 4)

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
plot(NEWS)

# I order to build my Smooth Transition Local Projection estimator I need to build the
# "State" variable, the gamma (referring to the paper names)

gamma_raw = (MTR - ATR) / (1-ATR)
plot(gamma_raw)

# after having created my state variable I need to standardize it, in order to obtain Z scores
# will take values between -3 +3, this is good for our buil in logistic function in the 
# smooth transition LP
gamma = (gamma_raw - mean(gamma_raw, na.rm = T)) / sd(gamma_raw, na.rm = T)
plot(gamma)

# Now I Create the Logistic Transition Variable (F_z)
# We use gamma = 3 as the smoothness parameter (standard in literature)
#### PER ADESSO LASCIO COSì, IN SEGUITO CERCO DI USARE Granger and Terasvirta (1993) 
# LOGISTIC FUNCTION PER LO SMOOTH TRANSITION
gamma_param <- 3
F_z <- exp(-gamma_param * gamma) / (1 + exp(-gamma_param * gamma))

####### Creating the regressors for my specification

# What I need
# Delta Yt+h (dependent variable)
# Delta Gt+h (endogenous variable that need to be instumented)
# Instruments BP and RZ
# State Variable: Progressivity variable (gamma)
# Controls: lagged of growth rates of log GDP, lagged growth rates of gov, 
# lagged Tax Rates (MTR), and non linear and non linear time trends (t,t^2,t^3,t^4)

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
  BP,
  RZ = NEWS,
  
  # the state variable
  gamma_z = gamma, # standardized gamma
#  F_z = F_z,        # progressive regime weight (0,1) Logistic function
  
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


df_lp <- data.frame(LP_data)

library(lpirfs)

# 1. Prepare Data
df_lp <- as.data.frame(LP_data)

# 2. Controls List (Your manual lags + trends)
controls_list <- c("gdp_d_1", "gdp_d_2", "gdp_d_3", "gdp_d_4",
                   "gdp_d_5", "gdp_d_6", "gdp_d_7", "gdp_d_8",
                   "gov_d_1", "gov_d_2", "gov_d_3", "gov_d_4",
                   "gov_d_5", "gov_d_6", "gov_d_7", "gov_d_8",
                   "MTR_1",   "MTR_2",   "MTR_3",   "MTR_4",
                   "MTR_5",   "MTR_6",   "MTR_7",   "MTR_8",
                   "t", "t_2", "t_3", "t_4")

# 3. Run STLP-IV

results_stlp <- lp_nl_iv(
  
  # --- Variables ---
  endog_data      = df_lp[, "gdp", drop = FALSE],
  lags_endog_nl   = 0,                             # Keep this at 0
  
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
  hor             = 20.   # Number of h Horizon of the projecting
)

# Think Locally.. Project Globally!

# 4. Plot
plot(results_stlp)
