# Uploading the Ferrier & Navvarro dataset
data <- read.csv("Dataset/DATA_MACRO_FN.csv")

GDP           <- ts(data$GDP, start = c(1913, 2), frequency = 4)
GOV           <- ts(data$GOV, start = c(1913, 2), frequency = 4)
NEWS          <- ts(data$NEWS, start = c(1913, 2), frequency = 4)
UNEMP         <- ts(data$UNEMP, start = c(1913, 2), frequency = 4)
TB3           <- ts(data$TB3, start = c(1913, 2), frequency = 4)
RDEF          <- ts(data$RDEF, start = c(1913, 2), frequency = 4)
MTR           <- ts(data$MTR, start = c(1913, 2), frequency = 4)
ATR           <- ts(data$ATR, start = c(1913, 2), frequency = 4)
MITR          <- ts(data$MITR, start = c(1913, 2), frequency = 4)
AITR          <- ts(data$AITR, start = c(1913, 2), frequency = 4)
ATR_PSZ       <- ts(data$ATR_PSZ, start = c(1913, 2), frequency = 4)
ATR_B90_PSZ   <- ts(data$ATR_B90_PSZ, start = c(1913, 2), frequency = 4)
ATR_B50_PSZ   <- ts(data$ATR_B50_PSZ, start = c(1913, 2), frequency = 4)
ATR_T10_PSZ   <- ts(data$ATR_T10_PSZ, start = c(1913, 2), frequency = 4)
ATR_T05_PSZ   <- ts(data$ATR_T05_PSZ, start = c(1913, 2), frequency = 4)
ATR_T01_PSZ   <- ts(data$ATR_T01_PSZ, start = c(1913, 2), frequency = 4)
psoc          <- ts(data$psoc, start = c(1913, 2), frequency = 4)
pfed          <- ts(data$pfed, start = c(1913, 2), frequency = 4)
rinvfx        <- ts(data$rinvfx, start = c(1913, 2), frequency = 4)
rnri          <- ts(data$rnri, start = c(1913, 2), frequency = 4)
wgenofarm     <- ts(data$wgenofarm, start = c(1913, 2), frequency = 4)
hours         <- ts(data$hours, start = c(1913, 2), frequency = 4)

# RZ Shock (Defence News Shock)
plot(NEWS)

# I want to create the Blanchard Perotti Shock
# BP Shock
# plot(BP)

plot(GDP)
plot(GOV)
plot(MTR)
plot(ATR)

# Transforming the data
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



plot(data_reg$GOV)
plot(data_reg$ln_GOV)
plot(data_reg$d_gov)

plot(data_reg$GDP)
plot(data_reg$ln_GDP)
plot(data_reg$d_gdp)


BP_model <- lm(d_gov ~ gov_L1 + gov_L2 + gov_L3 + gov_L4 + 
                 gdp_L1 + gdp_L2 + gdp_L3 + gdp_L4, 
               data = data_reg,
               na.action = na.exclude)

# storing the Shock
BP <- residuals(BP_model)
BP           <- ts(BP, start = c(1913, 2), frequency = 4)

plot(BP, 
     main = "Blanchard-Perotti Spending Shock", 
     ylab = "Shock", 
     xlab = "Year",
     col = "blue",
     lwd = 1.5)
abline(h = 0, col = "red", lty = 2)


plot(BP, 
     main = "BP ans RZ Shocks",
     ylab = "Shock",
     xlab = "Year",
     col = "blue",
     lwd = 1.5)
lines(NEWS,
      col = "red")
abline(h = 0, col = "black", lty = 2)
legend("topright", legend = c("BP", "RZ"), 
       col=c("blue", "red"), lty = 1, cex=0.7)

####################################################################################

# I order to build my Smooth Transaction Local Projection estimator I need to build the
# "State" variable, the gamma (referring to the paper names)

gamma_raw = (MTR - ATR) / (1-ATR)
plot(gamma_raw)

# after having created my state variable I need to standardize it, in order to obtain Z scores
# will take values between -3 +3, this is good for our buil in logistic function in the 
# smooth transaction LP
gamma = (gamma_raw - mean(gamma_raw, na.rm = T)) / sd(gamma_raw, na.rm = T)
plot(gamma)


############







