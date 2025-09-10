###################################################################################################
# Code for the Financial Volatility project
# Authors: Jana Cvetinovska and Jacopo Rauseo
# HSG

####################################################################################################

# Packages
library(fBasics) 
library(zoo)
library(fGarch) 
library(tseries)
library(car)
library(systemfit)
library(mvtnorm)
library(quadprog)
library(VGAM)
library(sandwich)
library(readxl)
library(ggplot2)
library(rugarch)
library(strucchange)
library(MSGARCH)
library(lmtest)
library(MCS)


####################################################################################################
# setwd("~/Master_HSG/Spring 2025/Financial Volatility/Project")

# Import the data
Oil_daily <- read_excel("Data.xlsx", 
                        sheet = 1, skip = 1)


###################################################################################################
################################# PLOT of the TS ##################################################
###################################################################################################


# Recession period according to NBER
recession_start <- as.Date("2020-02-01")
recession_end   <- as.Date("2020-04-30")

# Oil PLOT
ggplot(Oil_daily, aes(x = as.Date(Date), y = Last_Price)) +
  geom_rect(aes(xmin = recession_start, xmax = recession_end,
                ymin = -50, ymax = Inf),
            fill = "grey", alpha = 0.4) +
  geom_line() +
  
  # Add arrow and label for Ukraine War
  annotate("segment", x = as.Date("2021-05-01"), xend = as.Date("2022-02-01"),
           y = 95, yend = 88,
           arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  annotate("text", x = as.Date("2019-08-01"), y = 100,
           label = "Russia-Ukraine War", color = "red", hjust = 0) +
  
  # Add annotation for COVID Lockdown
  annotate("segment", x = as.Date("2018-10-01"), xend = as.Date("2020-02-01"),
           y = -45, yend = 10,
           arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
  annotate("text", x = as.Date("2018-01-02"), y = -50,
           label = "COVID-19 recession", color = "blue", hjust = 0) +
  
  # Add annotation for COVID Vaccine
  annotate("segment", x = as.Date("2022-03-01"), xend = as.Date("2021-03-30"),
           y = 24, yend = 50,
           arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
  annotate("text", x = as.Date("2021-02-01"), y = 20,
           label = "COVID Vaccine announced as safe", color = "blue", hjust = 0) +
  
  labs(title = "Daily Oil Price",
       x = "Date", y = "Price (USD/bbl)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



# Daily Oil Price
Oil_daily$Date <- as.Date(Oil_daily$Date)
op_day <- zoo(Oil_daily$Last_Price, order.by = Oil_daily$Date)


# Deal with negative price of oil
#transform series to non-negative
op_day <- na.omit(op_day)
op_day[op_day < 0] <- 0.000000000000001



###################################################################################################
#################################### Simple Returns #############################################
###################################################################################################

or_day <- diff(op_day)/op_day[1:(length(op_day)-1)]

# Take the returns: use head to get all the data except the last obs
or_day <- diff(op_day)/head(op_day, -1)


# Plots of the returns
plot(or_day, main = "Oil: Daily Returns", xlab= "Time", ylab="Return")

# Statistics for the simple returns
basicStats(or_day)



###################################################################################################
####################################### Log returns ###################################################
###################################################################################################

# Log returns
log_ord <- log(1+or_day) #oil

log_ord = na.omit(log_ord)


# Plots of the returns
plot(log_ord, main = "Oil: Daily Log Returns", xlab= "Time", ylab="Return")

# Statistics for the 
cbind(basicStats(log_ord), basicStats(log_grd))
jarque.bera.test(coredata(log_ord))

# ACF plot
acf(coredata(log_ord), main = "ACF Daily Oil Log Returns")

# PACF plot
pacf(coredata(log_ord), main = "PACF Daily Oil Log Returns")



###################################################################################################
####################################### FUNCTIONS ###################################################
###################################################################################################

# asymptotic acf function
gamma=function(x,h)
{
  n=length(x)
  h=abs(h)
  x=x-mean(x)
  gamma=sum(x[1:(n-h)]*x[(h+1):n])/n
}


rho=function(x,h)
{
  rho=gamma(x,h)/gamma(x,0)
}

#Asymptotical cov. matrix of the sample autocorrelations
gamma.asy=function(x,h) 
{
  n=length(x)
  h=abs(h)
  x=x-mean(x)
  x2=x^2
  gamma.asy<-matrix(,h,h)  
  for (i in 1:h)
  {
    for (j in i:h)
    {
      gamma.asy[i,j]=gamma.asy[j,i]=sum(x[(j-i+1):(n-i)]*x[1:(n-j)]*x2[(j+1):n])/n
    }
  }
  rho.asy=1/gamma(x,0)^2*gamma.asy
  list(gamma.asy=gamma.asy,rho.asy=rho.asy)
}

# ARCH test
LM=function(x,h)
{
  n=length(x)
  x2=x^2-mean(x^2)
  dat=matrix(,n-h,h+1)
  for (i in 1:(h+1))
  {
    dat[,i]=x2[(h+2-i):(n-i+1)]
  }
  a=lm(dat[,1]~dat[,2:(h+1)])
  r2=summary(a)$r.squared
  print(r2 * n)
  print(1-pchisq(r2*n,h))
}

#Corrected portmanteau test under GARCH assumption
corr.Box.test=function(x,h) 
{
  n<-length(x)
  a=gamma.asy(x,h)
  acf.val=sapply(c(1:h),function(h) rho(x,h))
  val=n*(acf.val%*%solve(a$rho.asy)%*%acf.val)
  print(paste("Corrected statistic: " , val))  
  print(paste("p-value: " ,1-pchisq(val,h)))
  print(paste("lag: " , h))
}


# plot with corrected bound for ARCH effect
n1.acf=function(x,main=NULL,method="NP")
{
  n=length(x)
  nlag=as.integer(min(10*log10(n),n-1))
  acf.val=sapply(c(1:nlag),function(h) rho(x,h))
  x2=x^2
  var= 1+(sapply(c(1:nlag),function(h) gamma(x2,h)))/gamma(x,0)^2
  band=sqrt(var/n)
  minval=1.2*min(acf.val,-1.96*band,-1.96/sqrt(n))
  maxval=1.2*max(acf.val,1.96*band,1.96/sqrt(n))
  acf(x,xlab="Lag",ylab="Sample autocorrelations",ylim=c(minval,maxval),main=main)
  lines(c(1:nlag),-1.96*band,lty=1,col="red")
  lines(c(1:nlag),1.96*band,lty=1,col="red")
}

#Estimate an asymmetric GARCH(1,1) model with Student's t innovations
my.loglike.t=function(theta) 
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
  
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start
  
  my.sigma=c(sqrt(my.sigmasq[1]),rep(0,n))
  
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }
  
  for (i in 2:(n+1))
  {
    my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)
    
    #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4](data[i-1]-my.mean[i-1])^2((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    #my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])((data[i-1]-my.mean[i-1])&gt;0)-theta[4]^2(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }
  
  #my.sigmasq=my.sigma^2
  #my.sigmasq=exp(log.sigmasq)
  
  1/2*sum(log(my.sigmasq[2:(n+1)](theta[6]-2)/theta[6])) - sum(log(dt((data[2:(n+1)]-my.mean[2:(n+1)])/sqrt(my.sigmasq[2:(n+1)](theta[6]-2)/theta[6]),df=theta[6])))+10^(10)(theta[6]<2)+10^(10)(theta[6]>10)
}


sigmasq.model=function(theta)
{
  n=length(returns)
  x.start= mean(returns)
  sigmasq.start= var(returns)
  
  data=c(x.start,returns)
  my.sigmasq= rep(0,n+1)
  my.sigmasq[1]=sigmasq.start
  
  my.sigma=sqrt(my.sigmasq)
  
  my.mean=rep(0,n+1)
  for(j in 2:(n+1))
  {
    my.mean[j]=theta[1] #Constant conditional mean
  }
  
  
  for (i in 2:(n+1))
  {
    my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4]^2*my.sigmasq[i-1] #GARCH(1,1)
    
    #my.sigmasq[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])^2 + theta[4](data[i-1]-my.mean[i-1])^2((data[i-1]-my.mean[i-1])<=0)+theta[5]^2*my.sigmasq[i-1] #GJR-GARCH(1,1)
    
    #my.sigma[i]=theta[2]^2 + theta[3]^2*(data[i-1]-my.mean[i-1])((data[i-1]-my.mean[i-1])&gt;0)-theta[4]^2(data[i-1]-my.mean[i-1])*((data[i-1]-my.mean[i-1])<=0)+ theta[5]^2*my.sigma[i-1] #TGARCH(1,1)
    
    #my.sigma[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))+theta[5]^2*my.sigma[i-1] #PGARCH(1,1) with d=1
    #my.sigmasq[i]=theta[2]^2+theta[3]^2*(abs((data[i-1]-my.mean[i-1]))-theta[4]*(data[i-1]-my.mean[i-1]))^2+theta[5]^2*my.sigmasq[i-1] #PGARCH(1,1) with d=2
  }
  
  my.sigmasq=my.sigma^2
  #my.sigmasq=exp(log.sigmasq)
  
  
  list(my.sigmasq = my.sigmasq[2:(n + 1)],my.mean=my.mean[2:(n+1)])
}


fit_custom_garch <- function(returns) {
  assign("returns", returns, envir = .GlobalEnv)  # required for your current function definition
  par.start <- c(rep(0.5, 5), 4)
  
  opt_result <- tryCatch(
    nlm(my.loglike.t, par.start, iterlim = 1000, print.level = 0),
    error = function(e) return(NULL)
  )
  
  if (is.null(opt_result)) return(NULL)
  
  loglik <- -opt_result$minimum
  k <- 6
  aic <- 2 * k - 2 * loglik
  
  return(list(
    params = opt_result$estimate,
    loglik = loglik,
    aic = aic
  ))
}



# FUNCTION FOR THE PREDICTIONS FOR GARCH, EGARCH, GJRGARCH

predictions <- function(m, garch,a,b) {
  ########
  # Inputs
  #### m = length of the window
  #### garch = name of the garch for the unrgach specifications
  #### a,b = orders of the garch
  # Outputs
  #### prediction plot 
  #### dataset with predictions
  #### specifiactions of the model
  ########
  
  Time <- length(log_ord)
  N <- Time - m + 1  # Number of rolling windows
  
  spec <- ugarchspec(
    variance.model = list(model = garch, garchOrder = c(a,b)),
    mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
    distribution.model = "std"
  )
  
  # df for storing the coefficients
  garch_coeff <- data.frame(
    alpha1 = rep(NA, N),
    beta1 = rep(NA, N),
    one_step_pred = rep(NA, N)
  )
  
  for (i in 1:N) {
    wind <- log_ord[i:(i + m - 1)]  
    fit <- ugarchfit(spec, data = wind, solver.control = list(trace = 0))  
    par <- coef(fit)
    
    garch_coeff[i, "alpha1"] <- par["alpha1"]
    garch_coeff[i, "beta1"] <- par["beta1"]
    
    pred <- ugarchforecast(fit, n.ahead = 1)
    garch_coeff[i, "one_step_pred"] <- sigma(pred)  
  }
  
  
  x <- time(log_ord)
  pred_garch <- zoo(garch_coeff[,"one_step_pred"],x[m:Time-1])
  
  
  plot(pred_garch, ,xlab = "Time", ylab = "Conditional Volatility", 
       main = paste("One-Day ahead forecast ", garch))
  
  return(list(
    coeff_df = garch_coeff,
    forecast = pred_garch,
    spec = spec
  ))
}


# Hansen's SPA test for superior predictive ability
spa <- function(per, bench = 1, m = 9, obs = 926, q = 0.25, iter = 1, periodogram = TRUE, k = seed) {
  # per: matrix of losses (rows = observations, cols = models)
  # bench: index of benchmark model
  # m: number of models
  # obs: number of observations
  # q: parameter for bootstrap chain (0 < q < 1)
  # iter: number of bootstrap iterations
  # periodogram: use periodogram (TRUE) or lag window (FALSE) for variance estimation
  set.seed(k)
  e <- bench  # benchmark index
  d <- matrix(NA, nrow = obs, ncol = m - 1)  # loss differentials matrix
  s <- 0
  for (i in seq_len(m)) {
    if (i != e) {
      s <- s + 1
      d[, s] <- per[, e] - per[, i]
    }
  }
  
  w <- numeric(m - 1)  # variance estimators
  
  for (k in 1:(m - 1)) {
    # Covariance function or spectrum of the differential loss
    acf_cov <- acf(d[, k], lag.max = obs - 1, type = "covariance", plot = FALSE)$acf
    
    if (!periodogram) {
      # Weighted sum of autocovariances (lag window estimator)
      lags <- 1:(obs - 1)
      weights <- ((obs - lags)/obs) * (1 - q)^lags + (lags/obs) * (1 - q)^(obs - lags)
      w[k] <- sqrt(acf_cov[1] + 2 * sum(weights * acf_cov[-1]))
    } else {
      # Periodogram variance estimator
      w[k] <- sqrt(spectrum(d[, k], plot = FALSE)$spec[1])
    }
  }
  
  test_stats <- sqrt(obs) * colMeans(d) / w
  stat <- max(0, max(test_stats))
  
  # Bootstrap
  stat_boot <- numeric(iter)
  
  for (r in 1:iter) {
    tau <- integer(obs)
    tau[1] <- sample(obs, 1)
    
    for (i in 2:obs) {
      s_r <- runif(1)
      tau[i] <- if (s_r < q) sample(obs, 1) else min(tau[i - 1] + 1, obs)
    }
    
    d_boot <- d[tau, , drop = FALSE]
    
    e_boot <- matrix(NA, nrow = obs, ncol = m - 1)
    for (k in 1:(m - 1)) {
      mean_dk <- mean(d[, k])
      threshold <- -sqrt(w[k]^2 / obs * 2 * log(log(obs)))
      adjustment <- ifelse(mean_dk >= threshold, mean_dk, 0)
      e_boot[, k] <- d_boot[, k] - adjustment
    }
    
    stat_boot[r] <- max(0, max(sqrt(obs) * colMeans(e_boot) / w))
  }
  
  p_value <- mean(stat_boot > stat)
  
  list(p.value = p_value, stat.boot = stat_boot, stat = stat)
}

# ORDER SELECTION FOR GARCH
order_selection_garch <- function(df){
  # Fit different oderers
  garch_1_1 = garchFit(~garch(1,1),data=df, include.mean = TRUE, cond.dist = "std", trace=F)
  G1_1 <- garch_1_1@fit$ics["AIC"]
  
  garch_1_2 = garchFit(~garch(1,2),data=df,, include.mean = TRUE, cond.dist = "std", trace=F)
  G1_2 <- garch_1_2@fit$ics["AIC"]
  
  garch_2_1 = garchFit(~garch(2,1),data=df,, include.mean = TRUE, cond.dist = "std", trace=F)
  G2_1 <-garch_2_1@fit$ics["AIC"]
  
  garch_2_2 = garchFit(~garch(2,1),data=df,, include.mean = TRUE, cond.dist = "std", trace=F)
  G2_2 <-garch_2_2@fit$ics["AIC"]
  
  garch_1_0 = garchFit(~garch(1,0),data=df,, include.mean = TRUE, cond.dist = "std", trace=F)
  G1_0 <-garch_1_0@fit$ics["AIC"]
  
  garch_2_0 = garchFit(~garch(2,0),data=df,, include.mean = TRUE, cond.dist = "std", trace=F)
  G2_0 <-garch_1_0@fit$ics["AIC"]
  
  garch_2_3 = garchFit(~garch(2,3),data=df,, include.mean = TRUE, cond.dist = "std", trace=F)
  G2_3 <-garch_2_3@fit$ics["AIC"]
  
  garch_3_2 = garchFit(~garch(3,2),data=df,, include.mean = TRUE, cond.dist = "std", trace=F)
  G3_2 <-garch_3_2@fit$ics["AIC"]
  
  garch_3_3 = garchFit(~garch(3,3),data=df,, include.mean = TRUE, cond.dist = "std", trace=F)
  G3_3 <-garch_3_3@fit$ics["AIC"]
  
  
  return(rbind(G1_0,G2_0,G1_1,G1_2, G2_1, G2_2, G2_3, G3_2, G3_3))
}


# Student t distribution E[z^4]
E4 <- function(nu) {
  if (nu <= 4) {
    return(NA)  # fourth moment does not exist
  }
  return(3 * nu^2 / ((nu - 2) * (nu - 4)))
}



# Lag function for zoo
lag_zoo <- function(x, k = 1) {
  lag(x, -k, TRUE)
}

# LEVERAGE TEST FUNCTION
lm_test <- function(x){
  # Create required variables
  eps_squared <- coredata(x)^2
  eps_lag <- lag_zoo(coredata(x), 1)
  
  # Construct the variables
  S_minus <- as.numeric(coredata(x) < 0)
  eps_minus <- S_minus * coredata(x)
  eps_plus  <- as.numeric(coredata(x) > 0) * coredata(x)
  
  # Align and drop NA due to lag
  df <- na.omit(data.frame(
    eps2 = eps_squared,
    S_minus = lag(S_minus, 1),
    eps_minus = lag(eps_minus, 1),
    eps_plus = lag(eps_plus, 1)
  ))
  
  # Fit the regressions
  reg_sign_bias <- lm(eps2 ~ S_minus, data = df)
  reg_neg_size_bias <- lm(eps2 ~ eps_minus, data = df)
  reg_pos_size_bias <- lm(eps2 ~ eps_plus, data = df)
  
  # Output the summaries
  cat("Sign Bias Test:\n")
  print(summary(reg_sign_bias))
  
  cat("\nNegative Size Bias Test:\n")
  print(summary(reg_neg_size_bias))
  
  cat("\nPositive Size Bias Test:\n")
  print(summary(reg_pos_size_bias))
}







###################################################################################################
###################################################################################################
###################################################################################################
# 4. ESTIMATION 
###################################################################################################
###################################################################################################
###################################################################################################


##################################################################################################
# 3. Structural breaks 
###################################################################################################




# No breaks in the mean
ts <- breakpoints(coredata(log_ord) ~ 1)
summary(ts)


# Define the square series for the 
squares = log_ord^2

# Perform a breakpoints test
bp <- breakpoints(coredata(squares) ~ 1)

# Summary of breakpoints
summary(bp)

# Plot the breakpoints
plot(bp)
lines(squares)
bp$breakpoints
index(squares)[bp$breakpoints]

cut1 <- as.Date("2020-02-21")
cut2 <- as.Date("2021-03-30")


# Windows definition
log_ord_1 <- window(log_ord, end = cut1)
log_ord_2 <- window(log_ord, start = cut1 + 1, end = cut2)
log_ord_3 <- window(log_ord, start = cut2 + 1)


###################################################################################################
##################################################################################################
# 4.0. CONDITIONAL MEAN
###################################################################################################
###################################################################################################


##################################################################################################
# Jarque-Bera for the windows
#  Corrected Box Ljung
###################################################################################################

# Testing normality
print("Jarque-Brera Test Window 1")
jarque.bera.test(coredata(log_ord_1))
print("Jarque-Brera Test Window 2")
jarque.bera.test(coredata(log_ord_1))
print("Jarque-Brera Test Window 3")
jarque.bera.test(coredata(log_ord_1))


# Corrected Box Ljung
print("Corrected Box.Ljung Window 1")
corr.Box.test(coredata(log_ord_1),5)
corr.Box.test(coredata(log_ord_1),10)

print("Corrected Box.Ljung Window 2")
corr.Box.test(coredata(log_ord_2),5)
corr.Box.test(coredata(log_ord_2),10)

print("Corrected Box.Ljung Window 3")
corr.Box.test(coredata(log_ord_3),5)
corr.Box.test(coredata(log_ord_3),10)


# ARCH effect test
LM(log_ord_1,10)
LM(log_ord_2,10)
LM(log_ord_3,10)


# ACF plot
acf(coredata(log_ord_1), main = "ACF Daily Oil Log Returns Window 1")
acf(coredata(log_ord_2), main = "ACF Daily Oil Log Returns Window 2")
acf(coredata(log_ord_3), main = "ACF Daily Oil Log Returns Window 3")

# PACF plot
pacf(coredata(log_ord_1), main = "PACF Daily Oil Log Returns Window 1")
pacf(coredata(log_ord_2), main = "PACF Daily Oil Log Returns Window 2")
pacf(coredata(log_ord_3), main = "PACF Daily Oil Log Returns Window 3")



# Corrected ACF/PACF plots 
n1.acf(coredata(log_ord),main=c("Daily Log Oily Prices Corrected Bounds")) 
n1.acf(coredata(log_ord_1),main=c("Daily Log Oily Prices Corrected Bounds Window 1")) 
n1.acf(coredata(log_ord_2),main=c("Daily Log Oily Prices Corrected Bounds Window 2")) 
n1.acf(coredata(log_ord_3),main=c("Daily Log Oily Prices Corrected Bounds Window 3")) 




###################################################################################################
##################################################################################################
# 4.0. SQUARES: CONDITIONAL VARIANCE
###################################################################################################
###################################################################################################


# Square ACF plot
squares_1 = log_ord_1^2
squares_2 = log_ord_2^2
squares_3 = log_ord_3^2

# ACF plot
acf(coredata(squares_1), main = "ACF Squared Daily Oil Log Returns Window 1")
acf(coredata(squares_2), main = "ACF Squared Daily Oil Log Returns Window 2")
acf(coredata(squares_3), main = "ACF Squared Daily Oil Log Returns Window 3")

# PACF plot
pacf(coredata(log_ord_1), main = "PACF Squared Daily Oil Log Returns Window 1")
pacf(coredata(log_ord_2), main = "PACF Squared Daily Oil Log Returns Window 2")
pacf(coredata(log_ord_3), main = "PACF Squared Daily Oil Log Returns Window 3")




# reject long memory test: also visual no long memory
a=c()
for (i in 1:length(log_ord))
{
  a=c(a,sum(log_ord[1:i]^2-mean(log_ord^2)))
}

1/sqrt(length(log_ord))*1/sqrt(var(log_ord^2))*(max(a)-min(a)) #R/S statistic with sample variance


1/sqrt(length(log_ord))*1/sqrt(lrvar(log_ord^2,
                                     type="Newey-West")*length(log_ord^2))*(max(a)-min(a)) #R/S statistic with HAC variance



#################################################################################################
#################################################################################################
#################################################################################################
# WINDOWS: FIT MODELS
#################################################################################################
#################################################################################################
#################################################################################################


#################################################################################################
# GARCH MODEL ESTIMATION
#################################################################################################


# Order selection
order_selection_garch(log_ord_1) #(1,1)
order_selection_garch(log_ord_2) #(1,2)
order_selection_garch(log_ord_3) #(1,1)



# Window 1: GARCH(1,1)
garch1 = garchFit(~garch(1,1),data=log_ord_1, include.mean = FALSE, cond.dist = "std", 
                  trace=F)
summary(garch1)
par(mfrow=c(2,2))
a=time(log_ord_1)
plot(zoo(garch1@residuals/garch1@sigma.t,a),xlab="",ylab="",main="GARCH(1,1) residuals")
qqnorm(garch1@residuals/garch1@sigma.t)
acf(garch1@residuals/garch1@sigma.t,main="GARCH(1,1) residuals")
acf(garch1@residuals^2/garch1@sigma.t^2,main="GARCH(1,1) squared residuals")


garch2 = garchFit(~garch(1,2),data=log_ord_2, include.mean = TRUE, cond.dist = "std", 
                  trace=F)
summary(garch2)
par(mfrow = c(2, 2))

a <- time(log_ord_2)
plot(zoo(garch2@residuals / garch2@sigma.t, a), xlab = "", ylab = "", main = "GARCH(1,2) residuals")
std_res <- residuals(garch2, standardize = TRUE)
nu <- garch2@fit$par["shape"]
qqplot(qt(ppoints(length(std_res)), df = nu), std_res,
       main = "Q-Q Plot: Residuals vs Student-t",
       xlab = "Theoretical Quantiles (Student-t)",
       ylab = "Standardized Residuals")
abline(0, 1, col = "red")
acf(garch2@residuals / garch2@sigma.t, main = "ACF: Standardized Residuals")
acf((garch2@residuals / garch2@sigma.t)^2, main = "ACF: Squared Std Residuals")


# Window 3: GARCH(1,1)
garch3 = garchFit(~garch(1,1),data=log_ord_3, include.mean = TRUE, cond.dist = "std", 
                  trace=F)
summary(garch3)
par(mfrow=c(2,2))

a <- time(log_ord_3)
plot(zoo(garch3@residuals / garch3@sigma.t, a), xlab = "", ylab = "", main = "GARCH(1,1) residuals")
std_res <- residuals(garch3, standardize = TRUE)
nu <- garch3@fit$par["shape"]
qqplot(qt(ppoints(length(std_res)), df = nu), std_res,
       main = "Q-Q Plot: Residuals vs Student-t",
       xlab = "Theoretical Quantiles (Student-t)",
       ylab = "Standardized Residuals")
abline(0, 1, col = "red")
acf(garch3@residuals / garch3@sigma.t, main = "ACF: Standardized Residuals")
acf((garch3@residuals / garch3@sigma.t)^2, main = "ACF: Squared Std Residuals")





# Plot the combined volatility with windows breaks
par(mfrow=c(1,1))
garch_vol <- zoo(
  c(garch1@sigma.t, garch2@sigma.t, garch3@sigma.t),
  c(time(log_ord_1), time(log_ord_2), time(log_ord_3))
)


plot(garch_vol, xlab = "Time", ylab = "Conditional Volatility", 
     main = "Conditional Volatility from the fitted GARCH")
abline(v = c(as.Date(cut1), as.Date(cut2)), col = "red", lty = 2, lwd = 2)




# Coefficients and moments condition
tab <-data.frame(coef(garch1)[3:4],coef(garch2)[3:4],coef(garch3)[3:4])
tab["beta2",] <-c(0,coef(garch2)[5],0)
tab["m1",] <- c(sum(tab[,1]),sum(tab[,2]),sum(tab[,3]))
tab["m2",] <- c(tab[4,1]^2 + (E4(garch1@fit$par["shape"])-1)*tab[1,1],
                tab[4,2]^2 + (E4(garch2@fit$par["shape"])-1)*tab[1,2],
                tab[4,2]^2 + (E4(garch3@fit$par["shape"])-1)*tab[1,3])

print(tab)




# Predictions
# garch1_pred <- predict(garch1,5) #Predictions until month t+5
# garch1_pred




#################################################################################################
# EGARCH MODEL ESTIMATION
#################################################################################################


# positive returns
a=c()
for (i in 1:length(log_ord))
{
  a=c(a,max(log_ord[i],0))
}
a

# correlation 
for (h in 1:40)
{
  print(h)
  print(cor(a[(1+h):(length(a))],log_ord[(1):(length(log_ord)-h)])) #Leverage effect
}




# Test for leverage
lm_test(log_ord)
lm_test(log_ord_1) #Negative coefficient: leverage 
lm_test(log_ord_2) #Negative coefficient in the sign bias test: leverage
lm_test(log_ord_3) #Same



# Set the specification for the egarch
spec <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "std")

spec2 <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 2)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "std")

# Windows
egarch1 <- ugarchfit(spec = spec, data = log_ord_1)
show(egarch1) 
coef(egarch1)  

egarch2 <- ugarchfit(spec = spec2, data = log_ord_2)
show(egarch2) 
coef(egarch2)

egarch3 <- ugarchfit(spec = spec, data = log_ord_3)
show(egarch3) 
coef(egarch3)



# Volatility EGARCH plot
par(mfrow=c(1,1))
egarch_vol <- zoo(
  c(sigma(egarch1), sigma(egarch2), sigma(egarch3)),
  c(time(log_ord_1), time(log_ord_2), time(log_ord_3))
)

plot(egarch_vol,xlab = "Time", ylab = "Conditional Volatility", 
     main = "Estimated volatility (E-GARCH)")
abline(v = c(as.Date(cut1), as.Date(cut2)), col = "red", lty = 2, lwd = 2)




# Controlling the residuals 
residuals_check <-  function(garch, series){
  res <- residuals(garch, standardize = TRUE)
  par(mfrow=c(2,2))
  a=time(log_ord_3)
  plot(zoo(res),xlab="",ylab="",main="EGARCH(1,1) residuals")
  qqnorm(res)
  acf(res,main="EGARCH(1,1) residuals")
  acf(res^2,main="EGARCH(1,1) squared residuals")
}


# Residuals check
residuals_check(egarch1, log_ord_1)
residuals_check(egarch2, log_ord_2)
residuals_check(egarch3, log_ord_3)


#News impact curve 
nic1 <- newsimpact(egarch1)
nic2 <- newsimpact(egarch2)
nic3 <- newsimpact(egarch3)


df_nic <- data.frame(
  shock = nic1$zx,
  variance1 = nic1$zy,
  variance2 = nic2$zy,
  variance3 = nic3$zy
)

df_nic

par(mfrow=c(1,1))
plot(df_nic$shock, df_nic$variance1, type = "l", col = "blue",
     ylab = "Conditional Variance", xlab = "Shock",
     main = "News Impact Curves")
lines(df_nic$shock, df_nic$variance3, col = "darkgreen")
legend("topright", legend = c("Window 1", "Window 3"),
       col = c("blue", "darkgreen"), lty = 1)



#################################################################################################
# GJRGARCH MODEL ESTIMATION
#################################################################################################



# Set the specification for the gjrgarch
t <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "std")

t2 <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 2)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "std")

# Windows
tgarch1 <- ugarchfit(spec = t, data = log_ord_1)
show(tgarch1) 
coef(tgarch1)  

tgarch2 <- ugarchfit(spec = t, data = log_ord_2)
show(tgarch2) 
coef(tgarch2)

tgarch3 <- ugarchfit(spec = t, data = log_ord_3)
show(tgarch3) 
coef(tgarch3)


# Volatility plot
par(mfrow=c(1,1))
tgarch_vol <- zoo(
  c(sigma(tgarch1), sigma(tgarch2), sigma(tgarch3)),
  c(time(log_ord_1), time(log_ord_2), time(log_ord_3))
)

plot(tgarch_vol,xlab = "Time", ylab = "Conditional Volatility", 
     main = "Estimated Volatility (GJR-GARCH)")
abline(v = c(as.Date(cut1), as.Date(cut2)), col = "red", lty = 2, lwd = 2)


residuals_check <-  function(garch, series){
  res <- residuals(garch, standardize = TRUE)
  par(mfrow=c(2,2))
  a=time(log_ord_3)
  plot(zoo(res),xlab="",ylab="",main="GJRGARCH(1,1) residuals")
  qqnorm(res)
  acf(res,main="GJRGARCH(1,1) residuals")
  acf(res^2,main="GJRGARCH(1,1) squared residuals")
}


# Residuals check
residuals_check(tgarch1, log_ord_1)
residuals_check(tgarch2, log_ord_2)
residuals_check(tgarch3, log_ord_3)



# News Impact curve
nic1 <- newsimpact(tgarch1)
nic2 <- newsimpact(tgarch2)
nic3 <- newsimpact(tgarch3)

df_nic <- data.frame(
  shock = nic1$zx,
  variance1 = nic1$zy,
  variance2 = nic2$zy,
  variance3 = nic3$zy
)

par(mfrow=c(1,1))
plot(df_nic$shock, df_nic$variance1, type = "l", col = "blue",
     ylab = "Conditional Variance", xlab = "Past Innovations",
     main = "News Impact Curves")
lines(df_nic$shock, df_nic$variance2, col = "red")
lines(df_nic$shock, df_nic$variance3, col = "darkgreen")
legend("topright", legend = c("Model 1", "Model 2", "Model 3"),
       col = c("blue", "red", "darkgreen"), lty = 1)










#################################################################################################
#################################################################################################
#################################################################################################
# MARKOV REGIME SWITCHING MODEL
#################################################################################################
#################################################################################################
#################################################################################################



train <- window(log_ord, start = as.Date("2018-01-01"), end = as.Date("2025-03-31"))
test <- window(log_ord, start = as.Date("2025-04-01"), end = as.Date("2025-05-07"))


# Compare 2 vs 3 regimes
spec2 <- CreateSpec(variance.spec = list(model = c("gjrGARCH", "gjrGARCH")), 
                    distribution.spec = list(distribution = rep("std", 2)),
                    switch.spec = list(do.mix = FALSE))
GRSgarch2 <- FitML(spec2, data = coredata(train))

spec3 <- CreateSpec(variance.spec = list(model = c("gjrGARCH", "gjrGARCH", "gjrGARCH")), 
                    distribution.spec = list(distribution = rep("std", 3)),
                    switch.spec = list(do.mix = FALSE))
GRSgarch3 <- FitML(spec3, data = coredata(train))


summary(GRSgarch2)
summary(GRSgarch3)
# Only two regimes



# Extract smoothed probabilities and convert to matrix
probs <- State(GRSgarch2, type = "smoothing")


probs_mat <- do.call(cbind, probs)  # Convert list to matrix

probs_mat

# Plot smoothed probabilities
matplot(probs_mat, type = "l", lty = 1, col = c("red", "blue"),
        ylab = "Probability", main = "Smoothed Regime Probabilities")

# Add legend
legend("topright", legend = paste("Regime", 1:2),
       col = c("red", "blue"), lty = 1)



vol <- Volatility(GRSgarch2)
plot(x = time(train),y = vol, type = "l", col = "black", lwd = 2,
     ylab = "Conditional Volatility", xlab = "Time",
     main = "Estimated Volatility (MSGARCH)")


pred <- predict(GRSgarch2, nahead = 25, do.return.draw = TRUE)
pred_mrsgarch <- as.numeric(pred$vol) 





#################################################################################################
#################################################################################################
#################################################################################################
# PREDICTIONS:
#################################################################################################
#################################################################################################
#################################################################################################

# MODEL: GARCH(1,1)
# ROLLING WINDOW ANALYSIS
# Window 15.04.2025 - 07.05.2025


# Parameters for the rolling window
m <- 1000  # Window size
Time <- length(log_ord) # length of the TS
N <- Time - m + 1  # Number of rolling windows


spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)


# Dataframe for storing the coefficients
garch11_coeff <- data.frame(
  alpha1 = rep(NA, N),
  beta1 = rep(NA, N),
  se_a = rep(NA, N),
  se_b = rep(NA, N),
  upper_bound_a = rep(NA, N),
  lower_bound_a = rep(NA, N),
  upper_bound_b = rep(NA, N),
  lower_bound_b = rep(NA, N),
  one_step_pred = rep(NA, N)
)


for (i in 1:N) {
  wind <- log_ord[i:(i + m - 1)]  
  fit <- ugarchfit(spec, data = wind, solver.control = list(trace = 0))  
  par <- coef(fit)
  
  garch11_coeff[i, "alpha1"] <- par["alpha1"]
  garch11_coeff[i, "beta1"] <- par["beta1"]
  
  se_mat <- robust_se <- sqrt(diag(fit@fit$robust.cvar))
  garch11_coeff[i, "se_a"] <- se_mat[2]
  garch11_coeff[i, "se_b"]  <- se_mat[3]
  
  garch11_coeff[i,"upper_bound_a"] <- garch11_coeff[i,"alpha1"] + 2*garch11_coeff[i, "se_a"]
  garch11_coeff[i,"lower_bound_a"] <- garch11_coeff[i,"alpha1"] - 2*garch11_coeff[i, "se_a"]
  
  garch11_coeff[i,"upper_bound_b"]  <- garch11_coeff[i,"beta1"] + 2*garch11_coeff[i, "se_b"]
  garch11_coeff[i,"lower_bound_b"] <- garch11_coeff[i,"beta1"] - 2*garch11_coeff[i, "se_b"]
  
  pred <- ugarchforecast(fit, n.ahead = 1)
  garch11_coeff[i, "one_step_pred"] <- sigma(pred)  
}


# Plot of the predictions
x <- time(log_ord)
pred_garch_11 <- zoo(garch11_coeff[,"one_step_pred"],x[1000:Time-1])
plot(pred_garch_11, ,xlab = "Time", ylab = "Conditional Volatility", 
     main = "One-Day ahead forecast GARCH(1,1)")



# Plots of the coefficients 
ggplot(data = garch11_coeff[1:600,], aes(x = 1:600, y = alpha1)) +
  geom_line() +
  geom_line(aes(y = lower_bound_a), colour = "blue") +
  geom_line(aes(y = upper_bound_a), colour = "blue") +
  labs(title = "GARCH(1,1): Alpha coefficient and confidence bounds",
       x = "Window", y = "Values") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(data = garch11_coeff[1:600,], aes(x = 1:600, y = beta1)) +
  geom_line() +
  geom_line(aes(y = lower_bound_b), colour = "blue") +
  geom_line(aes(y = upper_bound_b), colour = "blue") +
  labs(title = "GARCH(1,1): Beta coefficient and confidence bounds",
       x = "Window", y = "Values") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



# MODEL: EGARCH(1,1) + GJRGARCH(1,1)
# ROLLING WINDOW ANALYSIS

egarch_pred <- predictions(m = 1000, garch = "eGARCH", 1,1)
gjrgarch_pred <- predictions(m = 1000, garch = "gjrGARCH", 1,1)

data_egarch <- egarch_pred$coeff_df
data_gjrgarch <- gjrgarch_pred$coeff_df


date <- time(log_ord)


# Dataframe with all the predictions and the proxy
evaluation_df <- data.frame("Date" = date[m:Time-1], "Squared_r" = log_ord[m:Time-1]^2, 
                            "GARCH" = garch11_coeff$one_step_pred^2,
                            "EGARCH" = data_egarch$one_step_pred^2, 
                            "GJRGARCH" = data_gjrgarch$one_step_pred^2)



# evaluation_df <- read.csv("predictions.csv")

# Redefine the length 
evaluation_df <- evaluation_df[(851-24):851,]



#################################################################################################
#################################################################################################
# 5. DIRECT COMPARISON
#################################################################################################
#################################################################################################


# Add the MRSGARCH and tree
evaluation_df[,"MRSGARCH"] <- pred_mrsgarch^2
# tree_df <- read.csv("predictions_tree.csv")
tree_df <- read.csv("predictions_tree.csv")
evaluation_df[,"TREE"] <- tree_df[1803:1827,"volatility"]^2


# write.csv(evaluation_df, file = "predictions3.csv", row.names = FALSE)

evaluation_df <- read.csv("predictions.csv")

# Days of predictions
predictions_range <- evaluation_df$Date

pred_res <- data.frame("GARCH" = c(0,0,0), "EGARCH" = c(0,0,0), "GJRGARCH" = c(0,0,0), 
                       "MRSGARCH" = c(0,0,0), "TREE" = c(0,0,0),
                       row.names = c("MSE","MAE","QLIKE"))



pred_res[1,"GARCH"] = mean((evaluation_df[,"Squared_r"]-evaluation_df[,"GARCH"])^2)
pred_res[1,"EGARCH"] = mean((evaluation_df[,"Squared_r"]-evaluation_df[,"EGARCH"])^2)
pred_res[1,"GJRGARCH"] = mean((evaluation_df[,"Squared_r"]-evaluation_df[,"GJRGARCH"])^2)
pred_res[1,"MRSGARCH"] = mean((evaluation_df[,"Squared_r"]-evaluation_df[,"MRSGARCH"])^2)
pred_res[1,"TREE"] = mean((evaluation_df[,"Squared_r"]-evaluation_df[,"TREE"])^2)


pred_res[2,"GARCH"] = 1/length(evaluation_df[,"Squared_r"])*sum(abs(evaluation_df[,"Squared_r"]-evaluation_df[,"GARCH"]))
pred_res[2,"EGARCH"] = mean(abs(evaluation_df[,"Squared_r"]-evaluation_df[,"EGARCH"]))
pred_res[2,"GJRGARCH"] = mean(abs(evaluation_df[,"Squared_r"]-evaluation_df[,"GJRGARCH"]))
pred_res[2,"MRSGARCH"] = mean(abs(evaluation_df[,"Squared_r"]-evaluation_df[,"MRSGARCH"]))
pred_res[2,"TREE"] = mean(abs(evaluation_df[,"Squared_r"]-evaluation_df[,"TREE"]))


pred_res[3,"GARCH"] = mean(log(evaluation_df[,"GARCH"]) + evaluation_df[,"Squared_r"]/evaluation_df[,"GARCH"])
pred_res[3,"EGARCH"] = mean(log(evaluation_df[,"EGARCH"]) + evaluation_df[,"Squared_r"]/evaluation_df[,"EGARCH"])
pred_res[3,"GJRGARCH"] = mean(log(evaluation_df[,"GJRGARCH"]) + evaluation_df[,"Squared_r"]/evaluation_df[,"GJRGARCH"])
pred_res[3,"MRSGARCH"] = mean(log(evaluation_df[,"MRSGARCH"]) + evaluation_df[,"Squared_r"]/evaluation_df[,"MRSGARCH"])
pred_res[3,"TREE"] = mean(log(evaluation_df[,"TREE"]) + evaluation_df[,"Squared_r"]/evaluation_df[,"TREE"])


print(pred_res)



#################################################################################################
# 5. SPA
#################################################################################################



# 25 x 5 matrix of predictions
pred_matrix <- as.matrix(evaluation_df[,3:7])  

# Loss for each step-ahead prediction
loss_mse <- (pred_matrix - evaluation_df$Squared_r)^2
loss_mae <- abs(pred_matrix - evaluation_df$Squared_r)

loss_qlike <- matrix(NA, nrow = nrow(pred_matrix), ncol = ncol(pred_matrix))

rownames(pred_matrix) <- NULL
colnames(pred_matrix) <- NULL
rownames(evaluation_df) <- NULL

for (j in 1:5){
  for (i in 1:25){
    loss_qlike[i,j] = log(pred_matrix[i,j]) + evaluation_df[i,"Squared_r"]/pred_matrix[i,j]
  }
}




# Store the results
p_val = data.frame("GARCH" = c(0,0,0), "EGARCH" = c(0,0,0), "GJRGARCH" = c(0,0,0), 
                               "MRSGARCH" = c(0,0,0), "TREE" = c(0,0,0),
                               row.names = c("MSE","MAE","QLIKE"))


for (i in 1:5) {
  result <- spa(per = loss_mse, bench = i, m = 5, obs = 25, q = 0.25, iter = 1000, 
                periodogram = TRUE, k = 122333)
  p_val[1,i] <- result$p.value
}

#MAE
for (i in 1:5) {
  result <- spa(per = loss_mae, bench = 1, m = 5, obs = 25, q = 0.25, iter = 1000, 
                periodogram = TRUE, k = 122333)
  p_val[2,i] <- result$p.value
}

#QLIKE
for (i in 1:5) {
  result <- spa(per = loss_qlike, bench = 1, m = 5, obs = 25, q = 0.25, iter = 1000, 
                periodogram = TRUE, k = 122333)
  p_val[3,i] <- result$p.value
}


# results
print(p_val)


#################################################################################################
# 5. SCA
#################################################################################################



# Run MCS
mcs_result <- MCSprocedure(
  Loss = loss_qlike,
  alpha = 0.05,      
  B = 5000,          
  statistic = "Tmax" 
)
summary(mcs_result)





# Predictions plot
ggplot(data = evaluation_df, aes(x = as.Date(Date))) +
  geom_line(aes(y = GARCH, color = "GARCH"), size = 1) +
  geom_line(aes(y = EGARCH, color = "EGARCH"), size = 1) +
  geom_line(aes(y = GJRGARCH, color = "GJRGARCH"), size = 1) +
  geom_line(aes(y = MRSGARCH, color = "MRSGARCH"), size = 1) +
  geom_line(aes(y = TREE, color = "TREE"), size = 1) +
  geom_line(aes(y = Squared_r, color = "Squared Returns"), size = 1) +
  scale_color_manual(
    name = "Model",
    values = c(
      "GARCH" = "black",
      "EGARCH" = "blue",
      "GJRGARCH" = "red",
      "MRSGARCH" = "green",
      "TREE" = "pink",
      "Squared Returns" = "purple"
    )
  ) +
  labs(title = "GARCH-type Model Predictions",
       x = "Date", y = "Predicted Variance") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )




#################################################################################################
#################################################################################################
#################################################################################################
# TREE
#################################################################################################
#################################################################################################
#################################################################################################




#Function to fit GARCH(1,1) inside the tree
fit_garch_aic <- function(returns) {
  na.omit(returns)
  fit <- garchFit(~ garch(1, 1), 
                  data = returns, 
                  include.mean = TRUE,
                  cond.dist = "std",
                  trace = FALSE)
  aic  <- as.numeric(fit@fit$ics["AIC"])
  
  return(list(model = fit, aic = aic))
}


#Function for tree GARCH
#Please make sure data is a dataframe for this
#returns_col is the name of the returns column as a string
#covariates is a list of covariate names as strings
tree_garch <- function(data , covariates, returns_col, min_node_size = 100){
  #Stop if not enough data in node for reliable GARCH
  if (nrow(data) < min_node_size){
    model <- fit_garch_aic(data[[returns_col]])
    return(list(
      type = "leaf",
      model = model,
      aic = model$aic,
      data = data
    ))
  }
  
  best_score = Inf #Create Variable that stores best score (Akaike Criterion)
  best_split = NULL #Create variable that stores best split
  best_sets = NULL #Create variable that stores best sets
  
  for(var in covariates){
    #Splits to try, here you randomly select 50 quantiles of each covariate to try
    thresholds <- unique(quantile(data[[var]], probs = seq(0, 1, length.out = 20), na.rm = TRUE))
    
    for (t in thresholds) {
      #For each treshold, take the values less or equal as the left subtree, 
      #the values greater or equal as right subtree
      left <- data[data[[var]] <= t, ]
      right <- data[data[[var]] > t, ]
      #If there are too few observations in either the left or right subtree skip
      #that iteration of the loop and continue to the next
      if (nrow(left) < min_node_size || nrow(right) < min_node_size) next
      #Store Akaike criteria for each subtree
      aic_left <- fit_garch_aic(left[[returns_col]])$aic
      aic_right <- fit_garch_aic(right[[returns_col]])$aic
      total_aic <- aic_left + aic_right
      print(str(total_aic))
      #If the akaike is an improvement over the previous best, update the best
      #Also store the best split (covariate and treshold) and best subtree
      if (total_aic < best_score) {
        best_score <- total_aic
        best_split <- list(var = var, threshold = t)
        best_sets <- list(left = left, right = right)
      }
    }
  }
  #If no optimal split is found, stop the recursion
  if (is.null(best_split)) {
    model <- fit_garch_aic(data[[returns_col]])
    return(list(
      type = "leaf",
      model = model,
      aic = model$aic,
      data = data
    ))
  }
  #Before splitting the current node, fit a GARCH model on 
  #unsplit data to be used for pruning
  parent_model <- fit_garch_aic(data[[returns_col]])
  #Continue the algorithm onto the new subtrees
  left_subtree <- tree_garch(best_sets$left, covariates, returns_col, min_node_size)
  right_subtree <- tree_garch(best_sets$right, covariates, returns_col, min_node_size)
  #Output of tree function, gives the split by covariates and treshold, the left and right subtrees, 
  #the akaike criterion before splitting and after splitting 
  return(list(
    type = "node",
    split_var = best_split$var,
    threshold = best_split$threshold,
    left = left_subtree,
    right = right_subtree,
    aic = parent_model$aic,
    model = parent_model,
    data = data
  ))
}
#Example function call: 
#tree <- tree_garch(data, covariates = c("X1", "X2"), returns_col = "log_returns", min_size = 100)

#df_1 = read.csv("daily_data.csv")
df_1 <- read_excel("Data.xlsx", 
                   sheet = 2)

df_aux <- data.frame(
  date = index(log_ord),
  oil = coredata(log_ord)
)

df_1$Date <- as.Date(df_1$Date, format = "%Y-%m-%d")
names(df_1)[names(df_1) == 'Date'] <- 'date'
library(dplyr)

df_gas <- data.frame(
  date = index(log_grd),
  gas = coredata(log_grd)
)
# Assuming df1 and df2 both have a "date" column
df2 <- left_join(df_aux, df_1, by = "date")
df <- left_join(df2,df_gas,by = "date")
df = na.omit(df)
#train_data = filter(df, as.Date(df$date) <= "2025-01-01")
#train_data = df
train_data  <- df[1:(nrow(df) - 25), ]
tree <- tree_garch(train_data, covariates = c("GOLD", "SP500","DXY","10YRBOND","SP500_ENERGY","BLOOM_COMODITY_INDEX","GEOPOLITICAL_RISK_INDEX"), returns_col = "oil", min_node_size =300)
print(tree)
cat("Split variable:", tree$split_var, "\n")
cat("Threshold:", tree$threshold, "\n")
left_subtree <- tree$left
right_subtree <- tree$right

# Repeat analysis on left subtree
left_subtree$type
left_subtree$split_var  # if node
print_tree <- function(node, depth = 0) {
  indent <- paste(rep("  ", depth), collapse = "")
  if (node$type == "leaf") {
    cat(indent, "- Leaf node: AIC =", round(node$aic, 2), "\n")
  } else {
    cat(indent, "- Node: split_var =", node$split_var, ", threshold =", round(node$threshold, 3), ", AIC =", round(node$aic, 2), "\n")
    print_tree(node$left, depth + 1)
    print_tree(node$right, depth + 1)
  }
}
print(tree)

prune_tree_aic <- function(node, delta = 0) {
  # If leaf node, return as is
  if (node$type == "leaf") {
    return(node)
  }
  
  # Recursively prune children first
  node$left <- prune_tree_aic(node$left, delta)
  node$right <- prune_tree_aic(node$right, delta)
  
  # Calculate sum of children's AICs
  n_left <- nrow(node$left$data)
  n_right <- nrow(node$right$data)
  total_n <- n_left + n_right
  
  children_aic_sum <- node$left$aic + node$right$aic
  
  parent_aic <- node$aic
  
  # If parent's AIC is better (or close within delta), prune children
  if (parent_aic <= children_aic_sum + delta) {
    # Turn node into leaf: remove children, keep parent's model and data
    node$type <- "leaf"
    node$model <- node$model
    node$aic <- parent_aic
    node$left <- NULL
    node$right <- NULL
  }
  
  return(node)
}
pruned_tree <- prune_tree_aic(tree, delta = 0)  # delta is a tolerance parameter
print_tree(pruned_tree)

get_leaf_nodes <- function(tree) {
  if (tree$type == "leaf") return(list(tree))
  c(get_leaf_nodes(tree$left), get_leaf_nodes(tree$right))
}

leaf_nodes <- get_leaf_nodes(pruned_tree)
total_tree_aic <- pruned_tree$aic

print_tree(pruned_tree)

total_tree_aic

# Recursive function to assign conditional volatility
assign_volatility <- function(node) {
  if (node$type == "leaf") {
    preds <- volatility(node$model$model)  # Extract conditional volatility
    node$data$volatility <- preds
    return(node$data)
  } else {
    left_data <- assign_volatility(node$left)
    right_data <- assign_volatility(node$right)
    return(rbind(left_data, right_data))
  }
}
library(ggplot2)

# Assign volatility estimates
vol_data <- assign_volatility(pruned_tree)

# Order by date (optional)
vol_data <- vol_data[order(vol_data$date), ]
train_vol = vol_data$volatility
# Plot conditional volatility
ggplot(vol_data, aes(x = date, y = volatility)) +
  geom_line(color = "black") +
  labs(
    title = "Conditional Volatility Tree GARCH",
    x = "Time",
    y = "Conditional Volatility"
  ) +
  scale_y_continuous(limits = c(0, NA)) +  # Start y-axis at 0
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center and bold title
    axis.title.x = element_text(margin = margin(t = 10)),   # Add space below x-axis title
    axis.title.y = element_text(margin = margin(r = 10))    # Add space beside y-axis title
  )

####### Get GARCH Coefficients
leaf_nodes <- get_leaf_nodes(pruned_tree)

trace_leaf_paths <- function(node, path = list(), results = list()) {
  if (node$type == "leaf") {
    results[[length(results) + 1]] <- list(
      path = path,
      leaf = node
    )
    return(results)
  }
  
  # Add current split to path
  split_info <- list(var = node$split_var, threshold = node$threshold)
  
  # Traverse left and right subtrees
  results <- trace_leaf_paths(node$left,  append(path, list(paste0(node$split_var, " <= ", round(node$threshold, 4)))), results)
  results <- trace_leaf_paths(node$right, append(path, list(paste0(node$split_var, " > ",  round(node$threshold, 4)))), results)
  
  return(results)
}


# Step 1: Trace paths from the root to each leaf
leaf_paths <- trace_leaf_paths(pruned_tree)

# Step 2: Extract GARCH coefficients and metadata for each leaf + its path
leaf_coefs <- lapply(leaf_paths, function(leaf_info) {
  leaf <- leaf_info$leaf
  path <- leaf_info$path
  
  coefs <- coef(leaf$model$model)
  omega <- coefs["omega"]
  alpha <- coefs["alpha1"]
  beta <- coefs["beta1"]
  gamma <- coefs["gamma1"]  
  shape <- coefs["shape"] 
  
  n <- length(leaf$data$oil)
  eps <- leaf$data$oil
  
  # Approximate unconditional variance for TGARCH
  uncond_var <- omega / (1 - alpha - beta - gamma / 2)
  
  list(
    coefficients = coefs,
    n_obs = n,
    omega = omega,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    uncond_var = uncond_var,
    eps = eps,
    path = path,
    shape = shape
  )
})

# Step 3: Print them
for (i in seq_along(leaf_coefs)) {
  cat("Leaf", i, "- n =", leaf_coefs[[i]]$n_obs, "\n")
  cat("Path:", paste(leaf_coefs[[i]]$path, collapse = " → "), "\n")
  print(leaf_coefs[[i]]$coefficients)
  cat("\n")
}

# Compute the conditional volatility using GARCH(1,1)
compute_garch_volatility <- function(eps, omega, alpha, beta) {
  n <- length(eps)
  sigma2 <- numeric(n)
  
  # Initialize with unconditional variance
  sigma2[1] <- omega / (1 - alpha - beta)
  
  for (t in 2:n) {
    sigma2[t] <- omega +
      alpha * eps[t - 1]^2 +
      beta * sigma2[t - 1]
  }
  
  return(sqrt(sigma2))  # Return conditional volatility
}


# Filter test set
test_data <- train_data  <- df[(nrow(df) - 25):nrow(df), ]

# Assume test_data has a full series of past returns to compute rolling volatility
# Initialize output vector
predicted_vol <- rep(NA, nrow(test_data))

# Loop through test_data
for (i in 2:nrow(test_data)) {
  row <- test_data[i, ]
  eps_series <- test_data[1:i, "oil"]  # returns up to current time
  
  # Apply tree splits with updated parameters including 'shape'
  if (row$10YRBOND <= 3.9703 && row$DXY <= 92.896) {
    mu     = 2.960216e-03
    omega  = 7.155247e-05
    alpha  = 1.509411e-01
    beta   = 6.535480e-01
    shape  = 6.192218e+00
  } else if (row$10YRBOND <= 3.9703 && row$DXY > 92.896 && row$GOLD <= 1492.51) {
    mu     = 1.722160e-03
    omega  = 2.401725e-05
    alpha  = 8.554585e-02
    beta   = 9.159468e-01
    shape  = 2.694742e+00
  } else if (row$10YRBOND <= 3.9703 && row$DXY > 92.896 && row$GOLD > 1492.51 && row$SP500_ENERGY <= 583.7595) {
    mu     = 5.176838e-04
    omega  = 9.213131e-05
    alpha  = 2.132494e-01
    beta   = 7.572345e-01
    shape  = 3.386012e+00
  } else if (row$10YRBOND <= 3.9703 && row$DXY > 92.896 && row$GOLD > 1492.51 && row$SP500_ENERGY > 583.7595) {
    mu     = 7.392844e-04
    omega  = 6.959928e-05
    alpha  = 4.570180e-02
    beta   = 8.435972e-01
    shape  = 1.000000e+01
  } else {  # 10YRBOND > 3.9703
    mu     = 1.422442e-04
    omega  = 1.243412e-05
    alpha  = 6.086418e-02
    beta   = 9.015197e-01
    shape  = 1.000000e+01
  }
  
  
  
  eps <- eps_series - mu
  vol_series <- compute_garch_volatility(eps, omega, alpha, beta)
  
  # Store only the latest volatility estimate
  predicted_vol[i] <- tail(vol_series, 1)
}



# Add result to test_data
test_data$predicted_vol <- predicted_vol
train_vol = vol_data$volatility

# Create unified data frames with consistent structure
train_df <- data.frame(date = vol_data$date, volatility = train_vol, type = "Train")
test_df  <- data.frame(date = test_data$date, volatility = test_data$predicted_vol, type = "Test")

# Combine them
combined_vol <- rbind(train_df, test_df)

library(ggplot2)

ggplot(combined_vol, aes(x = date, y = volatility, color = type)) +
  geom_line() +
  labs(
    title = "Train vs Forecasted Volatility",
    y = "Conditional Volatility",
    x = "Date",
    color = "Data Type"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggplot(subset(combined_vol, type == "Test"), aes(x = date, y = volatility^2)) +
  geom_line() +
  labs(
    title = "Forecasted Variance",
    y = "Conditional Variance",
    x = "Date"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Compute real and predicted values
real <- test_data$oil^2
predicted <- combined_vol[combined_vol$type == "Test", "volatility"]
predicted = predicted^2
# Calculate MSE
mse <- mean((real - predicted)^2, na.rm = TRUE)
mse

log_pred  <- log(predicted, base = exp(1)) 
qlike = mean(log_pred + (real/predicted), na.rm = TRUE)
qlike

write.csv(combined_vol, "predictions_tree.csv", row.names = FALSE)








#################################################################################################
#################################################################################################
#################################################################################################
# APPENDIX:
#################################################################################################
#################################################################################################
#################################################################################################


#################################################################################################
# GARCH WITH NORMAL INNOVATIONS
#################################################################################################




# Window 1: GARCH(1,1)
garch1 = garchFit(~garch(1,1),data=log_ord_1, include.mean = FALSE, cond.dist = "norm", 
                  trace=F)
summary(garch1)
par(mfrow=c(2,2))
a=time(log_ord_1)
plot(zoo(garch1@residuals/garch1@sigma.t,a),xlab="",ylab="",main="GARCH(1,1) residuals")
qqnorm(garch1@residuals/garch1@sigma.t)
acf(garch1@residuals/garch1@sigma.t,main="GARCH(1,1) residuals")
acf(garch1@residuals^2/garch1@sigma.t^2,main="GARCH(1,1) squared residuals")
jarque.bera.test(garch1@residuals) #reject

# Window 2: GARCH(1,2)
garch2 = garchFit(~garch(1,2),data=log_ord_2, , include.mean = TRUE, cond.dist = "norm", 
                  trace=F)
summary(garch2)
par(mfrow=c(2,2))
a=time(log_ord_2)
plot(zoo(garch2@residuals/garch2@sigma.t,a),xlab="",ylab="",main="GARCH(1,2) residuals")
qqnorm(garch2@residuals/garch2@sigma.t)
acf(garch2@residuals/garch2@sigma.t,main="GARCH(1,2) residuals")
acf(garch2@residuals^2/garch2@sigma.t^2,main="GARCH(1,2) squared residuals")
jarque.bera.test(garch2@residuals) #reject

# Window 3: GARCH(1,1)
garch3 = garchFit(~garch(1,1),data=log_ord_3, include.mean = TRUE, cond.dist = "norm", 
                  trace=F)
summary(garch3)
par(mfrow=c(2,2))
a=time(log_ord_3)
plot(zoo(garch3@residuals/garch3@sigma.t,a),xlab="",ylab="",main="GARCH(1,1) residuals")
qqnorm(garch3@residuals/garch3@sigma.t)
acf(garch3@residuals/garch3@sigma.t,main="GARCH(1,1) residuals")
acf(garch3@residuals^2/garch3@sigma.t^2,main="GARCH(1,1) squared residuals")
jarque.bera.test(garch3@residuals) #reject




#################################################################################################
# URGARCH PACKAGE
#################################################################################################


#Uncomment any of the plot commands to get a menu of available plots

spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std",  # Student's t-distribution
  #fixed.pars = list(shape = 4)
)
fit <- ugarchfit(spec = spec, data = log_ord)
show(fit)            # Summary output
coef(fit)            # Model coefficients
#plot(fit)            # Diagnostic plots

#EGARCH
spec <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std",  # Student's t-distribution
  #fixed.pars = list(shape = 5)
)
fit <- ugarchfit(spec = spec, data = log_ord)
show(fit)            # Summary output
coef(fit)            # Model coefficients
#plot(fit)            # Diagnostic plots

#TGARCH
spec_tgarch <- ugarchspec(
  variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0)),
  distribution.model = "std",  # Student's t-distribution
  #fixed.pars = list(shape = 5)
)
fit <- ugarchfit(spec = spec_tgarch, data = coredata(log_ord))
show(fit)            # Summary output
coef(fit)            # Model coefficients
#plot(fit)            # Diagnostic plots




###############################################################################################################
###############################################################################################################

###############################################################################################################
###############################################################################################################



# Import the data
data <- read_excel("Data.xlsx", 
                        sheet = 3, skip = 0)



#select a short window
data <- data[data$Date > as.Date("2021-01-01"),]


# Daily Oil Price
data$Date <- as.Date(data$Date)

oil <- zoo(data$log_oil, order.by = data$Date)
gas <- zoo(data$log_gas, order.by = data$Date)


oil <- na.omit(oil)
gas <- na.omit(gas)


# Individual GARCH models

spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std",  # Student's t-distribution
  #fixed.pars = list(shape = 4)
)


garch_oil <- ugarchfit(spec = spec, data = oil)
garch_gas <- ugarchfit(spec = spec, data = gas)


#standardized residulas
z_oil <- residuals(garch_oil, standardize = TRUE)
z_gas <- residuals(garch_gas, standardize = TRUE)

z_df <- data.frame("z_oil" = coredata(z_oil), "z_gas" = coredata(z_gas))


#volatility vector(s)
h_oil <- coredata(sigma(garch_oil))
h_gas <- coredata(sigma(garch_gas))
h = data.frame(h_oil, h_gas)

# Fix correlation matrix
R_hat <- cor(z_df)
R_hat


h_11 <- c(rep(NA,nrow(h)))
h_22 <- c(rep(NA,nrow(h)))

for (i in 1:nrow(h)){
  D_t <- diag(h[i,])
  H_t <- D_t * R_hat * D_t
  h_11[i] = H_t[1,1]
  h_22[i] = H_t[2,2]
  
}

date <- time(gas)
h_22 <- zoo(h_22, date)
h_11 <- zoo(h_11, date)

plot(h_22, main = "Gas - CCC-GARCH(1,1)", xlab= "Time", ylab="Conditional variance")
plot(h_11, main = "Oil - CCC-GARCH(1,1)", xlab= "Time", ylab="Conditional variance")



######## DCC garch



# Define univariate GARCH specification (same for both series)
uspec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0)),
  distribution.model = "std"
)

# Combine into a multivariate DCC-GARCH spec
spec <- dccspec(
  uspec = multispec(replicate(2, uspec)),
  dccOrder = c(1, 1),
  distribution = "mvnorm"
)

returns <- merge(oil, gas)

# Fit the model
fit <- dccfit(spec, data = returns)


R_t <- rcor(fit)

# Plot dynamic correlation over time
cor_ts <- zoo(R_t[1, 2, ], order.by = index(returns))
plot(cor_ts, type = "l", main = "Dynamic Conditional Correlation: Oil vs Gas", ylab = "Correlation")

# Conditional covariances
H_t <- rcov(fit)  # [asset1, asset2, time]

sigma_t <- sigma(fit)

hdcc_11 <- zoo(sigma_t$oil, date)
hdcc_22 <- zoo(sigma_t$gas, date)



plot(hdcc_22^2, main = "Gas - DCC-GARCH(1,1)", xlab= "Time", ylab="Conditional variance")
plot(hdcc_11^2, main = "Oil - DCC-GARCH(1,1)", xlab= "Time", ylab="Conditional variance")



#############################################################################################################
#############################################################################################################
# Realized VOLATILITY
#############################################################################################################
#############################################################################################################

library(dplyr)
library(lubridate)


RV <- read_excel("RV.xlsx", sheet = 2)
RV$Date <- as.Date(RV$Date)

# 1 min
one_min <- RV %>% 
  group_by(Date) %>% 
  summarise(RV = sum(`r^2`, na.rm = TRUE)) 
one_min  


# 3 min
RV <- read_excel("RV.xlsx", sheet = 3)
RV$Date <- as.Date(RV$Date)

three_min <- RV %>% 
  group_by(Date) %>% 
  summarise(RV = sum(`r^2`, na.rm = TRUE)) 
three_min <- three_min[-7,]


# 5-min
RV <- read_excel("RV.xlsx", sheet = 4)
RV$Date <- as.Date(RV$Date)

five_min <- RV %>% 
  group_by(Date) %>% 
  summarise(RV = sum(`r^2`, na.rm = TRUE)) 
five_min 


# 10-min
RV <- read_excel("RV.xlsx", sheet = 5)
RV$Date <- as.Date(RV$Date)

ten_min <- RV %>% 
  group_by(Date) %>% 
  summarise(RV = sum(`r^2`, na.rm = TRUE)) 
ten_min 

# 30-min
RV <- read_excel("RV.xlsx", sheet = 6)
RV$Date <- as.Date(RV$Date)

halfh_min <- RV %>% 
  group_by(Date) %>% 
  summarise(RV = sum(`r^2`, na.rm = TRUE)) 
halfh_min 


# 1h
RV <- read_excel("RV.xlsx", sheet = 7)
RV$Date <- as.Date(RV$Date)

one_hour <- RV %>% 
  group_by(Date) %>% 
  summarise(RV = sum(`r^2`, na.rm = TRUE)) 
one_hour 


# 2h
RV <- read_excel("RV.xlsx", sheet = 8)
RV$Date <- as.Date(RV$Date)

two_h <- RV %>% 
  group_by(Date) %>% 
  summarise(RV = sum(`r^2`, na.rm = TRUE)) 
two_h 


df <- data.frame("date" = one_min[,1])
df[,"min1"] <- one_min[,2]
df[,"min3"] <- three_min[,2]
df[,"min5"] <- five_min[,2]
df[,"min10"] <- ten_min[,2]
df[,"min30"] <- halfh_min[,2]
df[,"h1"] <- one_hour[,2]
df[,"h2"] <- two_h[,2]


df <- df[,-1]

vec <- c(rep(NA,7))




library(tidyr)
library(dplyr)
library(ggplot2)

# 1. Compute mean RV for each frequency
mean_rv <- df %>%
  summarise(across(everything(), ~mean(.x, na.rm = TRUE))) %>%
  pivot_longer(cols = everything(), names_to = "frequency", values_to = "mean_rv")

# 2. Convert frequency labels to numeric time in minutes for plotting
mean_rv <- mean_rv %>%
  mutate(freq_min = case_when(
    frequency == "min1" ~ 1,
    frequency == "min3" ~ 3,
    frequency == "min5" ~ 5,
    frequency == "min10" ~ 10,
    frequency == "min30" ~ 30,
    frequency == "h1" ~ 60,
    frequency == "h2" ~ 120
  ))

# 3. Plot the signature plot
ggplot(mean_rv, aes(x = freq_min, y = mean_rv)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = mean_rv$freq_min) +
  labs(
    x = "Sampling Frequency (minutes)",
    y = "Mean Realized Volatility",
    title = "Signature Plot (Mean RV by Sampling Frequency)"
  ) +
  theme_minimal()



































































