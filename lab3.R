# Setting Working Directory
setwd("~/Documents/Gothenburg University/Courses/MSA680 - Data Science for Bioscience/Labs/lab3")

#Load libraries
library(tidyverse)
library(ggplot2)
library(GGally)
library(survival)
library(survminer)

# Reading in data
df <- read.csv(file="whas500.csv",header=T, sep=";")
df$FSTAT <- as.factor(df$FSTAT)
df$GENDER <- as.factor(df$GENDER)

# Plotting histogram (all subjects)
ggplot(df, aes(x=LENFOL)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=30)+
  geom_density(alpha=.2, fill="#FF6666") +
  ggtitle("Histogram of Length of Follow up (days)") +
  labs(x="Length of Follow up (days)", y="Density")

# Plotting histogram (subjects who died)
ggplot(df[df$FSTAT == 1,], aes(x=LENFOL)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", binwidth=30)+
  geom_density(alpha=.2, fill="#FF6666") +
  ggtitle("Histogram of Length of Follow up (days) for Subjects who Died") +
  labs(x="Length of Follow up (days)", y="Density")

# Plotting cdf for death probability (dataset filtered for those who died)
ggplot(df[df$FSTAT == 1,], aes(x=LENFOL)) +
  stat_ecdf(geom = "step", pad=FALSE) +
  ggtitle("CDF of Probability of Death")+
  labs(x="Length of Follow up (days)", y="Probability of death")

# Calculating probability of surviving less than 200 days
df |> group_by(FSTAT) |> 
  summarize(less_than_200 = sum(LENFOL<200), at_least_200 = sum(LENFOL >=200), n = n())
# So p=114/215, roughly 53%

# Plotting the 1 - ecdf survival function estimate
ggplot(df[df$FSTAT == 1,], aes(x=LENFOL)) +
  geom_step(aes(y = 1 - ..y..), stat = "ecdf", pad=FALSE) +
  ggtitle("1 - ecdf Survival Function")+
  labs(x="Length of Follow up (days)", y="Probability of survival")

# Calculating probability of surviving at least 1000 days
df |> group_by(FSTAT) |> 
  summarize(less_than_1000 = sum(LENFOL<1000), at_least_1000 = sum(LENFOL >=1000), n = n())
# So p=36/215, roughly 16.7%

# First look at the data
df %>% group_by(FSTAT) %>% summarise(N=n(),
                                         Mean=mean(LENFOL),
                                         SD=sd(LENFOL),
                                         Min=min(LENFOL),
                                         Max=max(LENFOL))

df %>% group_by(FSTAT) %>% summarise(N=n(),
                                     Mean=mean(AGE),
                                     SD=sd(AGE),
                                     Min=min(AGE),
                                     Max=max(AGE))

df %>% group_by(FSTAT) %>% summarise(N=n(),
                                     Mean=mean(BMI),
                                     SD=sd(BMI),
                                     Min=min(BMI),
                                     Max=max(BMI))

df %>% group_by(FSTAT) %>% summarise(N=n(),
                                     Mean=mean(HR),
                                     SD=sd(HR),
                                     Min=min(HR),
                                     Max=max(HR))

# Pairs plots
ggpairs(df, columns = c(2,3,4,7,18), ggplot2::aes(colour = FSTAT))

##########
# Non-parametric analysis
##########
df$FSTAT <- as.numeric(df$FSTAT)
df$GENDER <- as.numeric(df$GENDER)

# Creating survival object
surv.data <- with(df, Surv(LENFOL, FSTAT))

# Kaplan-Meier estimate
fit <- survfit(surv.data~1, data = df)
summary(fit)
# Median time at t=1627, SE=0.02822
# [0.4471, 0.558] for 95% CI probability of survival
# [1527, NA] for 95% CI probability of median time

# Kaplan-Meier estimate, log-log confidence interval
fit2 <- survfit (surv.data~1, data = df, conf.type="log-log")
fit2
# SE=0.02822, [0.4429, 0.553]

# Plotting the KM survival curve
ggsurvplot(fit,
          risk.table = TRUE,
          surv.median.line = "hv",
          ggtheme = theme_bw())

# Plotting the cumulative hazard curve
ggsurvplot(fit,
          fun = "cumhaz",
          risk.table = TRUE,
          ggtheme = theme_bw ())

# Separate KM curves by gender
fit3 <- survfit (surv.data~GENDER, data = df)
summary(fit3)
ggsurvplot(fit3,
           risk.table = TRUE,
           surv.median.line = "hv",
           ggtheme = theme_bw())

ggsurvplot(fit3,
           fun = "cumhaz",
           risk.table = TRUE,
           ggtheme = theme_bw ())

# Log-rank test to compare survival between genders
surv.test <- survdiff(surv.data~GENDER, data = df )
surv.test
# Significant difference at p=0.005

##########
#Parametric analysis
##########
# A random example
# Setting up the dataframe
random_df <- data.frame (
  Patient = seq(1, 6),
  Survtime = c(6, 7, 10, 15, 19, 25),
  Censor = c(1, 0, 1, 1, 0, 1),
  Group = c("C", "C", "T", "C", "T", "T")
)

# Constructing survival object, setting up Cox PH fit
random_surv.data <- with(random_df, Surv(Survtime , Censor))
random_fit <- coxph(random_surv.data~Group, data=random_df)
summary(random_fit)
# GroupT coef=-1.3261

#Computing by hand...
loglikelihood <- function(beta) {
  return(beta - log(3*exp(beta)+3)-log(3*exp(beta)+1)-log(2*exp(beta)+1))
}
ggplot(data.frame(x = c(-4, 1)), aes(x = x)) + 
  geom_function(fun = loglikelihood, color = "red", size = 1) + 
  labs(x = "Beta", y = "Log-likelihood") +
  theme_bw()

result <- optim(par = 0, fn = loglikelihood, method = "L-BFGS-B",
                control = list(fnscale = -1),
                lower = -4, upper =1)
# Beta partial likelihood estimate = -1.326129, same as above!
# Value of partial log-likelihood = -3.671981


# Going back to the WHAS dataset...
fit4 <- coxph(surv.data~GENDER + AGE , data = df)
summary(fit4)
# Gender not significant, but age is

# Correcting by mean age, then graphing separate regression curves
sex.dt <- with(df,
              data.frame(GENDER = c("Male", " Female " ),
                        AGE = rep(mean(AGE, na.rm = TRUE), 2)))

surv.fit.gender <- survfit (fit4, data = df, newdata = sex.dt)

ggsurvplot(surv.fit.gender, 
           conf.int = TRUE,
           legend.labs = c("Male","Female"),
           ggtheme = theme_minimal())

# Trying all covariates in a full model...
# DSTAT excluded since wasn't sure whether it was death stat...
fit5 <- coxph(surv.data~ AGE + HR + SYSBP + BMI + CVD + AFB + SHO + 
                 CHF + AV3 + MIORD + MITYPE, data = df)
summary(fit5)
# Seems like age, hr, BMI, SHO, CHF significant...

fit6 <- coxph(surv.data~ AGE + HR + BMI + 
                SHO + CHF, data = df)
summary(fit6)

extractAIC(fit4)
extractAIC(fit5)
extractAIC(fit6)

