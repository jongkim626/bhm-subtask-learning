################################################################################
###
### A Bayesian Hierarchical Model for a microgenetic analysis of subtask learning
###
###	Jong Kim (jongkim626@gmail.com)
###	Last Updated:  1 March 2019
###
################################################################################
###
###   Table of Contents
###
###   1. Learning Curves at the Global Level (Figure 1)
###   2. Analytics of Subtasks Data
###   3.  Toward linear mixed effects model
###   4.  Set up the environment for Bayesian computations
###   5.  The Models
###############################################################################

## load library

library(ggplot2)
library(rstan)

## set working directory
setwd("~/rdata/subtask-learning")

################################################################################
###
###	1.  Learning Curves at the Global Level (Figure 1)
###
################################################################################






## load data 
herbal.data <- read.csv("~/rdata/subtask-learning/herbal-data.csv", header = T)
attach(herbal.data)
names(herbal.data)
head(herbal.data)

## a plot for learning at the global performance 
ggplot(herbal.data, aes(x = Day, y = Time, shape = Performance)) +
  geom_line(data = subset(herbal.data, Performance %in% c("Novice",
                                                           "Expert20",
                                                           "Expert50",
                                                           "Expert70")), size=1) +
  geom_point(data = subset(herbal.data, Performance %in% c("Novice",
                                                            "Expert20",
                                                           "Expert50",
                                                           "Expert70")), size=3) +
  #HumanData
  geom_line(data = subset(herbal.data, Performance %in% "HumanData"), size = 2.5, color="red") +
  geom_point(data = subset(herbal.data, Performance %in% "HumanData"), size = 4) +
  geom_errorbar(data = subset(herbal.data, Performance %in% "HumanData"), 
                aes(ymin = Time-c(55.60, 31.23, 33.81, 26.14), 
                    ymax=Time+c(55.60, 31.23, 33.81, 26.14)), 
                width = 0.3, linetype="dotted") +
  #KLM
  geom_line(data = subset(herbal.data, Performance %in% "KLM"), size = 1, color="blue") +
  geom_point(data = subset(herbal.data, Performance %in% "KLM"), size = 3) +
  
  scale_x_continuous(breaks = c(1:10)) +
  scale_y_continuous(limits = c(100, 1500), breaks = seq(0, 1500, 200)) +
  ylab (label = "Task Completion Time (s)") +
  xlab (label = "Trials (Day)")+
  ggtitle("Learning Curves by Human and Model vs. the KLM Prediction") +
  theme(plot.title = element_text(hjust=0.5)) + 
  scale_shape_discrete(name="Performance",
                       breaks=c("Novice", "Expert20", "Expert50", "Expert70", "HumanData", "KLM"),
                       labels=c("Model - Novice", "Model 20% Expert", "Model 50% Expert", "Model 70% Expert", "Human", "KLM Prediction"))

  
################################################################################
####
####  2. Subtask Learning
####
################################################################################



###
###  2.1. Setting up the data
###


## Preparing data

subtask.learning.k <- read.table("~/rdata/subtask-learning/sl-data.txt", header=T)

attach(subtask.learning.k)
names(subtask.learning.k)
str(subtask.learning.k) # subtask is int, need to change it into factor

subtask.learning.k$subtask <- factor(subtask.learning.k$subtask, levels = c(1:14),
                                   labels =c("S1:FileOpen",
                                             "S2:SaveAs",
                                             "S3:NormCalc",
                                             "S4:Sum",
                                             "S5:FreqCalc",
                                             "S6:Sum",
                                             "S7:Length",
                                             "S8:TotalLength",
                                             "S9:TypdChar",
                                             "S10:TotalTypdChar",
                                             "S11:InsRows",
                                             "S12:Name",
                                             "S13:Date",
                                             "S14:SavePrn"))

subtask.learning.k$day <- factor(subtask.learning.k$day)

subtask.learning.k$time <- as.integer(subtask.learning.k$time)



# check the variables
str(subtask.learning.k)  # now subtask, day, pID are treated as a factor 
head(subtask.learning.k)
summary(subtask.learning.k)
# there some of the observations have missing values for certain covariates. 
# while we do not subset the data to only include complete cases to demonstrate that 
# automatically drops these observations, it is generally good practice to manually 
# do so if required.
subtask.learning.k <- na.omit(subtask.learning.k)

## looking at the data
names(subtask.learning.k)
str(subtask.learning.k)
levels(subtask.learning.k)
dim(subtask.learning.k)
head(subtask.learning.k)


## before Stan call
rstan_options(auto_write = TRUE) # autosave for no more recompilation
options(mc.cores = parallel::detectCores()) # parallel execution of MC


###############################################################################
###
###   Fixed Effect Model
###
################################################################################
stanDat <- list(time= subtask.learning.k$time,
                day = as.integer(subtask.learning.k$day),
                N = nrow(subtask.learning.k))

fixEfFit <- stan(file = "BHM-kimR19-fixed.stan", 
                 data = stanDat, 
                 iter = 2000, chains = 4)

# The function generates four chanis of samples.
# A Markov chanin is a stochastic process in which random values are sequentially
# generated. Each of the four chains contains 2000 samples of each parameter.

# Each chain is independent from the others.
# if all chains are converged to the same region of the parameter space, 
#  it is more likely that they converge to the posterior distributions. 

## plot traceplot

traceplot(fixEfFit, pars = c("beta","sigma_e"), inc_warmup = FALSE)
## choose beta and sigma parameters and  omit the warmup samples

# If traceplot looks like a fat hairy caterpillar (Lunn et al, 2012) that 
# does not ben, this suggests that the chains have converged to the poterior
# distributions.  

# if the Rhat statistic = 1, if the chain is converged


## print, the model summary, print(fixEfFit)

print(fixEfFit, pars = c("beta","sigma_e"), probs = c(0.025, 0.5, 0.975))

## Each parameter has the Rhat statistic, the ratio of btwn-chain variance to
## within-chain variance

## Make a table including credible intervals and Rhat statistic
## parameter mean 2.5% 97.5% Rhat
##

## examine quantiles of parameter of interests
## the chains have converged, now let's look at the posterior distribution
beta1 <- unlist(extract(fixEfFit, pars = "beta[2]"))
print(quantile(beta1, probs = c(0.025, 0.5, 0.975)))

## save
save(list="fixEfFit",file="fixEfFit.Rda",compress="xz")

################################################################################
##
##  Varying Intercept Mixed Model
##
################################################################################

### setting up the data for stan
## setting up the data for stan
stanDat <- list(pID = as.integer(subtask.learning.k$pID),
                subtask = as.integer(subtask.learning.k$subtask),
                time = subtask.learning.k$time,
                day = as.integer(subtask.learning.k$day),
                N = nrow(subtask.learning.k),
                J = nlevels(subtask.learning.k$pID),
                K = nlevels(subtask.learning.k$subtask))


## Sample from Posterior Distribution
varIntFit <- stan(file = "BHM-kimR19-VarInt.stan",
                 data = stanDat,
                 iter = 2000,
                 chains = 4)

## Summarize results
traceplot(varIntFit, pars = c("beta", "sigma_e", "sigma_u", "sigma_w"), 
          inc_warmup=FALSE)
print(varIntFit, pars = c("beta", "sigma_e", "sigma_u", "sigma_w"),
      probs = c(0.025, 0.5, 0.975))

beta1 <- unlist(extract(varIntFit, pars = "beta[2]"))
print(quantile(beta1, probs = c(0.025, 0.5, 0.975)))

u <- unlist(extract(varIntFit, pars = "u"))
print(quantile(u, probs = c(0.025, 0.5, 0.975)))

w <- unlist(extract(varIntFit, pars = "w"))
print(quantile(w, probs = c(0.025, 0.5, 0.975)))
## Posterior Probability of Beta1 being less than 0
mean(beta1 < 0)
mean(u<0)
mean(w<0)
###############################################################################
###
###   Varying Intercept and Varying Slopes Mixed Model
###
###############################################################################

## Setting up the data for stan
stanDat <- list(time = subtask.learning.k$time, 
                pID = as.integer(subtask.learning.k$pID),
                subtask = as.integer(subtask.learning.k$subtask), 
                day = as.integer(subtask.learning.k$day), 
                N = nrow(subtask.learning.k),
                J = nlevels(subtask.learning.k$pID),
                K = nlevels(subtask.learning.k$subtask))




# compile and fit the model
ranIntSlpNoCorFit <- stan(file = "BHM-kimR19-VarIntVarSlp.stan", 
                          data = stanDat,
                          iter = 2000, chains =4)

## traceplot
traceplot(ranIntSlpNoCorFit, pars = c("beta", "sigma_e", "sigma_u", "sigma_w"), 
          inc_warmup=FALSE)


print(ranIntSlpNoCorFit, pars = c("beta", "sigma_e", "sigma_u", "sigma_w"),
      probs = c(0.025, 0.5, 0.975))

beta1 <- unlist(extract(ranIntSlpNoCorFit, pars = "beta[2]"))

print(quantile(beta1, probs = c(0.025, 0.5, 0.975)))



###
###  Consideration of Correlation 
###

### setting up data for stan
stanDat <- list(time = subtask.learning.k$time, 
                pID = as.integer(subtask.learning.k$pID),
                subtask = as.integer(subtask.learning.k$subtask), 
                day = as.integer(subtask.learning.k$day), 
                N = nrow(subtask.learning.k),
                J = nlevels(subtask.learning.k$pID),
                K = nlevels(subtask.learning.k$subtask))

## fit model
VarIntVarSlpCorFit <- stan(file = "BHM-kimR19-VarIntVarSlpCorFit.stan", 
                          data = stanDat,
                          iter = 2000, chains =4)



## Use the L matrices to compute the corrleation matrices 
# L matrices
L_u <- extract(VarIntVarSlpCorFit, pars = "L_u")$L_u
L_w <- extract(VarIntVarSlpCorFit, pars = "L_w")$L_w

# correlation parameters
cor_u <- apply(L_u, 1, function(x) tcrossprod(x)[1,2])
cor_w <- apply(L_w, 1, function(x) tcrossprod(x)[1,2])

print(signif(quantile(cor_u, probs = c(0.025, 0.5, 0.975)), 2))
print(mean(cor_u))
print(signif(quantile(cor_w, probs = c(0.025, 0.5, 0.975)), 2))
print(mean(cor_w))

## traceplot
traceplot(VarIntVarSlpCorFit, pars = c("beta", "sigma_e", "sigma_u", "sigma_w"), 
          inc_warmup=FALSE)


print(VarIntVarSlpCorFit, pars = c("beta", "sigma_e", "sigma_u", "sigma_w"),
      probs = c(0.025, 0.5, 0.975))

## posterior probability of beta 1 being less than 0;
beta0 <- unlist(extract(VarIntVarSlpCorFit, pars="beta[1]"))
print(quantile(beta0, probs = c(0.025, 0.5, 0.975)))
beta1 <- unlist(extract(VarIntVarSlpCorFit, pars="beta[2]"))
print(quantile(beta1, probs = c(0.025, 0.5, 0.975)))
mean(beta1 <0)

#postsig <- rstan::extract(VarIntVarSlpCorFit, pars=c("beta","z_u", "z_w", "sigma_e", "sigma_u","sigma_w"))
#ref <- melt(postsig)
#colnames(ref)
#ggplot(data=ref, aes(x=u[1], y=u[2]))+geom_density()



################################################################################
####   End of the code for Bayesian Inference 
################################################################################



###
###
###  Using rstanarm 
###
###

## I need more details of model specifications by using lme4 style notations

# (1|day) + (1|subtask), varying intercepts 
# x + (x | g)
# 1+x + (1+ x | g) 

# 1. 
# (1 + A | B)
# a varying slope and varying intercept model for B
# to allow the slope of A variable to vary by B variable 

#2
# (1 | A) + (1 | B)
# Varying intercepts 

# days, subjects, subtasks
# multiple measurements

# (1 + subtasks|days)



# time ~ days + subtask + subjects + (1 + subtasks | days)

# time ~ days + subtask + subjects + (1 + subjects | days)

# time ~ days + subtask + subjects + (1 + subtasks | days) + (1 + subjects | days) 




###############################################################################
###
###	  Basic plots: Scatter, mean, KLM prediction
###
################################################################################

##  Scatter plot
p <- ggplot(subtask.learning.k, aes(x = day, y = time)) +
  geom_point(colour="black", size=0.5, na.rm=F) +
  facet_wrap(~subtask, ncol=7) 

print(p)

## Scatter plot + linear line
p1 <- p + labs(title="Subtask Completion Time", x ="Trials(Day)", y="Time(s)") +
  geom_smooth(method=lm, na.rm = F, formula =y~x, se=F)

print(p1)

###
### 2.3. Stat function for plot (retrieved from cookbook-r)
###

# Gives count, mean, standard deviation, standard error of the mean, and 
# confidence interval (default 95%).
# data: a data frame.
# measurevar: the name of a column that contains the var. to be summariezed
# groupvars: a vector containing names of columns that contain grouping var.
# na.rm: a boolean indicating whether to ignore NA's
# conf.interval: the percent range of the conf. interval (default is 95%)

summarySE <- function(data=subtask.learning.k, measurevar, 
                      groupvars=subtask, na.rm=T,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## A data set (subtask, time, sd, se, ci)
new.data <- summarySE(subtask.learning.k, measurevar = "time", 
                 groupvars = c("subtask", "day"))
###  new.data
##             subtask day  N      time         sd         se        ci
## 1        S1:FileOpen   1 30  39.03733  13.625058  2.4875838  5.087680
## 2        S1:FileOpen   2 30  33.61400  12.961345  2.3664070  4.839846
## 3        S1:FileOpen   3 30  28.04633  11.211934  2.0470098  4.186605
## 4        S1:FileOpen   4 30  22.52733  10.849987  1.9809276  4.051452


##	Plot for mean with SE
p.mean <- ggplot(new.data, aes(x=day, y=time))+
  facet_wrap(~subtask, ncol=7)+
  geom_errorbar(aes(ymin=time-se, ymax=time+se), width=.3) +
  geom_line() +
  geom_point(size=2)+
  ggtitle("Subtask Completion Time")+
  theme(plot.title= element_text(hjust=0.5))+
  xlab("Trials (day)")+
  ylab("Time (s)")

print(p.mean)

##	Plot with KLM-GOMS Prediction
klm.data <- data.frame(z= c(rep(19.96,4), rep(16.27,4), rep(101.26,4), 
                              rep(18.09,4), rep(106.43,4), rep(22.79,4), 
                              rep(143.83,4), rep(18.09,4), rep(141.42,4), 
                              rep(18.09,4), rep(21.03,4), rep(6.93,4),
                              rep(16.68,4), rep(15.80,4)), 
                         new.data)
##klm.data
##         z           subtask day  N      time         sd         se        ci
## 1   19.96       S1:FileOpen   1 30  39.03733  13.625058  2.4875838  5.087680
## 2   19.96       S1:FileOpen   2 30  33.61400  12.961345  2.3664070  4.839846
## 3   19.96       S1:FileOpen   3 30  28.04633  11.211934  2.0470098  4.186605
## 4   19.96       S1:FileOpen   4 30  22.52733  10.849987  1.9809276  4.051452
## 5   16.27         S2:SaveAs   1 23  44.79913  25.185488  5.2515373 10.891022
## 6   16.27         S2:SaveAs   2 23  45.70870  23.832127  4.9693420 10.305785
## 7   16.27         S2:SaveAs   3 23  24.79391  13.090497  2.7295573  5.660755
## 8   16.27         S2:SaveAs   4 23  21.44652   9.492224  1.9792656  4.104746
## 9  101.26       S3:NormCalc   1 30 262.13867  74.808910 13.6581759 27.934106

p.mean.klm <- p.mean +
    geom_hline(aes(yintercept = z), color="red", linetype="dashed", klm.data)


##### Plot for mean, SE, and KLM #### (Figure 3 in ICCM16)  ####################
print(p.mean.klm)   ############################################################
################################################################################

##<<------------------------ some explorations ---------------------------------

##	Scatter plot in log-log
p.log <- p + scale_x_log10() + scale_y_log10() +
  #coord_trans(x="log10", y="log10") +
  labs(x ="logDay", y="logTime")
#labs(title="Subtask Completion Time in log-log Coordinates", 
#     x ="logDay", y="logTime")

print(p.log)

## with linear line
p.log1 <- p.log+stat_smooth(method=lm, na.rm = F, se=F)
print(p.log1)

##	Add equations, r^2 (This is not done yet)
# a data set of coefficients from lme analysis
#a=intercept, b =slope

lme.coeff <- data.frame(a = c(1.596300, 1.652294, 2.401673,
                              1.891016, 2.354749, 1.563019,
                              2.216870, 1.570739, 2.332884,
                              1.461958, 1.885475, 1.438642,
                              1.442136, 1.699970),
                        b = c(-0.4336350, -0.5494246, -0.7054928,
                              -1.2025624, -0.6239880, -0.6751516,
                              -0.4764150, -0.6212238, -0.5114052,
                              -0.4342365, -0.8917299, -0.9272653,
                              -0.6211262, -0.5312155),
                        type = c("S1:FileOpen",
                                 "S2:SaveAs",
                                 "S3:NormCalc",
                                 "S4:Sum",
                                 "S5:FreqCalc",
                                 "S6:Sum",
                                 "S7:Length",
                                 "S8:TotalLength",
                                 "S9:TypdChar",
                                 "S10:TotalTypdChar",
                                 "S11:InsRows",
                                 "S12:Name",
                                 "S13:Date",
                                 "S14:SavePrn"))

require(plyr)
eq <- substitute(italic(y) == a + b%.%itallic(x),
                 list(a=format(lme.coeff[1], digits=2),
                      b=format(lme.coeff[2], digits=2)))
as.character(as.expression(eq))
lme.coeff$group <-c(1:14)

#fix bug#########################################
eq1 <- ddply(lme.coeff, .(group, type), eq)

pt <- p.log1 + 
  geom_text(data=eq,aes(x = 2, y = 420,label=V1), 
                        parse = TRUE, inherit.aes=FALSE, size=3, 
                        colour = "gray")

#####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++############
  



########up to here ######

q1 <- q + geom_smooth(method = lm, na.rm = T)
print(q1)

q2 <- p.log + geom_abline(aes(intercept=a, slope=b), 
                      color="blue", linetype="dashed",na.rm=T,lme.coeff)+
  facet_wrap(~type,ncol=7)
print(q2)

# this is wrong, wrong.... 
q3 <- p.log + geom_line()
print(q3)

# regression line only
r <- ggplot(data=)

test <- ggplot(new.data, aes(x=day, y=time, group=subtask)) +
  facet_wrap(~subtask, ncol=7)+geom_point()+
  coord_trans(x="log10", y="log10")

## this is wrong...
p.log+geom_abline(aes(intercept = a, slope = b), data = lme.coeff)

####-------------------------------------------------------------------------->>


###
###	2.3.  Toward linear mixed effects model (Kim & Ritter, 2016)
###  


## library
library(lme4)
library(devtools)
library(arm)
##  Check with missing values in time column
which(is.na(subtask.learning.k$time))
# > 48 missing values

##  Check the assumption

# plot with ordinary least square regression line
plot(subtask.learning.k$day, subtask.learning.k$time, pch=16, 
     col=rgb(0,0,204,102,maxColorValue=255))
olsLine <- lm(subtask.learning.k$time ~ subtask.learning.k$day)
abline(olsLine, col="red")

## check the coorelation
summary((olsLine))
# THE SUMMARY OF THE BASIC O.L.S. REGRESSION SUGGESTS THAT 
# THERE IS A STATISICALLY SIGNIFICANT CORRELATION BETWEEN time and day

## check the normality assumption
qqnorm(residuals(olsLine))
qqline(residuals(olsLine))

# --> Q-Q Plot suggest the residuals are not normally distributed

## log transform
subtask.learning.k$logDay <- log10(subtask.learning.k$day)
subtask.learning.k$logTime <- log10(subtask.learning.k$time)
ols2 <- lm(subtask.learning.k$logTime ~ subtask.learning.k$logDay)
qqnorm(residuals(ols2))
qqline(residuals(ols2))


##	Exploring some preliminary linear mixed models 

## fixed effects model
mod <- lm(logTime~logDay, data = subtask.learning.k)
display(mod) # not working    
mod.glm <- glm(logTime~logDay, data=subtask.learning.k)
display(mod.glm) # not working 
AIC(mod.glm)

mod.glm1 <- glm(logTime~logDay+subtask, data=subtask.learning.k)
display(mod.glm1)
AIC(mod.glm1)
## the subtask effect greatly improves the model fit.  
anova(mod.glm, mod.glm1, test = "F")

mod.glm2 <- glm(logTime~logDay+pID, data=subtask.learning.k)
display(mod.glm2)
AIC(mod.glm2)

table(subtask.learning.k$pID, subtask.learning.k$subtask)
## a perfectly balanced design with 4 observations in each combination of 
## participants and subtasks

## interaction between participants and subtasks
mod.glm3 <- glm(logTime~logDay+pID:subtask, data=subtask.learning.k)
display(mod.glm3)
AIC(mod.glm3)

## what if we want to undrestand both the effect of the participants and 
## the effect of the class, as well as the effect of the participants and subtasks?
## this is not easily done with the standard glm.
mod.glm4 <- glm(logTime~logDay+pID*subtask-1, data=subtask.learning.k)
display(mod.glm4)
AIC(mod.glm4)

### Fit a varying intercept model with lmer

## (1|subtask) tells lmer to fit a linear model with a varying intercept group 
## effect using the variable "subtask"

## (1|pID) tells lmer to fit a linear model with a varying intercept group
## effect using the variable "pID"


## varying intercepts by pID and varying intercepts by subtask
mod0 <-lmer(logTime~(1|pID)+(1|subtask), data = subtask.learning.k, REML = F)
summary(mod0)
coef(mod0)
coef(mod0)$subtask[1]

## varying intercepts by pID
mod1 <- lmer(logTime~(1|pID), data = subtask.learning.k)
summary(mod1)
coef(mod1)
coef(mod1)$pID[1]


## including fixed effects
#install.packages("dplyr")
#library(dplyr)

mod2 <- lmer(logTime~logDay+(1|pID), data = subtask.learning.k, REML = F)
summary(mod2)
coef(mod2)

mod3 <- lmer(logTime~logDay+(1|pID)+(1|subtask),
             data = subtask.learning.k, REML = F)
summary(mod3)
coef(mod3)

anova(mod0, mod3)

## varying intercepts for pID, and varying intercepts and slopes for subtask
## to allow the slope of the logDay variable to vary by subtask
mod4 <- lmer(logTime~logDay+(1|pID)+(1+logDay|subtask),
             data = subtask.learning.k, REML = F) 
summary(mod4)
coef(mod4)

## varying slops and varying intercepts for pID and subtask
mod5 <- lmer(logTime~logDay+(1+logDay|pID)+(1+logDay|subtask),
             data = subtask.learning.k, REML = F)

summary(mod5)
coef(mod5)

## compare
anova(mod4, mod5)
## Data: subtask.learning.k
## Models:
## mod4: logTime ~ logDay + (1 | pID) + (1 + logDay | subtask)
## mod5: logTime ~ logDay + (1 + logDay | pID) + (1 + logDay | subtask)
##      Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)    
## mod4  7 -1116.0 -1078.2 565.02  -1130.0                             
## mod5  9 -1134.7 -1086.2 576.37  -1152.7 22.701      2  1.177e-05 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## lmerTest
#install.packages("lmerTest")
library(lmerTest)

mod6 <- lmer(logTime~logDay+pID+(1+logDay|subtask),
             data = subtask.learning.k, REML = F)

summary(mod6)
coef(mod6)

anova(mod6, mod4)

## result 2

mod7 <- lmer(logTime~logDay+(1+logDay|pID)+(1|subtask),
             data = subtask.learning.k, REML = F)

summary(mod7)
anova(mod4, mod7)



qqnorm(residuals(mod1))
qqline(residuals(mod1))

summary(mod1)
## no p-value in "summary"

## coeffcients
coef(mod1)

## item effects
subtask.learning.k$subtask = factor(subtask.learning.k$subtask)

##
##  Plot A box plot describing the distributions of the task completion times 
##  showing minimum, maximum, median values,
##

ggplot(subtask.learning.k, aes(subtask, y=time, colour=subtask)) +
  geom_boxplot()+scale_x_discrete(breaks=c("S1:FileOpen",
                                           "S2:SaveAs",
                                           "S3:NormCalc",
                                           "S4:Sum",
                                           "S5:FreqCalc",
                                           "S6:Sum",
                                           "S7:Length",
                                           "S8:TotalLength",
                                           "S9:TypdChar",
                                           "S10:TotalTypdChar",
                                           "S11:InsRows",
                                           "S12:Name",
                                           "S13:Date",
                                           "S14:SavePrn"),
                                  labels=c("1","2","3","4","5","6","7",
                                           "8","9","10","11","12","13","14"))
                                  
                                  #labels=c("S1","S2","S3","S4","S5","S6","S7",
                                  #         "S8","S9","S10","S11","S12","S13","S14"))

mod7 <- lmer(logTime~logDay+(1|pID)+(1|subtask), data = subtask.learning.k)
summary(mod7)

mod2 <- lmer(logTime~logDay+(1|pID), data = subtask.learning.k)
summary(mod2)
anova(mod7, mod2, refit=F)
coef(mod7)
ranef(mod7)

#varying slopes
d_bysubtask = na.omit(subtask.learning.k) %>%
  group_by(subtask, day) %>%
  summarise(mean_time = mean(time))

ggplot(d_bysubtask, aes(x=day, y=mean_time, 
                         colour=subtask, group=subtask, label=subtask)) +
  geom_line() + geom_point(shape=21, fill="white") +geom_text()

mod7b <- lmer(logTime~logDay+(1|pID)+(1 + logDay|subtask), data = subtask.learning.k)
summary(mod7b)
coef(mod7b)

## Result 2
anova(mod7, mod7b, refit=F)

# interpretation
# 1. To assess the significance of practice trials (day) as a predictor, we 
# looked at the t-value of the fixed effects. The t-value of the slope estimate 
# is large enough, we can estimate that the predictor is significant
# since our dataset is fairly large with 1680 observations. 
# 
# 2. In our model, the intercepts and slopes are not correlated by
# the covariates of either participants and subtasks.  ????



##	Understanding operators (slopes, mops, keys)

#Slope <- lme.coeff[2]

#(Intercept)     logDay
#S1:FileOpen          1.596632 -0.4345600
#S2:SaveAs            1.657985 -0.5653091
#S3:NormCalc          2.401562 -0.7052524
#S4:Sum               1.891906 -1.2038078
#S5:FreqCalc          2.354736 -0.6240242
#S6:Sum               1.564564 -0.6782297
#S7:Length            2.217044 -0.4769718
#S8:TotalLength       1.569021 -0.6169571
#S9:TypdChar          2.332999 -0.5118074
#S10:TotalTypdChar    1.461769 -0.4342316
#S11:InsRows          1.885252 -0.8910907
#S12:Name             1.438460 -0.9266803
#S13:Date             1.441696 -0.6203430
#S14:SavePrn          1.701905 -0.5368610

smk <- data.frame(Slope=c(-0.4336350,-0.5494246,-0.7054928,-1.2025624, 
                          -0.6239880, -0.6751516, -0.4764150, -0.6212238,
                          -0.5114052, -0.4342365, -0.8917299, -0.9272653,
                          -0.6211262, -0.5312155), 
                  MOps=c(3,3,20,4,20,4,39,4,40,4,2,2,4,3), 
                  KOps=c(33,26,158,27,169,37,194,27,186,27,39,9,24,25))
attach(smk)
smk
# add subtask
smk$subtask <- c("S1:FileOpen",
                 "S2:SaveAs",
                 "S3:NormCalc",
                 "S4:Sum",
                 "S5:FreqCalc",
                 "S6:Sum",
                 "S7:Length",
                 "S8:TotalLength",
                 "S9:TypdChar",
                 "S10:TotalTypdChar",
                 "S11:InsRows",
                 "S12:Name",
                 "S13:Date",
                 "S14:SavePrn")

str(smk)

##
## Figure 5 in ICCM16
##

ggplot(smk, aes(MOps, Slope, KOps)) +
    geom_point()

smk.data <- read.csv("~/rdata/subtask-learning/slope.csv", header = T)

smk.data$Subtask <- factor(smk.data$Subtask, 
                           levels = c("S1","S2","S3","S4","S5","S6","S7",
                                      "S8","S9","S10","S11","S12","S13","S14"),
                           labels =c("1","2","3","4","5","6","7",
                                     "8","9","10","11","12","13","14"))

ggplot(smk.data, aes(Count, Slope, color=Ops, label=Subtask, shape=Ops)) +
  geom_point(size=4) +
  geom_text(vjust = 0, nudge_y = 0, colour="black", check_overlap = F, size=4.5)+
  labs(x ="Operator Count", y="Slope")





################################################################################
###
###	  3. Subtask Learning Analytics-- Bayesian Approach
###
###
################################################################################

##
##  3.1. Set up the environments-- packages and library
##

#devtools::install_github("rmcelreath/glmer2stan", force = T)
#Sys.setenv(MAKEFLAGS = "-j4") 
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)
#install.packages(c("Rcpp", "rstan"), type = "source")
#install.packages('lmerTest')
#install.packages("shinystan")

# verify toolchain is working
fx <- inline::cxxfunction( signature(x = "integer", y = "numeric" ) , '
	return ScalarReal( INTEGER(x)[0] * REAL(y)[0] ) ;
                           ' )
fx( 2L, 5 ) 
# the result should be 10

# library
library(devtools)
library(shinystan)
library(glmer2stan)
library(rstan)
library(lmerTest)

# for initial data inspection, using built-in plotting functions
#library(nlme)
library(lme4)
library(lattice)
library(ggplot2)

##
##   Load data
##

## Load the Dismal task learning data (Kim & Ritter, 2015)
## The data is 30 participants of keyboard users doing the Dismal task
## that consists of 14 subtasks. This will be up in the github some time later

## set the workning directory
setwd("~/rdata/subtask-learning/")

microDat <- subtask.learning.k

attach(microDat)
names(microDat)


###
###		 The Models 
###

##
##	No Random Effects Model (just regular linear model)
##

## Model 1 
m1_lm <- lm(logTime~logDay,data=microDat)
confint(m1_lm)
summary(m1_lm)

ggplot(microDat,aes(x=day,y=time))+
  geom_point()+
  guides(color=F)+
  geom_smooth(method=lm,se = F)

ggplot(microDat,aes(x=logDay,y=logTime))+
  geom_point()+
  guides(color=F)+
  geom_smooth(method=lm,se = F)

## with nesting of intercepts and slopes within subjects
ggplot(microDat,aes(x=day,y=time,color=pID,group=pID))+
  geom_point()+
  guides(color=F)+
  geom_smooth(method=lm,se = F)

## Model 2
m2_lm <- lm(logTime~logDay+subtask, data=microDat)
confint(m2_lm)
summary(m2_lm)

## model 3
## m3_lm <- lm(logTime~logDay+subtask+pID, data=microDat2)
## confint(m3_lm)
## summary(m3_lm)


## woops... I should use variables of "subj" and "subt"
## Convert pID and subtask to sequential integers
microDat$subj <- as.integer(as.factor(microDat$pID))
microDat$subt <- as.integer(as.factor(microDat$subtask))

## model 4
m4_lm <- lm(logTime~logDay+subt+subj, data=microDat)
confint(m4_lm)
summary(m4_lm)


##
##	Random effects model 
##


## Random slopes and random intercepts

## model 1
m1_lme4 <- lmer( logTime ~ logDay + (day | pID), microDat, REML=FALSE )
summary(m1_lme4)
confint(m1_lme4)

## AIC: Akaike Information Criterion, to compare the models 
AIC(m1_lm,m1_lme4)

## model 2
## (1|subtask) tells lmer to fit a linear model with a varying-intercept group effect
## using the variable "subtask"
m2_lme4 <- lmer( logTime ~ logDay + subtask + pID + (1|subtask) )

confint(m2_lme4)
summary(m2_lme4)
AIC(m3_lm, m2_lme4)
#        df      AIC
#m3_lm   33 18397.27
#m2_lme4 34 16674.39

AIC(m4_lm, m2_lme4)

m2_lme4_1 <- lmer( logTime ~ logDay + subt + subj + (1|subt), data=microDat2)
confint(m2_lme4_1)
summary(m2_lme4_1)
AIC(m2_lme4, m2_lme4_1)



# be careful- it is taking much time to compute 
#m3_lme4 <- lmer( time ~ day + subtask + pID + (1|subtask) + (1|pID) )
#confint(m3_lme4)
#summary(m3_lme4)
#AIC(m2_lme4, m3_lme4)

###  model
## varying slope of "day" by subtasks
## (1+day|subtask)
m4_lme4 <- lmer( logTime ~ logDay + subtask + pID + (1+day|subtask) )

confint(m4_lme4)
summary(m4_lme4)

m4_lme4_1 <- lmer( logTime ~ logDay + subt + subj + (1+logDay|subt), microDat )
confint(m4_lme4_1)
summary(m4_lme4_1)

AIC(m2_lme4_1, m4_lme4_1)


AIC(m1_lme4, m4_lme4_1)
AIC(m1_lme4, m2_lme4_1)


##
##	Using R and Stan
##

##	 Stan does not deal with non-numeric variables
## factors must converted to contrast codes (glmer2stan helps)
## grouping variables must be manually converted to integers

## Convert pID and subtask to sequential integers
microDat$subj <- as.integer(as.factor(microDat$pID))
microDat$subt <- as.integer(as.factor(microDat$subtask))

nwarm = 1000 # burn-in period, these samples are not included in estimation
niter = 2000 # number of steps per chain, more is better (but takes longer)
chains = 4 # number of chains, usually at least 2


## using m3_lm:  m3_lm <- lm(time~day+subtask+pID, data=microDat)

m3_lm_g2s <- lmer2stan(logTime ~ logDay + subt + subj, data=microDat,
                       calcWAIC=T,
                       warmup=nwarm,
                       iter=niter,
                       chains=chains)

print(m3_lm_g2s) # standard stan output
stanmer(m3_lm_g2s) # cleaned up stan output
plot(m3_lm_g2s) # *looks like shit
traceplot(m3_lm_g2s) ## more sampling


## using m2_lme4: m2_lme4 <- lmer( time ~ day + subtask + pID + (1|subtask) )

m2_lme4_g2s <- lmer2stan(logTime~logDay+subt+subj+(1|subt), data=microDat,
                calcWAIC=T,
                warmup=nwarm,
                iter=niter,
                chains=chains)

print(m2_lme4_g2s) # standard stan output
stanmer(m2_lme4_g2s) # cleaned up stan output
plot(m2_lme4_g2s) # *looks like shit
traceplot(m2_lme4_g2s) ## more sampling

## using m4_lme4: m4_lme4 <- lmer( time ~ day + subtask + pID + (1+day|subtask) )

m4_lme4_g2s <- lmer2stan(logTime~logDay+subt+subj+(1+logDay|subt), data=microDat,
                         calcWAIC=T,
                         warmup=nwarm,
                         iter=niter,
                         chains=chains)
print(m4_lme4_g2s) # standard stan output
stanmer(m4_lme4_g2s) # cleaned up stan output
plot(m4_lme4_g2s) # *looks like shit
traceplot(m4_lme4_g2s) ## more sampling


################################################################################


##	by-subject
m1_g2s <- lmer2stan( logTime ~ logDay + (logDay | subj), data=microDat,
                    calcWAIC=T,
                    warmup=nwarm, 
                    iter = niter, 
                    chains=chains) 

print(m1_g2s) # standard stan output
stanmer(m1_g2s) # cleaned up stan output
plot(m1_g2s) # *looks like shit
traceplot(m1_g2s) ## more sampling

##	by-subtask
m2_g2s <- lmer2stan( logTime ~ logDay + (logDay | subt), data=microDat,
                     calcWAIC=T,
                     warmup=nwarm, 
                     iter = niter, 
                     chains=chains) 
print(m2_g2s) # standard stan output
stanmer(m2_g2s) # cleaned up stan output
plot(m2_g2s) # *looks like shit
traceplot(m2_g2s) ## more sampling

################################################################################
##




