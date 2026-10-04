#################################################################################
######### CODE IMPLEMENTING THE DISCRETIZED-LOGIT-NORMAL DISTRIBUTION ###########
############################  accompanying the paper ############################
# A discrete version of the logit-normal distribution for modeling ordinal data #
#################################################################################

library(logitnorm)                # for the continuous logit-normal distribution
library(bbmle)                    # for maximum likelihood estimation
library(ggplot2)                  # for graphs of a bivariate PMF

library(reshape2)
# function for plotting a bivariate pmf
plot.pmf <- function(pmf)
{
  par(mar = c(0.1, 0.1, 0.1, 0.1))
  df<-melt(pmf)
  names(df) <- c("X1", "X2", "Probability")
  ggplot(df, aes(x = X1, y = X2, fill = Probability)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Probability)), 
            color = "black", size = 4) +
  scale_fill_gradient(low = "white", high = "navy") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12,  margin = margin(t = -10)),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(angle = 0, vjust  = 0.5,  hjust  = 1, margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 4)),
    plot.margin = margin(1, 1, 1, 1, "pt"),          # Reduce outer margins
    legend.margin = margin(0, 0, 0, 0, "pt"),        # Reduce legend margins
  ) +
  scale_x_continuous(expand = c(0, 0)) +             # Remove horizontal padding
  labs(
    #title = "",
    #subtitle = "",
    x = expression(X[1]),
    y = expression(X[2]),
    fill = "Probability"
  ) +
  coord_fixed() 
}
########################################################
############# CONTINUOUS LOGIT-NORMAL ##################
########################################################

#############################################
### plots of the pdf of the continuous LN  ##
# for various combinations of mu and sigma2 #
################# FIGURE 1 ##################
################# S T A R T #################
#############################################
op <- par()
par(mfrow=c(3,1), mai=c(0.3,0.3,0.1,0.1), lwd=2)
# mu=0
mu <- 0
sigma2 <- c(0.25, 0.5, 1, 2, 4)
n <- length(sigma2)
plot(c(0,0), c(1,1), type="n", xlim=c(0,1), ylim=c(0,4.5))
for(i in 1:n)
{
  plot(function(x) dlogitnorm(x, mu, sqrt(sigma2)[i]), add=TRUE, xlim=c(0,1), lty=i, col=i)
}
legend_text <- lapply(sigma2, function(s) bquote(sigma^2 == .(s))) 
legend("topleft", legend = legend_text,text.col=1:n, col = 1:n, lty=1:n,bty="n", cex=1.25)
legend("topright", legend = bquote(mu == .(mu)), bty="n", cex=1.25, inset = c(0.05, 0))
# mu=0.5
mu <- 1/2
plot(c(0,0), c(1,1), type="n", xlim=c(0,1), ylim=c(0,4.5))
for(i in 1:n)
{
  plot(function(x) dlogitnorm(x,mu,sqrt(sigma2)[i]), add=TRUE,xlim=c(0,1), lty=i, col=i)
}
legend_text <- lapply(sigma2, function(s) bquote(sigma^2 == .(s))) 
legend("topleft", legend = legend_text,text.col=1:n, col = 1:n, lty=1:n,bty="n", cex=1.25)
legend("topright", legend = bquote(mu == .(mu)), bty="n", cex=1.25, inset = c(0.05, 0))
# mu=1
mu <- 1
plot(c(0,0), c(1,1), type="n", xlim=c(0,1), ylim=c(0,4.5))
for(i in 1:n)
{
  plot(function(x) dlogitnorm(x,mu,sqrt(sigma2)[i]), add=TRUE,xlim=c(0,1), lty=i, col=i)
}
legend_text <- lapply(sigma2, function(s) bquote(sigma^2 == .(s))) 
legend("topleft", legend = legend_text, text.col=1:n,col = 1:n, lty=1:n, bty="n", cex=1.25)
legend("topright", legend = bquote(mu == .(mu)), bty="n", cex=1.25, inset = c(0.05, 0))
par(op)
#########################################
# plots of the pdf of the continuous LN #
################# E N D #################

# Numerical check of the formula of the expected value of the inverse
# of the continuous logit-normal random variable
mu    <- 0
sigma <- 1
x     <- rlogitnorm(10000, mu=mu, sigma=sigma)
mean(1/x)                      # MC mean
exp(-mu+sigma^2/2)+1           # Paper: Moments of the logit-normal distribution
# computing the expected value directly from the definition using integrate
integrate(function(x) 1/x*dlogitnorm(x,mu,sigma), lower=0.0001, upper=1)

# computing the expected value of a logit-normal rv using integrate
mu    <- 1
sigma <- seq(0.1,3,0.1)
for (i in 1:length(sigma))
{
res <- integrate(function(x) x*dlogitnorm(x,mu,sigma[i]),lower=0,upper=1)
print(res)
}

#################################################
###### Discrete Logit-Normal distribution #######
####### with parameters mu, sigma2, and k #######
#################################################
# pmf: ddlogitnorm returns the probabilities for all values from 1 to k
# We use sigma squared as the second parameter
ddlogitnorm <- function(mu, sigma2, k)
{
  i <- 1:k
  prob <- plogitnorm(i/k, mu, sqrt(sigma2)) - plogitnorm((i-1)/k, mu, sqrt(sigma2))
  prob
}
# cdf: pdlogitnorm returns the cumulative probabilities for all values from 1 to k
pdlogitnorm <- function(mu, sigma2, k)
{
  i <- 1:k
  prob <- plogitnorm(i/k, mu, sqrt(sigma2))
  prob
}
# quantile function
qdlogitnorm <- function(prob, mu, sigma2, k)
{
  ceiling(k*(1+exp(-(mu+sqrt(sigma2)*qnorm(prob))))^(-1))
}
# pseudo-random generation
rdlogitnorm <- function(n, mu, sigma2, k)
{
  u <- runif(n)
  qdlogitnorm(u, mu, sigma2, k)
}
# log-likelihood function; for ML estimation
log.lik.dLN <- function(mu, sigma2, k)
{
  # x here, data=list(x=...) in mle2
  elle <- -sum(log(ddlogitnorm(mu, sigma2, k)[x]))
  ifelse(elle==Inf, 2*length(x), elle)
}
# exact moments
M.dLN <- function(mu, sigma2, k)
{
  E <- k - sum(pdlogitnorm(mu, sigma2, k)[-k])
  V <- sum((1:k)^2*ddlogitnorm(mu, sigma2, k))-E^2
  return(c(E,V))
}
E.dLN <- Vectorize(function(mu, sigma2, k=k) M.dLN(mu,sigma2,k)[1])
V.dLN <- Vectorize(function(mu, sigma2, k=k) M.dLN(mu,sigma2,k)[2])
M.dLN(0, 1.5, 7)
E.dLN(0, 1.5, 7)
V.dLN(0, 1.5, 7)

########################################################
# plot of pmf for various combinations of mu and sigma #
########################  START  #######################
######################   FIGURE 2  #####################
########################################################
op <- par()
par(mfrow=c(3,1), mai=c(0.3,0.3,0.1,0.1), lwd=2)
#  mu=0
mu <- 0
sigma2 <- c(0.25, 0.5, 1, 2, 4)
n <- length(sigma2)
k <- 7 # number of categories
plot(c(1,k), c(0,0), type="n", xlim=c(1,k), ylim=c(0,0.5))
for(i in 1:n)
{
  points(1:k, ddlogitnorm(mu,sigma2[i],k), type="b", xlim=c(0,1), col=i,pch=i)
}
legend_text <- lapply(sigma2, function(s) bquote(sigma^2 == .(s))) 
legend("topleft", legend = legend_text,text.col=1:n, col = 1:n,bty="n",pch=1:n,cex=1.25)
legend("topright", legend = bquote(mu == .(mu)),bty="n",cex=1.25, inset = c(0.05, 0))
# mu=0.5
mu <- 0.5
plot(c(1,k),c(0,0),type="n",xlim=c(1,k),ylim=c(0,0.5))
for(i in 1:n)
{
  points(1:k, ddlogitnorm(mu,sigma2[i],k), type="b", xlim=c(0,1), col=i,pch=i)
}
legend_text <- lapply(sigma2, function(s) bquote(sigma^2 == .(s))) 
legend("topleft", legend = legend_text,text.col=1:n, col = 1:n,bty="n", pch=1:n, cex=1.25)
legend("topright", legend = bquote(mu == .(mu)),bty="n",cex=1.25, inset = c(0.05, 0))
# mu=1
mu <- 1
sigma <- c(0.5, 0.75, 1, 1.5, 2)
n <- length(sigma)
plot(c(1,k), c(0,0), type="n", xlim=c(1,k), ylim=c(0,0.5))
for(i in 1:n)
{
  points(1:k, ddlogitnorm(mu,sigma2[i],k), type="b", xlim=c(0,1), col=i,pch=i)
}
legend_text <- lapply(sigma2, function(s) bquote(sigma^2 == .(s))) 
legend("topleft", legend = legend_text, text.col=1:n, col = 1:n, bty="n", pch=1:n, cex=1.25)
legend("topright", legend = bquote(mu == .(mu)), bty="n", cex=1.25, inset = c(0.05, 0))
par(op)
########################################################
########################   END   #######################
########################################################

##########################################################
# contour plots of E and V for the discrete logit-normal #
###################### ---> FIGURE 3 #####################
##########################################################
k <- 5
# Create a sequence of equally spaced values for mu and sigma2 in the range *(0,3) and (0,9)*
v.mu    <- seq(0.005, 3, by = 0.005)
v.sigma <- seq(0.005, 9, by = 0.005)
# Create a grid of points
zE <- outer(v.mu, v.sigma, E.dLN, k=k)
zV <- outer(v.mu, v.sigma, V.dLN, k=k)
# Plot the level curves of E and V
op<-par()
par(mai=c(0.7,0.7,0.25,0.15), mgp=c(2.25,1,0), mfrow=c(1,2), cex=1)
contour(v.mu, v.sigma, zE, xlab = expression(mu), ylab = expression(sigma^2),
        main = "Level curves of expected value, k=5", nlevels=10)
contour(v.mu, v.sigma, zV, xlab = expression(mu), ylab = expression(sigma^2),
        main = "Level curves of variance, k=5", nlevels=10)
par(op)
##########################################################
k <- 5
# changing the range of sigma2
v.mu    <- seq(0.005, 3, by = 0.005)
v.sigma <- seq(0.005, 1/2, by = 0.005)
# Create a grid of points
zE <- outer(v.mu, v.sigma, E.dLN, k=k)
zV <- outer(v.mu, v.sigma, V.dLN, k=k)
# Plot the level curves of E and V
op<-par()
par(mai=c(0.7,0.7,0.25,0.15), mgp=c(2.25,1,0), mfrow=c(1,2),cex=1)
contour(v.mu, v.sigma, zE, xlab = expression(mu), ylab = expression(sigma^2),
        main = "Level curves of expected value, k=5",nlevels=10)
contour(v.mu, v.sigma, zV, xlab = expression(mu), ylab = expression(sigma^2),
        main = "Level curves of variance, k=5",nlevels=10)
# two "poles" are visible at log(3/2) and log(4/1)
par(op)
##########################################################
################# PARAMETER  ESTIMATION ##################
##########################################################
# quadratic loss function for the method of moments
f.mom <- function(par,k,x)
{
  mu     <- par[1]
  sigma2 <- par[2]
  mom    <- M.dLN(mu, sigma2, k)
  (mean(x) - mom[1])^2 + (mean(x^2) - mom[1]^2-mom[2])^2
}
#########################################################
########### function for parameter estimation ###########
#########################################################
est.dLN <- function(x, k=max(x), method="ML")
{
  if(method=="ML")                               # maximum likelihood estimation
  {
  res <- mle2(minuslogl=log.lik.dLN, fixed=list(k=k), data=list(x=x),
              start=list(mu=mean(x)-(max(x)+1)/2, sigma2=1), method="L-BFGS-B",
              lower=c(-5,0.01), upper=c(5,25))
  return(res)
  }
  else if(method=="MM")                                      # method of moments
  {
    res <- optim(par=c(1,1), k=k, fn=f.mom, x=x, control=list(abstol=1e-10))
    return(res)
  }
  else if(method=="MP")                                  # method of proportions
  {
    f1 <- mean(x==1)
    fk <- mean(x==k)
    sigma.hat <- 2*log(k-1)/(qnorm(1-fk)-qnorm(f1))
    mu.hat    <- log(k-1)*(1-2*qnorm(1-fk)/(qnorm(1-fk)-qnorm(f1)))
    return(c(mu.hat=mu.hat,sigma2.hat=sigma.hat^2))
  }
}
# Examples
z <- c(1, 1, 1, 2, 2, 3)
est.dLN(z, method="ML")
est.dLN(z, method="MP")
est.dLN(z, method="MM")
z <- c(1, 1, 1, 2, 2, 3, 3, 3)
est.dLN(z, method="MM") # the estimated mu is not exactly zero!

set.seed(12345)
z <- rdlogitnorm(n=1000, mu=1, sigma2=2^2, k=5)
table(z)                                # observed frequencies
ddlogitnorm(mu=1, sigma=2^2, k=5)*1000  # expected frequencies
est.dLN(z, k=5, method="ML")
est.dLN(z, k=5, method="MM")
est.dLN(z, k=5, method="MP")
# the estimates of mu and sigma are quite close to each other across the methods
# and close to the true values 1 and 4!

# "uniform" samples,
# with k=5...
x <- rep(1:5)
est.5 <- est.dLN(x, method="ML")@coef
ddlogitnorm(est.5[1], est.5[2], 5)
# ... and k=7
x <- rep(1:7)
est.7 <- est.dLN(x, method="ML")@coef
ddlogitnorm(est.7[1], est.7[2], 7)
################################################################################

#####################################################
########## R E G R E S S I O N   M O D E L ##########
#####################################################
# directly on mu and ln(sigma^2)
n <- 1000
set.seed(12345)
x1 <- rnorm(n)
x2 <- runif(n,-1,1)
mu.0    <- 0
sigma.0 <- 1
lambda1 <- 1; lambda2 <- .75; omega1 <- -.25; omega2 <- .5
mu.v    <- mu.0 + lambda1 *x1 + lambda2 *x2
sigma2.v <- exp(sigma.0^2 + omega1 *x1  + omega2 *x2)
plot(mu.v, sigma2.v)
y <- rdlogitnorm(n, mu.v, sigma2.v, k=7)
plot(table(y))

ddiscreteLN_pmf <- Vectorize(function(x, k, mu, Lsigma,log=FALSE) {
  p <- ddlogitnorm(mu, exp(Lsigma), k)[x]
  return(ifelse(log==FALSE,p,log(p)))
})
x <- cbind(x1, x2)
mle_fit <- mle2(y~ddiscreteLN_pmf(mu, Lsigma, k=7),
                ## specify the distribution of the response y
                data=data.frame(y, x),             
                ## need to specify data as data.frame
                parameters=list(mu~x, Lsigma~x),   
                ## linear model for mu and log(sigma)
                start=list(mu=0, Lsigma=0))
summary(mle_fit)

###################################################
# reparametrization in terms of mu and exp(sigma) #
###  of the pmf of the discrete logit-normal  #####
###################################################

ddiscreteLN_pmf <- Vectorize(function(x, k, mu, log.sigma2,log=FALSE) {
  p <- ddlogitnorm(mu, exp(log.sigma2), k)[x]
  return(ifelse(log==FALSE, p, log(p)))
})

#################################################
# reparametrization in terms of delta and gamma #
###  of the pmf of the discrete logit-normal  ###
#################################################

# ddiscreteLN_pmf <- Vectorize(function(x, k, eta, delta,log=FALSE) {
#   p <- ddlogitnorm(eta*exp(delta), exp(delta), k)[x]
#   return(ifelse(log==FALSE, p, log(p)))
# })

#############################################################
################# D A T A   A N A L Y S I S #################
#############################################################

# (1) PISA #
library(likert)
library(CUB) # for the CUB model and Leti dissimilarity index
data(pisaitems)
# ST39Q07
x <- pisaitems$ST39Q07
x <- na.omit(x)
table(x)
levels <- c("Never", "A few times a year", "About once a month", 
            "Several times a month", "Several times a week")
x<-factor(x, ordered="TRUE", levels=levels)
x<-as.numeric(x)
table(x)
est.dLN(x, method="MP") # Method of Proportion
est.dLN(x, method="MM") # Method of Moments
res <- est.dLN(x, method="ML")
res@coef
# std errors
sqrt(diag(res@vcov))
prob <- ddlogitnorm(res@coef[1], res@coef[2], 5)
prob
AIC(res)
op<-par()
par(mai=c(0.4,0.4,0.1,0.1), mgp=c(1,0.75,0), cex=1.1)
bp <- barplot(table(x)/sum(table(x)), ylim=c(0,0.325), cex.axis=1.1)
points(bp, prob, pch=19, col="red")
box()
par(op)
dissim(prob/sum(prob), table(x)/sum(table(x)))
cbind(prob/sum(prob), table(x)/sum(table(x)))
# cumulative logit model
library(ordinal)
fit_clm <- clm(as.factor(x) ~ 1, link = "logit")
AIC(fit_clm) # cumulative link model is better

# (2) UTAH dataset #
library(readxl)
# reading from the Excel file in my directory
dataset <- read_excel("UTAH.xlsx")
########## MAIN INDUSTRIES ##########
y <- dataset[[7]] # Major industries ENVIRONMENT
index.na <- which(y=="MD")
# delete NA in y
dataset <- dataset[-index.na,]
dim(dataset)
y <- dataset[[7]]
y <- as.numeric(y)

res <- est.dLN(y, method="ML")
res@fullcoef
sqrt(diag(vcov(res)))
p <- ddlogitnorm(res@coef[1], res@coef[2], k=max(y))
f<- table(y)/sum(table(y))
k <- max(y)

op<-par()
par(mai=c(0.4,0.4,0.1,0.1), mgp=c(1,0.75,0), cex=1.1)
barplot(f)
points(1:8+seq(-0.3,1.1,0.2), p, col="red", type="b", pch=19)
par(op)

fit_clm <- clm(as.factor(y) ~ 1, link = "logit")
AIC(fit_clm) # cumulative link model is better
AIC(res)
BIC(fit_clm) 
BIC(res)     # the proposed mode model is better


# including covariates
levels <-c("Less than $25,000","$25,000 to $34,999","$35,000 to $49,999","$50,000 to $74,999","$75,000 to $99,999","$100,000 to $149,999","Greater than $150,000","Prefer not to answer")
X1 <- factor(dataset[[1]],ordered=FALSE,levels=levels)
X2 <- factor(dataset[[2]],levels=c("Democrat","Independent","Libertarian","Republican","Other","No preference"))
X <- data.frame(X1,X2)
# in y the dependent variable, which we assume to follow the discrete LN distribution
# in X the covariates
mle_fit <- mle2(y~ddiscreteLN_pmf(mu,log.sigma2,k=max(y)), ## specify the distribution of the response y
                data=data.frame(y,X),                      ## need to specify as data frame
                parameters=list(mu~unlist(X1)+unlist(X2),log.sigma2~unlist(X1)+unlist(X2)),
                ## linear model for mu and log.sigma2
                start=list(mu=0,log.sigma2=0))
summary(mle_fit) # OK

# NEW! Collapsing categories
# X1: <= 49,999; 50,000-99,999, >=100,000
# X2: "Democrat", "Republican", all the other categories
y <- dataset[[7]]
y <- as.numeric(y)
levels <-c("Less than $25,000","$25,000 to $34,999","$35,000 to $49,999","$50,000 to $74,999","$75,000 to $99,999","$100,000 to $149,999","Greater than $150,000","Prefer not to answer")
X1 <- factor(dataset[[1]],ordered=FALSE,levels=levels)
X2 <- factor(dataset[[2]],levels=c("Democrat","Independent","Libertarian","Republican","Other","No preference"))
library(forcats)
index <- which(X1=="Prefer not to answer")
X1 <- X1[-index]
X1 <- droplevels(X1)
X1 <- fct_collapse(
  X1,
  "Less than 50,000" = c("Less than $25,000","$25,000 to $34,999","$35,000 to $49,999"),
  "50,000-99,999" =c("$50,000 to $74,999","$75,000 to $99,999"),
  "At least 100,000"=c("$100,000 to $149,999","Greater than $150,000")
)
X1 <- relevel(X1, ref ="Less than 50,000")
X2 <- X2[-index]
X2 <- fct_collapse(
  X2,
  "Other"=c("Independent","Libertarian","Other","No preference"),
 "Democrat"="Democrat",
 "Republican"="Republican"
)
X2 <- relevel(X2, ref ="Other")
y <- y[-index]
X <- data.frame(X1,X2)
# in y the dependent variable, which we assume to follow the discrete LN distribution
# in X the covariates
# mle_fit <- mle2(y~ddiscreteLN_pmf(eta,delta,k=max(y)),     ## specify the distribution of the response y
#                 data=data.frame(y,X),                      ## need to specify as data frame
#                 parameters=list(eta~unlist(X1)+unlist(X2),delta~unlist(X1)+unlist(X2)),
                  ## linear model for eta and gamma
#                 start=list(eta=0,delta=0))
mle_U_fit <- mle2(y~ddiscreteLN_pmf(mu, log.sigma2, k=max(y)),## specify the distribution of the response y
                data=data.frame(y,X),                      ## need to specify as data frame
                parameters=list(mu~unlist(X1)+unlist(X2),log.sigma2~unlist(X1)+unlist(X2)),
                ## linear model for mu and log.sigma2
                start=list(mu=0, log.sigma2=0))
summary(mle_U_fit) # OK
AIC(mle_U_fit)


# expectation and variance of the discrete logit-normal distribution
# as functions of eta and gamma
E.dLNt <- function(eta, gamma, k)
{
  mu <- eta*exp(gamma); sigma=exp(gamma)
  E.dLN(mu,sigma,k)
}
V.dLNt <- function(eta, gamma, k)
{
  mu <- eta*exp(gamma); sigma=exp(gamma)
  V.dLN(mu,sigma,k)
}


# (3) ISSP #
# International Social Survey Programme: Health and Health Care II - ISSP 2021 #
# using the file ZA8000_v2-0-0.csv downloaded from ISSP #
dataframe <- read.csv("ZA8000_v2-0-0.csv")

#  V2 (confidence) and v6 (inefficiency health system)
y <- dataframe$v2[dataframe$c_alphan=="DK"] 
y <- y[y>0]
table(y)
plot(table(y))
res <- est.dLN(y, method="ML")
summary(res)
p <- ddlogitnorm(res@coef[1], res@coef[2], k=max(y))
AIC(res)
res <- GEM(Formula(y~0|0|0), family="cub")
summary(res)
f <- table(y)/sum(table(y))
f
p; abs(p-f)
dissim(p,f)
pcub <- as.numeric(fitted(res)); abs(pcub-f)
dissim(f,pcub)
op<-par()
par(mai=c(0.5,0.5,0.1,0.2))
barplot(table(y)/length(y), ylim=c(0,0.65))
# Q2
legend=c("1 Complete confidence","2 A great deal of confidence","3 Some confidence",
         "4 Very little confidence","5 No confidence at all")
legend(3.5, 0.65, legend=legend, bty="n")
points(1:5+seq(-0.3,0.6,0.2), p, type="b", 16)
points(1:5+seq(-0.3,0.6,0.2), pcub, type="b", pch=13, col="red")
legend(5.5,0.3, legend=c("dLN","CUB"), col=c("black","red"),pch=c(16,13),bty="n")
box()
par(op)

data <- dataframe[dataframe$c_alphan=="DK",]
data <- data[,c("v2","SEX","AGE","URBRURAL","DK_ISCD","DK_INC")] # attention to the country selected! IT or DK
data <- na.omit(data)
data <- data[apply(data >= 0, 1, all), ]
names(data)[names(data) == "DK_ISCD"] <- "EDU"
names(data)[names(data) == "DK_INC"] <- "INC"
data$INC <- factor(data$INC,levels=c("70000","125000","175000","225000","275000","350000","450000","550000","650000","750000","850000","950000","1050000","1150000"),
                   labels=c("medium-low","medium-low","medium-low","medium-low","medium-low","medium-low","medium-low","medium-low","medium-high","medium-high","medium-high","medium-high","medium-high","medium-high"))
data$EDU <- factor(data$EDU, levels=c("0","100","244","344","353","354","550","640","650","740","750","840"),
                   labels=c("Sec.","Sec.","Sec.","Sec.","Sec.","Sec.","Tert.","Tert.","Tert.","Tert.","Tert.","Tert."))
data$SEX      <- factor(data$SEX, levels=c("1","2"), labels=c("M", "F"))
data$URBRURAL <- factor(data$URBRURAL,levels=c("1","2","3","4","5"), labels=c("Urban","Urban","Non-urban","Non-urban","Non-urban"))
data$AGE <- factor(ifelse(data$AGE <= 60, "<=60", ">60"))

y <- data$v2

ddiscreteLN_pmf <- Vectorize(function(x, k, mu, log.sigma2,log=FALSE) {
  p <- ddlogitnorm(mu, exp(log.sigma2), k)[x]
  return(ifelse(log==FALSE, p, log(p)))
})

mle_fit <- mle2(y~ddiscreteLN_pmf(mu,log.sigma2,k=5),     ## specify the distribution of the response y
                data=data,                                 ## need to specify as data frame
                parameters=list(mu~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC),
                               log.sigma2~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC)),
                ## linear model for mu and log.sigma2
                start=list(mu=0,log.sigma2=0))
summary(mle_fit)
AIC(mle_fit)

mle_fit_red <- mle2(y~ddiscreteLN_pmf(mu,log.sigma2,k=5),  ## specify the distribution of the response y
                data=data,                                 ## need to specify as data frame
                parameters=list(mu~unlist(SEX)+unlist(AGE)+unlist(URBRURAL),
                                log.sigma2~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)),
                ## linear model for mu and log.sigma2
                start=list(mu=0,log.sigma2=0))
summary(mle_fit_red)
AIC(mle_fit_red)

res <- est.dLN(y, method="ML")
p <- ddlogitnorm(res@coef[1], res@coef[2], k=max(y))
AIC(res)
### OK


data$URBRURAL  <- fct_collapse(
  data$URBRURAL,
  "city" = c(1,2),
  "country"=c(3,4)
)
data$EDU<- fct_collapse(
  data$EDU,
  "secondary" = c(1,2),
  "tertiary"=c(3,4)
)
data$INC<- fct_collapse(
  data$INC,
  "medium-low" = c(1,2),
  "medium-high"=c(3,4)
)
ddiscreteLN_pmf <- Vectorize(function(x, k, mu, log.sigma2,log=FALSE) {
  p <- ddlogitnorm(mu, exp(log.sigma2), k)[x]
  return(ifelse(log==FALSE, p, log(p)))
})

mle_fit <- mle2(y~ddiscreteLN_pmf(mu,log.sigma2,k=5), ## specify the distribution of the response y
                data=data,                            ## need to specify as data frame
                parameters=list(mu~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC),
                                log.sigma2~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC)),
                ## linear model for mu and log.sigma2
                start=list(mu=0,log.sigma2=0))
summary(mle_fit)
AIC(mle_fit)

# V6 (Q4c)
# 
y <- dataframe$v6[dataframe$c_alphan=="DK"] 
y <- y[y>0]
table(y)
plot(table(y))
res <- est.dLN(y, method="ML")
summary(res)
p <- ddlogitnorm(res@coef[1], res@coef[2], k=max(y))
AIC(res)
res<-GEM(Formula(y~0|0|0), family="cub")
summary(res)
f <- table(y)/sum(table(y))
f
p; abs(p-f)
dissim(p,f)
pcub <- as.numeric(fitted(res)); abs(pcub-f)
dissim(pcub,f)

data <- dataframe[dataframe$c_alphan=="DK",]
data <- data[,c("v6","SEX","AGE","URBRURAL","DK_ISCD","DK_INC")] # attention to the country selected! IT or DK
data <- na.omit(data)
data <- data[apply(data >= 0, 1, all), ]
names(data)[names(data) == "DK_ISCD"] <- "EDU"
names(data)[names(data) == "DK_INC"] <- "INC"
data$INC <- factor(data$INC,levels=c("70000","125000","175000","225000","275000","350000","450000","550000","650000","750000","850000","950000","1050000","1150000"),
                   labels=c("medium-low","medium-low","medium-low","medium-low","medium-low","medium-low","medium-low","medium-low","medium-high","medium-high","medium-high","medium-high","medium-high","medium-high"))
data$EDU <- factor(data$EDU, levels=c("0","100","244","344","353","354","550","640","650","740","750","840"),
                   labels=c("Sec.","Sec.","Sec.","Sec.","Sec.","Sec.","Tert.","Tert.","Tert.","Tert.","Tert.","Tert."))
data$SEX      <- factor(data$SEX, levels=c("1","2"), labels=c("M", "F"))
data$URBRURAL <- factor(data$URBRURAL,levels=c("1","2","3","4","5"), labels=c("Urban","Urban","Non-urban","Non-urban","Non-urban"))
data$AGE <- factor(ifelse(data$AGE <= 60, "<=60", ">60"))

y <- data$v6
op<-par()
par(mai=c(0.5,0.5,0.1,0.2))
barplot(table(y)/length(y), ylim=c(0,0.65))
# Q4c
legend=c("1 Strongly agree","2 Agree","3 Neither agree nor disagree","4 Disagree","5 Strongly disagree")
legend(3.25,0.65,legend=legend,bty="n")
points(1:5+seq(-0.3,0.6,0.2), p, type="b", 16)
points(1:5+seq(-0.3,0.6,0.2), pcub, type="b", pch=13, col="red")
legend(0,0.3, legend=c("dLN","CUB"), col=c("black","red"),pch=c(16,13),bty="n")
box()
par(op)

mle_fit_4c <- mle2(y~ddiscreteLN_pmf(mu,log.sigma2, k=5),  ## specify the distribution of the response y
                data=data,                                 ## need to specify as data frame
                parameters=list(mu~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC),
                                log.sigma2~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC)),
                ## linear model for mu and log.sigma2
                start=list(mu=0,log.sigma2=0))
summary(mle_fit_4c)
AIC(mle_fit_4c)

mle_fit_4c_red <- mle2(y~ddiscreteLN_pmf(mu,log.sigma2, k=5), ## specify the distribution of the response y
                   data=data,                                 ## need to specify as data frame
                   parameters=list(mu~unlist(SEX)+unlist(AGE)+unlist(URBRURAL),
                  log.sigma2~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)), 
                  ## linear model for mu and log.sigma2
                   start=list(mu=0,log.sigma2=0))
summary(mle_fit_4c_red)
AIC(mle_fit_4c_red)

# without unlist, it still works
mle_fit_4c <- mle2(y~ddiscreteLN_pmf(mu,log.sigma2,k=max(y)),  ## specify the distribution of the response y
                data=data,                                 ## need to specify as data frame
                parameters=list(mu~SEX+AGE+URBRURAL+EDU+INC,log.sigma2~SEX+AGE+URBRURAL+EDU+INC),
                ## linear model for mu and log.sigma2
                start=list(mu=0,log.sigma2=0))
summary(mle_fit_4c)
AIC(mle_fit_4c)

res <- est.dLN(y, method="ML")
p <- ddlogitnorm(res@coef[1], res@coef[2], k=max(y))
AIC(res)




data$URBRURAL  <- fct_collapse(
  data$URBRURAL,
  "city" = c(1,2),
  "country"=c(3,4)
)
data$EDU<- fct_collapse(
  data$EDU,
  "secondary" = c(1,2),
  "tertiary"=c(3,4)
)
# change here! in order to have about 50%-50%
data$INC<- fct_collapse(
  data$INC,
  "medium-low" = c(1,2),
  "medium-high"=c(3,4)
)
mle_fit_4c <- mle2(y~ddiscreteLN_pmf(mu,log.sigma,k=5), ## specify the distribution of the response y
                   data=data,                           ## need to specify as data frame
                   parameters=list(mu~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC),
                                   log.sigma~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC)),
                   ## linear model for mu and log.sigma
                   start=list(mu=0,log.sigma=0))
summary(mle_fit_4c)

# If I want to use sigma^2 (or log(sigma^2)) as the second parameter
# p-values, ell.max, AIC and BIC will not change
mle_fit_4c <- mle2(y~ddiscreteLN_pmf(mu,log.sigma2/2,k=5), ## specify the distribution of the response y
                data=data,                                 ## need to specify as data frame
                parameters=list(mu~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC),
                                log.sigma2~unlist(SEX)+unlist(AGE)+unlist(URBRURAL)+unlist(EDU)+unlist(INC)), 
                ## linear model for mu and log.sigma2
                start=list(mu=0,log.sigma2=0))
summary(mle_fit_4c)


################################################################################
################## B I V A R I A T E   Discrete Logit-Normal ###################
############################# with GAUSS copula ################################
################################################################################
library(mvtnorm)
# joint continuous cdf for the bivariate continuous logit-normal rv
# with Gauss copula
FBLN <- function(x, y, mu1, sigma2.1, mu2, sigma2.2, rho)
{
sigma <- matrix(c(1, rho, rho, 1), 2, 2)
pmvnorm(lower=c(-Inf,-Inf), upper=c(qnorm(plogitnorm(x,mu1,sqrt(sigma2.1))),
                                    qnorm(plogitnorm(y,mu2,sqrt(sigma2.2)))),
                                    mean=c(0,0), sigma=sigma)
}
# joint cdf for the bivariate discrete logit-normal rv
FdBLN <- function(i, j, mu1=0, sigma2.1=1, mu2=0, sigma2.2=1 ,rho=0, k=5)
{
  i <- floor(i); j<- floor(j)
  ifelse(i>=0 & j>=0,FBLN(i/k, j/k, mu1, sigma2.1, mu2, sigma2.2, rho), 0)
}
# joint pmf for the bivariate discrete logit-normal rv
ddBLN <- function(i, j, mu1=0, sigma2.1=1, mu2=0, sigma2.2=1, rho=0, k=5)
{
  FdBLN(i, j, mu1, sigma2.1, mu2, sigma2.2, rho, k)+
    FdBLN(i-1, j-1, mu1, sigma2.1, mu2, sigma2.2, rho, k)-
    FdBLN(i-1, j, mu1, sigma2.1, mu2, sigma2.2, rho, k)-
    FdBLN(i, j-1, mu1, sigma2.1, mu2, sigma2.2, rho, k)
}
# graphs of the joint pmf
op<-par()
par(mar = c(0.1, 0.1, 0.1, 0.1))
pmf <- outer(1:5, 1:5, Vectorize(function(x,y) ddBLN(x, y, rho=.4)))
plot.pmf(pmf)
pmf <- outer(1:5, 1:5, Vectorize(function(x,y) ddBLN(x, y, mu1=-1, mu2=1, rho=.4)))
plot.pmf(pmf)
# log-likelihood function (old version, based on the sample (x,y))
log.lik.dBLN.old <- function(mu1, sigma2.1, mu2, sigma2.2, rho, k=5, x, y)
{
  neg.ell <- -sum(log(mapply(function(x, y) ddBLN(x, y, mu1, sigma2.1, mu2, sigma2.2, rho, k), x, y)))
  neg.ell <- ifelse(is.nan(neg.ell)|neg.ell==Inf, length(x)^2, neg.ell)
  neg.ell
}
x <- rdlogitnorm(100,0,1,5)
y <- rdlogitnorm(100,0.5,2,5)
log.lik.dBLN.old(0,1,0,1,0,k=5,x,y)
# log-likelihood function for the bivariate discrete logit-normal
# x.cell and y.cell denote the distinct observed category pairs
# n.cell denotes their corresponding counts
log.lik.dBLN <- function(mu1, sigma2.1, mu2, sigma2.2, rho, k = 5, x.cell, y.cell, n.cell) {
  p <- mapply(
    function(xi, yi)
      ddBLN(xi, yi, mu1, sigma2.1, mu2, sigma2.2, rho, k),
    x.cell, y.cell
  )
  if (any(!is.finite(p)) || any(p <= 0))
    return(sum(n.cell)^2)
  
  -sum(n.cell * log(p))
}
# example
tab <- as.data.frame(
  table(x = factor(x, levels = 1:5),
        y = factor(y, levels = 1:5))
)
tab <- subset(tab, Freq > 0)
x.c <- as.integer(as.character(tab$x))
y.c <- as.integer(as.character(tab$y))
n.c <- tab$Freq
log.lik.dBLN(0,1,0,1,0,k=5,x.c,y.c,n.c)

# data analysis using ISSP2021
data <- dataframe[dataframe$c_alphan=="DK",] 
data <- data[,c("v2","v6")]
data <- data[apply(data >= 0, 1, all), ]
data
x <- data[,1]
y <- data[,2]
tab <- as.data.frame(
  table(x = factor(x, levels = 1:5),
        y = factor(y, levels = 1:5))
)
tab <- subset(tab, Freq > 0)
x.c <- as.integer(as.character(tab$x))
y.c <- as.integer(as.character(tab$y))
n.c <- tab$Freq
res.biv <- mle2(log.lik.dBLN,start=list(mu1=0, sigma2.1=1, mu2=0, sigma2.2=1, rho=cor(data)[1,2]),
                fixed=list(k=5), data= list(
                  x.cell = x.c,
                  y.cell = y.c,
                  n.cell = n.c
                ),
                method="L-BFGS-B",
                lower=c(mu1=-Inf, sigma2.1=1e-4, mu2=-Inf, sigma2.2=1e-4, rho=-1),
                upper=c(mu1=Inf, sigma2.1=Inf, mu2=Inf, sigma2.2=Inf, rho=1),
     control=list(trace=TRUE))
summary(res.biv)
AIC(res.biv)
# in case of independence
res.biv.indep <- mle2(log.lik.dBLN,start=list(mu1=0, sigma2.1=1, mu2=0, sigma2.2=1, rho=cor(data)[1,2]),
                fixed=list(k=5, rho=0), data= list(
                  x.cell = x.c,
                  y.cell = y.c,
                  n.cell = n.c
                ),
                method="L-BFGS-B",
                lower=c(mu1=-Inf, sigma2.1=1e-4, mu2=-Inf, sigma2.2=1e-4, rho=-1),
                upper=c(mu1=Inf, sigma2.1=Inf, mu2=Inf, sigma2.2=Inf, rho=1),
                control=list(trace=TRUE))
summary(res.biv.indep)
AIC(res.biv.indep)


# joint relative frequencies
fij <- table(x,y)/length(x)
# joint probabilities for the bivariate discrete logit-normal model
pij <- matrix(0, 5, 5)
for(i in 1:5)
{
for(j in 1:5)
{
pij[i,j] <- ddBLN(i, j, res.biv@coef[1], res.biv@coef[2], res.biv@coef[3], res.biv@coef[4],
                  res.biv@coef[5], k=5)
}
}
pij
1/2*sum(abs(pij-fij))
# diagnostics
library(copula)
set.seed(12345)
n <- length(x)
B <- 1000
TV <- numeric(B)
for(i in 1:B)
{
ub <- rCopula(n=length(x), normalCopula(param=res.biv@coef[5]))
cor(u)
xb <- qdlogitnorm(ub[,1],res.biv@coef[1], res.biv@coef[2], k=5)
yb <- qdlogitnorm(ub[,2],res.biv@coef[3], res.biv@coef[4], k=5)
tab.b <- as.data.frame(
  table(xb = factor(xb, levels = 1:5),
        yb = factor(yb, levels = 1:5))
)

tab.b <- subset(tab.b, Freq > 0)

x.c <- as.integer(as.character(tab.b$x))
y.c <- as.integer(as.character(tab.b$y))
n.c <- tab.b$Freq

res.biv.b <- mle2(log.lik.dBLN, start=list(mu1=0, sigma2.1=1, mu2=0, sigma2.2=1, rho=cor(xb,yb)),
                  fixed=list(k=5),
                  data  = list(
                    x.cell = x.c,
                    y.cell = y.c,
                    n.cell = n.c
                  ), method="L-BFGS-B",
                  lower=c(mu1=-Inf, sigma2.1=1e-4, mu2=-Inf, sigma2.2=1e-4, rho=-1),
                  upper=c(mu1=Inf, sigma2.1=Inf, mu2=Inf, sigma2.2=Inf, rho=1),
                  control=list(trace=TRUE))
pij.b <- matrix(0, 5, 5)
for(h in 1:5)
{
  for(k in 1:5)
  {
    pij.b[h,k] <- ddBLN(h, k, res.biv.b@coef[1], res.biv.b@coef[2], res.biv.b@coef[3],
                        res.biv.b@coef[4], res.biv.b@coef[5], k=5)
  }
}
TV[i] <- 1/2*sum(abs(pij - pij.b))
print(i)
print(TV[i])
}
summary(TV)

# adding the covariate AGE

# The bivariate pmf, with mu1 e mu2 depending on AGE (here 0/1)
ddBLN.AGE <- function(i, j, AGE,
                      beta10, beta11, sigma2.1,
                      beta20, beta21, sigma2.2,
                      rho, k = 5) {
  
  ddBLN(
    i = i, j = j,
    mu1 = beta10 + beta11 * AGE,
    sigma2.1 = sigma2.1,
    mu2 = beta20 + beta21 * AGE,
    sigma2.2 = sigma2.2,
    rho = rho, k = k
  )
}


# (negative) log-likelihood
log.lik.dBLN.AGE <- function(
    beta10, beta11, sigma2.1,
    beta20, beta21, sigma2.2,
    rho, k = 5, x.cell, y.cell, age.cell, n.cell) {
  
  p <- mapply(
    function(xi, yi, ai)
      ddBLN.AGE(
        i = xi, j = yi, AGE = ai,
        beta10 = beta10, beta11 = beta11,
        sigma2.1 = sigma2.1,
        beta20 = beta20, beta21 = beta21,
        sigma2.2 = sigma2.2,
        rho = rho, k = k
      ),
    x.cell, y.cell, age.cell
  )
  
  if (any(!is.finite(p)) || any(p <= 0))
    return(1e100)
  
  -sum(n.cell * log(p))
}


# data
data <- dataframe[
  which(dataframe$c_alphan == "DK"),
  c("v2", "v6", "AGE")
]
# data as factor
data$AGE <- factor(ifelse(data$AGE <= 60, 0, 1))
# only complete data
data <- data[
  complete.cases(data) &
    data$v2 %in% 1:5 &
    data$v6 %in% 1:5,
]
# frequencies also stratified by AGE
tab <- as.data.frame(
  table(
    x = factor(data$v2, levels = 1:5),
    y = factor(data$v6, levels = 1:5),
    AGE = factor(data$AGE, levels = 0:1)
  )
)

tab <- subset(tab, Freq > 0)

x.c   <- as.integer(as.character(tab$x))
y.c   <- as.integer(as.character(tab$y))
age.c <- as.integer(as.character(tab$AGE))
n.c   <- tab$Freq

# Estimation
res.biv.age <- mle2(
  log.lik.dBLN.AGE,
  start = list(
    beta10 = -0.633, beta11 = 0, sigma2.1 = 0.365,
    beta20 =  0.933, beta21 = 0, sigma2.2 = 2.142,
    rho = -0.423
  ),
  fixed = list(k = 5),
  data = list(
    x.cell = x.c,
    y.cell = y.c,
    age.cell = age.c,
    n.cell = n.c
  ),
  method = "L-BFGS-B",
  lower = c(
    beta10 = -Inf, beta11 = -Inf, sigma2.1 = 1e-4,
    beta20 = -Inf, beta21 = -Inf, sigma2.2 = 1e-4,
    rho = -0.999
  ),
  upper = c(
    beta10 = Inf, beta11 = Inf, sigma2.1 = Inf,
    beta20 = Inf, beta21 = Inf, sigma2.2 = Inf,
    rho = 0.999
  ),
  control = list(trace = TRUE, maxit = 2000)
)

summary(res.biv.age)
AIC(res.biv.age)

# diagnostics

tab.fit <- as.data.frame(
  table(
    x = factor(data$v2, levels = 1:5),
    y = factor(data$v6, levels = 1:5),
    AGE = factor(data$AGE, levels = 0:1)
  )
)

tab.fit$x   <- as.integer(as.character(tab.fit$x))
tab.fit$y   <- as.integer(as.character(tab.fit$y))
tab.fit$AGE <- as.integer(as.character(tab.fit$AGE))

b <- coef(res.biv.age)

tab.fit$p <- mapply(
  function(xi, yi, ai)
    ddBLN.AGE(
      i = xi, j = yi, AGE = ai,
      beta10 = b[["beta10"]],
      beta11 = b[["beta11"]],
      sigma2.1 = b[["sigma2.1"]],
      beta20 = b[["beta20"]],
      beta21 = b[["beta21"]],
      sigma2.2 = b[["sigma2.2"]],
      rho = b[["rho"]],
      k = 5
    ),
  tab.fit$x, tab.fit$y, tab.fit$AGE
)

# Size of the corresponding AGE group
tab.fit$n <- ave(tab.fit$Freq, tab.fit$AGE, FUN = sum)

tab.fit$expected <- tab.fit$n * tab.fit$p
tab.fit$observed.p <- tab.fit$Freq / tab.fit$n

# Numerical check: probabilities must sum up to 1
aggregate(p ~ AGE, data = tab.fit, FUN = sum)


tab.fit$diff <- tab.fit$Freq - tab.fit$expected

tab.fit$resid.pearson <- with(
  tab.fit,
  (Freq - expected) / sqrt(expected)
)

# Cells with largest residuals in absolute value
head(
  tab.fit[
    order(abs(tab.fit$resid.pearson), decreasing = TRUE),
    c("AGE", "x", "y", "Freq", "expected", "diff", "resid.pearson")
  ],
  12
)

# checking the margins across the two groups
aggregate(cbind(Freq, expected) ~ AGE + x, data = tab.fit, sum)
aggregate(cbind(Freq, expected) ~ AGE + y, data = tab.fit, sum)


# joint continuous cdf for the bivariate continuous logit-normal rv
# with Student's t copula
FBLNt <- function(x, y, mu1, sigma2.1, mu2, sigma2.2, rho, nu)
{
  sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  pmvt(lower=c(-Inf,-Inf), upper=c(qt(plogitnorm(x,mu1,sqrt(sigma2.1)),nu),
                                   qt(plogitnorm(y,mu2,sqrt(sigma2.2)),nu)),
       delta=c(0,0), df=nu, corr=sigma)
}
# joint cdf for the bivariate discrete logit-normal rv
FdBLNt <- function(i, j, mu1=0, sigma2.1=1, mu2=0, sigma2.2=1 ,rho=0, nu=Inf, k=5)
{
  i <- floor(i); j<- floor(j)
  ifelse(i>=0 & j>=0,FBLNt(i/k, j/k, mu1, sigma2.1, mu2, sigma2.2, rho, nu), 0)
}
# joint pmf for the bivariate discrete logit-normal rv
ddBLNt <- function(i, j, mu1=0, sigma2.1=1, mu2=0, sigma2.2=1, rho=0, nu=Inf, k=5)
{
  FdBLNt(i, j, mu1, sigma2.1, mu2, sigma2.2, rho, nu, k)+
    FdBLNt(i-1, j-1, mu1, sigma2.1, mu2, sigma2.2, rho, nu, k)-
    FdBLNt(i-1, j, mu1, sigma2.1, mu2, sigma2.2, rho, nu, k)-
    FdBLNt(i, j-1, mu1, sigma2.1, mu2, sigma2.2, rho, nu, k)
}
# log-likelihood function for the bivariate discrete logit-normal
# x.cell and y.cell denote the distinct observed category pairs
# n.cell denotes their corresponding counts
log.lik.dBLN <- function(mu1, sigma2.1, mu2, sigma2.2, rho, nu, k = 5, x.cell, y.cell, n.cell) {
  p <- mapply(
    function(xi, yi)
      ddBLNt(xi, yi, mu1, sigma2.1, mu2, sigma2.2, rho, nu, k),
    x.cell, y.cell
  )
  if (any(!is.finite(p)) || any(p <= 0))
    return(sum(n.cell)^2)
  
  -sum(n.cell * log(p))
}
# data analysis using ISSP2021
data <- dataframe[dataframe$c_alphan=="DK",] 
data <- data[,c("v2","v6")]
data <- data[apply(data >= 0, 1, all), ]
data
x <- data[,1]
y <- data[,2]
tab <- as.data.frame(
  table(x = factor(x, levels = 1:5),
        y = factor(y, levels = 1:5))
)
tab <- subset(tab, Freq > 0)
x.c <- as.integer(as.character(tab$x))
y.c <- as.integer(as.character(tab$y))
n.c <- tab$Freq
# nu = 6
res.biv <- mle2(log.lik.dBLN,start=list(mu1=0, sigma2.1=1, mu2=0, sigma2.2=1, rho=cor(data)[1,2], nu=5),
                fixed=list(k=5,nu=6), data= list(
                  x.cell = x.c,
                  y.cell = y.c,
                  n.cell = n.c
                ),
                method="L-BFGS-B",
                lower=c(mu1=-Inf, sigma2.1=1e-4, mu2=-Inf, sigma2.2=1e-4, rho=-1,nu=0.1),
                upper=c(mu1=Inf, sigma2.1=Inf, mu2=Inf, sigma2.2=Inf, rho=1, nu=Inf),
                control=list(trace=TRUE))
summary(res.biv)
AIC(res.biv)
pij <- matrix(0, 5, 5)
for(h in 1:5)
{
  for(k in 1:5)
  {
    pij[h,k] <- ddBLNt(h, k, res.biv@coef[1], res.biv@coef[2], res.biv@coef[3],
                       res.biv@coef[4], res.biv@coef[5], nu=6, k=5)
  }
}
1/2*sum(abs(pij - fij))



############################################################
###### J O I N T   G R A P H   O F    E   A N D   V  #######
####### for the discrete logit-normal distribution  ########
############################################################
# do not run #
k <- 5
op <- par()
par(mai=c(0.75,0.75,0.1,0.1), mgp=c(2.5,1,0), mfrow=c(1,1))
plot(c(1,k), c(0,4), type="n", xlab="E", ylab="V")
mu <- seq(from=-6, to=6, by=0.025)
sigma <- exp(mu)
t <- length(mu)
for(i in 1:t)
{
  for (j in 1:t)
  {
    points(E.dLN(mu[i], sigma[j], k), V.dLN(mu[i], sigma[j],k), pch=19, col="grey", cex=1) 
  }
}