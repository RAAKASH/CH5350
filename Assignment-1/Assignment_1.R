# Assignment 1

# Question 1

library(rmutil) # Required for int2() function

 pdf_fn <- function(x,y){
   res = (exp(-x/y) * exp(-y))/y
   return(res)
   } # bivariate pdf of random variable x,y

# (i)
K = 1 /(int2(pdf_fn,c(0,0),c(Inf,Inf),eps=10^-8))

# (ii)
marden_y <- function(x,y){ return(exp(-y))}

# (iii)
Pr = (int2(pdf_fn,c(0,0.2),c(1,0.4),eps=10^-8))

# (iv)
cond_exp <- function(x){
  y = 5
  return(x*exp(-x/y)/y) }

Exp = integrate(cond_exp,0,Inf) # Basically E(x|y) = y


# Q2 
# (a) No, Taking expectation we get E(k) = A*cos(2*pi*f0*k)

# (b)
rm(list= ls())
num = 10000
t = 1:num
for (i in 1:num) {
f = 1/20
k = 3
A = 10
l = i
N = 10000
#phi = seq(0.1,100,0.1)
k = seq(l+1,l+N,1)
ch1 = A*cos(2*pi*f*k)+ rnorm(length(k),0,5)
ch2 = A*cos(2*pi*f*(k-l))+ rnorm(length(k),0,5)
t[i] = (sum((ch1-mean(ch1))*(ch2-mean(ch2))))/(N)

}
plot(ch1[7700:8000],type = "l")
plot(1:length(t[7700:8000]),t[7700:8000],type = "l")


# Q3
n = 150000
set.seed(100)
X = rnorm(n,1,2)
Y = 2*X*X+4*X
sum((Y-mean(Y))*(X-mean(X)))/(n-1)
sum((Y-mean(Y))*(Y-mean(Y)))/(n-1)
sum((X-mean(X))*(X-mean(X)))/(n-1)
cov(cbind(X,Y))

  
## Q4
rm(list = ls())
#WN_1 = arima.sim(model =list(order = c(0, 0, 0)) , n = 100, mean = 100, sd =10)
n = c(50,100,200,500,1000,1500,2000,2500,3000,10000) # No of observations
#n = seq(50,10000,500)
r = 100 # No of realizations
#n = 200
#r = 1000
mod_mu = 1:length(n)
mod_var = 1:length(n)
set.seed(100)
Dat = matrix(data=NA,nrow=200,ncol=3)
mu_1   = matrix(data=NA,nrow=length(n),ncol=1)
mu_5   = mu_1
mu_17  = mu_1
var_1  = mu_1
var_5  = mu_1
var_17 = mu_1

for (i in 1:length(n)) {
  Dat    = matrix(data=NA,nrow=r,ncol=3)

   for (j in 1:r) {
     WN  = acf(x=rnorm(n[i]),lag.max = 19)
     Dat[j,1:3] = as.matrix(WN$acf[c(2,6,18)])
   }
  
  Dat = t(Dat) # transpose of matrix
  
  mu_1[i]   = mean(Dat[1,1:ncol(Dat)]) # mean at lag1
  mu_5[i]   = mean(Dat[2,1:ncol(Dat)]) # mean at lag5
  mu_17[i]  = mean(Dat[3,1:ncol(Dat)]) # mean at lag17
  var_1[i]  = var(Dat[1,1:ncol(Dat)])  # var  at lag1
  var_5[i]  = var(Dat[2,1:ncol(Dat)])  # var  at lag5
  var_17[i] = var(Dat[3,1:ncol(Dat)])  # var  at lag17
  
}

# Verifying the 95% Confidence Interval 
gauss <- function(x){
  mu = 1
  sigma = 2
  res = (exp(-((x-mu)^2)/2/sigma^2)/(2*pi*sigma^2)^0.5) 
  return(res)}
integrate(gauss,-2*1.96+1,2*1.96+1)

gauss_new <- function(x,mu_1,var_1){
  mu = mu_1
  sigma = var_1^0.5
  res = (exp(-((x-mu)^2)/2/sigma^2)/(2*pi*sigma^2)^0.5) 
  return(res)}

#************PDF at Lag 1****************#
# Note these statements are to be executed when n=200, r = 1000
hist(Dat[1,1:ncol(Dat)],breaks = 20,col = "red",xlab = "ACF across realizations",main=("PDF of ACF at Lag = 1"),xlim = c(-.2,.2))
lines(seq(-0.2,0.2,0.01),gauss_new(seq(-0.2,0.2,0.01),mu_1,var_1)*120/gauss_new(mu_1,mu_1,var_1))
#************PDF at Lag 5****************#
hist(Dat[2,1:ncol(Dat)],breaks = 20,col = "red",xlab = "ACF across realizations",main=("PDF of ACF at Lag = 5"),xlim =c(-.2,.2))
lines(seq(-0.2,0.2,0.01),gauss_new(seq(-0.2,0.2,0.01),mu_5,var_5)*120/gauss_new(mu_5,mu_5,var_5))
#************PDF at Lag 17****************#
hist(Dat[3,1:ncol(Dat)],breaks = 20,col = "red",xlab = "ACF across realizations",main=("PDF of ACF at Lag = 17"),xlim =c(-.2,.2))
lines(seq(-0.2,0.2,0.01),gauss_new(seq(-0.2,0.2,0.01),mu_17,var_17)*120/gauss_new(mu_17,mu_17,var_17))

# Plotting Confidence Intervals
# Note: Executed when n = in intervals of 50, r = 100

plot(x=1/n^0.5,y = mu_1-1.96*var_1,type = "l",ylim=range( c(min(mu_1-1.96*var_1),max(mu_1+1.96*var_1)) ),col="red",xlab = "1/N^0.5",ylab = "Probability Interval" )
lines(1/n^0.5,mu_1+1.96*var_1,col="blue")
legend(0.02, -0.02, legend=c("Lower Limit", "Upper Limit"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
title("Confidence Interval of 95% - Lag 1")


plot(x=1/n^0.5,y = mu_5-1.96*var_5,type = "l",ylim=range( c(min(mu_5-1.96*var_5),max(mu_5+1.96*var_5)) ),col="red",xlab = "1/N^0.5",ylab = "Probability Interval" )
lines(1/n^0.5,mu_5+1.96*var_5,col="blue")
legend(0.02, -0.02, legend=c("Lower Limit", "Upper Limit"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
title("Confidence Interval of 95% - Lag 5")


plot(x=1/n^0.5,y = mu_17-1.96*var_17,type = "l",ylim=range( c(min(mu_17-1.96*var_17),max(mu_17+1.96*var_17)) ),col="red",xlab = "1/N^0.5",ylab = "Probability Interval" )
lines(1/n^0.5,mu_17+1.96*var_17,col="blue")
legend(0.01, 0.025, legend=c("Lower Limit", "Upper Limit"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
title("Confidence Interval of 95% - Lag 17")

# Z scores
mu_1/var_1^0.5*n^0.5
mu_5/var_5^0.5*n^0.5
mu_17/var_17^0.5*n^0.5

## 3 b
A = 1
f = 0.15
N = 210
l = 10
snr = c(20,10,1,0.1 )

k = seq(l+1,N,1)
set.seed(50)
for (i in 1:length(snr)) {
  

ch1  = A*cos(2*pi*f*k) + ((0.5/snr[i])^0.5) * rnorm(N-l)
ac = acf(ch1,lag.max = 100)
plot(ac$acf,type = "l")
z = fft(ch1)
plot(Mod(z),type = "l")
z[which(Mod(z)==max(Mod(z[1:length(z)/2])))]
print(which(Mod(z)==max(Mod(z[1:length(z)/2]))))
}
((sum((ch1-mean(ch1))*(ch1-mean(ch1)))/N)^0.5)*((sum((ch2-mean(ch2))*(ch2-mean(ch2)))/N)^0.5)
z = fft(ch1)
plot(Mod(z),type = "l")
z[which(Mod(z)==max(Mod(z)))]
which(Mod(z)==max(Mod(z)))