#### -------------- Packages -------------- #### 

library(spatial)
library(ggplot2)
library(rgdal)
library(spdep)
library(Matrix)
library(MASS)
library(boot) #inv.logit function

set.seed(232)
setwd("C:/Users/noraa/OneDrive/Documents/BMAT/TMA4250/Project3")

#### -------------- load data -------------- #### 

load("Admin1Geography.RData")
load("Admin2Geography.RData")
Admin1Graph <- read.csv("Admin1Graph.txt", sep="")
Admin2Graph <- read.csv("Admin2Graph.txt", sep="")
DirectEstimates <- read.csv("DirectEstimates.txt", sep="")

# The code from functions.R

## plotAreaCol
# This functions plots the values and saves the figure to the specified file.
# It is advisable to not plot directly in R since it will take a long time.
# The arguments are:
#   fNamme: file name for saving figure
#   width: width of figure in inches
#   height: height of figure in inches
#   estVal: the k values to be plotted on geoMap
#   geoMap: the map containing k regions
#   leg: name to use on top of legend
#   colLim: control lower limit and upper limit of color scale (colLim = c(lowVal, highVal))
plotAreaCol = function(fName, estVal, geoMap, leg, width = 11, height = 11, colLim = NULL){
  if(is.null(colLim)){
    colLim = range(estVal)
  }
  
  # Set up data object for plotting
  nigeriaMapTmp = geoMap
  nigeriaMapTmp$MCV1 = estVal
  nigeria.df = merge(fortify(nigeriaMapTmp), as.data.frame(nigeriaMapTmp), by.x = "id", by.y = 0)
  nigeria.df$Longitude = nigeria.df$long
  nigeria.df$Latitude  = nigeria.df$lat
  
  # Plot
  map = ggplot() +
    geom_polygon(data = nigeria.df,
                 aes(x = Longitude, y = Latitude, group = group, fill = MCV1),
                 color = 'gray', size = .2)+
    scale_fill_viridis_c(direction = 1,
                         begin = 1,
                         end = 0,
                         limit = colLim,
                         name = leg) + 
    coord_fixed() + 
    theme(text = element_text(size=40),
          legend.key.height = unit(4, 'cm'),
          legend.key.width  = unit(1.75, 'cm'))
  ggsave(filename = fName,
         plot = map,
         width = width, 
         height = height)
}


#### -------------- Problem 1a -------------- #### 

# Make a function that constructs the 
# precision matrix up to scaling
construct_prec <- function(M){
  n_neig = colSums(M)
  R = -M
  diag(R) = n_neig
  return(R)
}

R_1 = construct_prec(Admin1Graph)
R_2 = construct_prec(Admin2Graph)


# Check dim and rank
dim(R_1)

dim(R_2)

qr(R_1)$rank

qr(R_2)$rank

# find proportion of non-zero elements
length(which(R_1 != 0))/length(R_1)^2
length(which(R_2 != 0))/length(R_2)^2

jpeg("sparse_admin1.jpg", height = 550, width = 550)
par(cex = 1.5,mar=c(3,3,3,4),mgp = c(2,1,0))

x = rev(seq(1,length(R_1)))
y = x
z = as.matrix(R_1)
z[which(z != 0)] = 1
image(x,y,z,
      xlab = 'row', ylab = 'col',
      col = c('white','red'))
dev.off()


jpeg("sparse_admin2.jpg", height = 550, width = 550)
par(cex = 1.5,mar=c(3,3,3,4),mgp = c(2,1,0))

x = rev(seq(1,length(R_2)))
y = x
z = as.matrix(R_2)
z[which(z != 0)] = 1
image(x,y,z,
      xlab = 'row', ylab = 'col',
      col = c('white','red'))

dev.off()


#### -------------- Problem 1b -------------- #### 

sim_firstorder_GMRF <- function(Q, eps = 1e-10){
  l = length(Q)
  Q_star = Q + eps*diag(l)
  L = chol(Q_star)
  z = rnorm(l)
  v = solve(L)%*%z
  x = v - mean(v)*rep(1,l)
  return(x)
}

sim_admin1 = sim_firstorder_GMRF(R_1)

plotAreaCol("1bsimadm1plot1.jpg", sim_admin1, nigeriaAdm1, "values", colLim = c(-3,3))

sim_norm1 = rnorm(length(R_1))

plotAreaCol("1bnormadm1plot1.jpg", sim_norm1, nigeriaAdm1, "values", colLim = c(-3,3))

sim_admin1 = sim_firstorder_GMRF(R_1)

plotAreaCol("1bsimadm1plot2.jpg", sim_admin1, nigeriaAdm1, "values", colLim = c(-3,3))

sim_norm1 = rnorm(length(R_1))

plotAreaCol("1bnormadm1plot2.jpg", sim_norm1, nigeriaAdm1, "values", colLim = c(-3,3))



#### -------------- Problem 1c -------------- #### 

sim_admin2 = sim_firstorder_GMRF(R_2)

plotAreaCol("1csimadm2plot1.jpg", sim_admin2, nigeriaAdm2, "values", colLim = c(-4,4))

sim_norm2 = rnorm(length(R_2))

plotAreaCol("1cnormadm2plot1.jpg", sim_norm2, nigeriaAdm2, "values", colLim = c(-4,4))

sim_admin2 = sim_firstorder_GMRF(R_2)

plotAreaCol("1csimadm2plot2.jpg", sim_admin2, nigeriaAdm2, "values", colLim = c(-4,4))

sim_norm2 = rnorm(length(R_2))

plotAreaCol("1cnormadm2plot2.jpg", sim_norm2, nigeriaAdm2, "values", colLim = c(-4,4))


#### -------------- Problem 1d -------------- #### 

sim_admin2_100 = matrix(NA,nrow = 775, ncol = 100)

for (i in (1:100)){
  sim_admin2_100[,i] = sim_firstorder_GMRF(R_2)
}

marg_var_100 = apply(sim_admin2_100,1,var)

plotAreaCol("1dvaradm2.jpg", marg_var_100, nigeriaAdm2, "var")

cov_vec = rep(NA, 775)
for (i in (1:775)){
  cov_vec[i] = cov(sim_admin2_100[150,],sim_admin2_100[i,])
}

cor_vec = rep(NA, 775)
for (i in (1:775)){
  cor_vec[i] = cor(sim_admin2_100[150,],sim_admin2_100[i,])
}

plotAreaCol("1dcovadm2.jpg", cov_vec, nigeriaAdm2, "cov", colLim = range(cor_vec))

plotAreaCol("1dcoradm2.jpg", cor_vec, nigeriaAdm2, "cor", colLim = range(cor_vec))




#### -------------- Problem 2a -------------- #### 

plotAreaCol("2aobsprob.jpg", inv.logit(DirectEstimates$Observation),nigeriaAdm1,"p", colLim = c(0,1))

#### -------------- Problem 2b -------------- #### 

# store observations in vector y
y = DirectEstimates$Observation

# store std_dev in vector V and make cov matrix
V = DirectEstimates$StdDev^2

V_covmat = diag(V)

# def parameters for x|y
sigma = 100

cond_mu = solve((1/sigma^2)*diag(37) + diag(1/V))%*%(diag(1/V)%*%y)

cond_cov = solve((1/sigma^2)*diag(37) + diag(1/V))

X_samples = mvrnorm(100, mu = cond_mu, Sigma = cond_cov)

P_samples = inv.logit(X_samples)


# compute our measures of interest
emp_median_P = apply(P_samples,2,median)

emp_mu_P = apply(P_samples,2,mean)

emp_sigma_P = apply(P_samples,2,sd)

emp_CV_P = emp_sigma_P/emp_mu_P


# plot
plotAreaCol("2bmedplot.jpg", emp_median_P, nigeriaAdm1, "median", colLim = c(0,1))

plotAreaCol("2bCVplot.jpg", emp_CV_P, nigeriaAdm1, "CV", colLim = c(0,0.4))

#### -------------- Problem 2c -------------- #### 

make_parm <- function(tau, V, R_1, y){
  Q = diag(1/V)
  Q_new = (tau*R_1 + Q)
  
  mu_new = solve(Q_new)%*%Q%*%y
  return(list(mu = mu_new, Q = Q_new))
}

new_parm = make_parm(1,V, R_1, y)


# Check rank

qr(new_parm$Q)$rank


# repeat 2b
X_samples_new = mvrnorm(100, mu = new_parm$mu, Sigma = solve(new_parm$Q))

P_samples_new = inv.logit(X_samples_new)


# compute our measures of interest
emp_median_P2 = apply(P_samples_new,2,median)

emp_mu_P2 = apply(P_samples_new,2,mean)

emp_sigma_P2 = apply(P_samples_new,2,sd)

emp_CV_P2 = emp_sigma_P2/emp_mu_P2


# plot
plotAreaCol("2cmedplot.jpg", emp_median_P2, nigeriaAdm1, "median", colLim = c(0,1))

plotAreaCol("2cCVplot.jpg", emp_CV_P2, nigeriaAdm1, "CV", colLim = c(0,0.4))



#### -------------- Problem 2d -------------- #### 

y_38 = c(y,0.5)

make_parm2 <- function(tau, V, R_1, y, A){
  Q = diag(1/V)
  Q_new = (tau*R_1 + t(A)%*%Q%*%A)
  
  mu_new = solve(Q_new)%*%t(A)%*%Q%*%y
  return(list(mu = mu_new, Q = Q_new))
}

which(colnames(Admin1Graph)=="Kaduna")

A = diag(rep(1,37))

A = rbind(A,c(rep(0,18),1,rep(0,18)))

parm_38 = make_parm2(1,c(V,0.1^2), R_1, y_38, A)

# repeat 2b
X_samples_38 = mvrnorm(100, mu = parm_38$mu, Sigma = solve(parm_38$Q))

P_samples_38 = inv.logit(X_samples_38)


# compute our measures of interest
emp_median_P38 = apply(P_samples_38,2,median)

emp_mu_P38 = apply(P_samples_38,2,mean)

emp_sigma_P38 = apply(P_samples_38,2,sd)

emp_CV_P38 = emp_sigma_P38/emp_mu_P38


# plot
plotAreaCol("2dmedplot.jpg", emp_median_P38, nigeriaAdm1, "median", colLim = c(0,1))

plotAreaCol("2dCVplot.jpg", emp_CV_P38, nigeriaAdm1, "CV", colLim = c(0,0.4))


#### -------------- Problem 2e -------------- #### 

parm_01 = make_parm(0.1, V, R_1, y)
parm_10 = make_parm(10, V, R_1, y)

# repeat 2c

X_01 = mvrnorm(100, mu = parm_01$mu, Sigma = solve(parm_01$Q))
X_10 = mvrnorm(100, mu = parm_10$mu, Sigma = solve(parm_10$Q))


P_01= inv.logit(X_01)
P_10= inv.logit(X_10)

# compute our measures of interest
median_P01 = apply(P_01,2,median)
median_P10 = apply(P_10,2,median)

mu_P01 = apply(P_01,2,mean)
mu_P10 = apply(P_10,2,mean)

sigma_P01 = apply(P_01,2,sd)
sigma_P10 = apply(P_10,2,sd)

CV_P01 = sigma_P01/mu_P01
CV_P10 = sigma_P10/mu_P10


# plot
plotAreaCol("2emedplot01.jpg", median_P01, nigeriaAdm1, "median", colLim = c(0,1))

plotAreaCol("2eCVplot01.jpg", CV_P01, nigeriaAdm1, "CV", colLim = c(0,0.4))

plotAreaCol("2emedplot10.jpg", median_P10, nigeriaAdm1, "median", colLim = c(0,1))

plotAreaCol("2eCVplot10.jpg", CV_P10, nigeriaAdm1, "CV", colLim = c(0,0.4))


#### -------------- Problem 2f -------------- #### 
# take one random X from a previous simulation
x = X_10[35,]
R_1 = as.matrix(R_1)
Q = diag(1/V)
y = DirectEstimates$Observation

# define our likelihood function
f <- function(tau){
  return(18*log(tau) - (tau/2) * (t(x)%*%R_1%*%x) - (1/2) * (t(y-x)%*%Q%*%(y-x)) - (1/2) * (log(det(tau*R_1 + Q))) + (1/2) * ((t(x-solve(tau*R_1 + Q)%*%Q%*%y))%*%(tau*R_1 + Q)%*%(x-solve(tau*R_1 + Q)%*%Q%*%y)))
}

tau_MLE = optimize(f, interval = c(0,10000), maximum = T)$maximum

# repeat 2c
parm_MLE = make_parm(tau_MLE, V, R_1, y)

X_MLE = mvrnorm(100, mu = parm_MLE$mu, Sigma = solve(parm_MLE$Q))

P_MLE= inv.logit(X_MLE)

# compute our measures of interest
median_PMLE = apply(P_MLE,2,median)

mu_PMLE = apply(P_MLE,2,mean)

sigma_PMLE = apply(P_MLE,2,sd)

CV_PMLE = sigma_PMLE/mu_PMLE


# plot
plotAreaCol("2fmedplotMLE.jpg", median_PMLE, nigeriaAdm1, "median", colLim = c(0,1))

plotAreaCol("2fCVplotMLE.jpg", CV_PMLE, nigeriaAdm1, "CV", colLim = c(0,0.4))
