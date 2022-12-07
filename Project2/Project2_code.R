##### Packages #####

library(MASS)
library(spatial)
library(ggplot2)
library(plyr)


# Set seed
set.seed(309)


df_cells <- ppinit("cells.dat")
df_redwood <- ppinit("redwood.dat")
df_pines <- ppinit("pines.dat")


##### -------------- Problem 1 -------------- #####

jpeg("df_cells.jpg")
plot(df_cells, main = 'Plot of cells.dat', xlab = 'x', ylab = 'y')
dev.off()

jpeg("df_redwood.jpg")
plot(df_redwood, main = 'Plot of redwood.dat', xlab = 'x', ylab = 'y')
dev.off()

jpeg("df_pines.jpg")
plot(df_pines, main = 'Plot of pines.dat', xlab = 'x', ylab = 'y')
dev.off()



##### Kfn #####

####--- Cells ----####

ppregion(0,1,0,1)
L_cells = Kfn(df_cells, fs=1)

## Plot
jpeg("L_cells.jpg")
plot(L_cells, type = 'l', lwd = 2, xlab = 'Distance r', ylab = 'L(r)')
abline(0,1, col = 'red')
legend('bottomright', cex = 0.8, c("L(r) = r", "Empirical L-function"), lty=c(1,1), col = c("red", "black"), lwd = c(1,2))
dev.off()

####--- Redwood ----####

# Compute redwood L-function
ppregion(0,1,-1,0)
L_redwood = Kfn(df_redwood,fs=1)

## Plot
jpeg("L_redwood.jpg")
plot(L_redwood, type = 'l', lwd = 2, xlab = 'Distance r', ylab = 'L(r)')
abline(0,1, col = 'red')
legend('bottomright', cex = 0.8, c("L(r) = r", "Empirical L-function"), lty=c(1,1), col = c("red", "black"), lwd = c(1,2))
dev.off()


####--- Pines ----####

# Compute pine L-function
ppregion(0,9.6,0,10)
L_pines = Kfn(df_pines,fs=10)

## Plot
jpeg("L_pines.jpg")
plot(L_pines, type = 'l', lwd = 2, xlab = 'Distance r', ylab = 'L(r)')
abline(0,1, col = 'red')
legend('bottomright', cex = 0.8, c("L(r) = r", "Empirical L-function"), lty=c(1,1), col = c("red", "black"), lwd = c(1,2))
dev.off()




##### 1c #####


# Construct a function to simulate the homogeneous Poisson point process
# conditional on the number of points (n_points)

hom_pois_sim <- function(n_points, x_min = 0, x_max = 1, y_min = 0, y_max = 1){
  ppregion(xl=x_min,xu=x_max,yl=y_min,yu=y_max)
  x_dist = x_max-x_min
  y_dist = y_max-y_min
  x_coord = x_dist*runif(n_points)+x_min  # x-coordinates
  y_coord = y_dist*runif(n_points)+y_min  # y coordinates
  
  simulated = list('x' = x_coord, 'y' = y_coord, 'area' = c(x_min,x_max, y_min, y_max))
  class(simulated) = 'pp'
  return(simulated)
}

# Make an function that calculates the empirical L-function
L_emp <- function(length_n, length_k,fs, xl, xu, yl, yu, sim = 100){
  L_emp_mat = matrix(NA, nrow = sim, ncol = length_k)
  
  # Loop through to find empirical L-function values from the realizations 
  ppregion(xl,xu,yl,yu)
  for(i in 1:sim){
    points = hom_pois_sim(length_n, xl, xu, yl, yu) # length gives number of pts
    
    L_func = Kfn(points, fs, sim) # Construct L
    L_emp_mat[i,] = L_func$y
  }
  
  # Save x-vector to plot
  x_vec_L = L_func$x
  
  # Calculate a 90% CI
  L_emp_upper = apply(L_emp_mat, MARGIN = 2, FUN = quantile, probs=0.95)
  L_emp_lower = apply(L_emp_mat, MARGIN = 2, FUN = quantile, probs=0.05)
  
  return(list(CI_u = L_emp_upper, 
              CI_l = L_emp_lower, 
              x = x_vec_L))
}

# ------------------ Cell point pattern -----------------

L_emp_cells = L_emp(length(df_cells$y),70, fs = 1, 0,1,0,1)


# Plot
jpeg("emp_L_cells.jpg")
plot(L_emp_cells$x, L_cells$y, type = 'l' , lwd = 2, col = "blue2", xlab = 'Distance r', ylab = 'L(r)')
lines(L_emp_cells$x,L_emp_cells$CI_l, type = "l", lty = 2, col = "red", lwd = 1.5)
lines(L_emp_cells$x,L_emp_cells$CI_u, type = "l", lty = 2, col = "red", lwd = 1.5)
legend('bottomright', cex = 0.7, c("Upper/lower 90% prediction interval", "Empirical L-function from the data"), lty=c(2,1), col = c("red", "blue"), lwd = c(1.5,2))
dev.off()


# ------------------ Redwood point pattern -----------------

L_emp_redwood = L_emp(length(df_redwood$y),70,fs = 1, 0,1,-1,0)


# Plot
jpeg("emp_L_redwood.jpg")
plot(L_emp_redwood$x, L_redwood$y, type = 'l' , lwd = 2, col = "blue2", xlab = 'Distance r', ylab = 'L(r)')
lines(L_emp_redwood$x,L_emp_redwood$CI_l, type = "l", lty = 2, col = "red", lwd = 1.5)
lines(L_emp_redwood$x,L_emp_redwood$CI_u, type = "l", lty = 2, col = "red", lwd = 1.5)
legend('bottomright', cex = 0.7, c("Upper/lower 90% prediction interval", "Empirical L-function from the data"), lty=c(2,1), col = c("red", "blue"), lwd = c(1.5,2))
dev.off()

# ------------------ Pine point pattern -----------------

L_emp_pines = L_emp(length(df_pines$y), 70,fs = 10, 0,9.6,0,10, sim = 101)


# Plot
jpeg("emp_L_pines.jpg")
plot(L_emp_pines$x[-70], L_pines$y, type = 'l' , lwd = 2, col = "blue2", xlab = 'Distance r', ylab = 'L(r)')
lines(L_emp_pines$x[-70],L_emp_pines$CI_l[-70], type = "l", lty = 2, col = "red", lwd = 1.5)
lines(L_emp_pines$x[-70],L_emp_pines$CI_u[-70], type = "l", lty = 2, col = "red", lwd = 1.5)
legend('bottomright', cex = 0.7, c("Upper/lower 90% prediction interval", "Empirical L-function from the data"), lty=c(2,1), col = c("red", "blue"), lwd = c(1.5,2))
dev.off()




##### -------------- Problem 2 -------------- #####


library(fields)
library(akima)

# plot the data

# read in the files
df_obs <- read.table('obspines.txt', header = T)
df_prob <- read.table('obsprob.txt', header = T)

x = unique(df_obs$x)
y = unique(df_obs$y)
z = matrix(df_obs$N_obs, nrow = length(x), ncol = length(y))

jpeg("obs_plot.jpg", height = 550, width = 550)
par(cex = 1.5, mar=c(3,3,2,2),mgp = c(2, 1, 0))

image.plot(x,y,z, 
           col = viridis(4),
           legend.args = list( text = "M",
                              cex = 1.5,
                              side = 3,
                              line = 0.6))
dev.off()


x = unique(df_prob$x)
y = unique(df_prob$y)
z = matrix(df_prob$alpha, nrow = length(x), ncol = length(y))


jpeg("prob_plot.jpg", height = 550, width = 550)
par(cex = 1.5, mar=c(3,3,2,2),mgp = c(2, 1, 0))

image.plot(x,y,z, 
           col = viridis(63),
           legend.args = list( text = expression(alpha),
                               cex = 1.5,
                               side = 3,
                               line = 0.6))
dev.off()



### estimate lambda_hat ###
lambda_hat = sum(df_obs$N_obs/df_prob$alpha)/300^2

# We simulate for the full grid, since it is the same when we have 
# a homogeneous Poisson process

jpeg("2csim1.jpg")
par(mar = c(3,3,3,4), mgp = c(2,1,0))
sim1 = hom_pois_sim(rpois(1,lambda_hat*10^2*900), 0,300,0,300)
plot(sim1, pch=19, cex = 0.4, xlab = 'x', ylab = 'y', col = '#3949AB')
dev.off()

jpeg("2csim2.jpg")
par(mar = c(3,3,3,4), mgp = c(2,1,0))
sim1 = hom_pois_sim(rpois(1,lambda_hat*10^2*900), 0,300,0,300)
plot(sim1, pch=19, cex = 0.4, xlab = 'x', ylab = 'y', col = '#3949AB')
dev.off()

jpeg("2csim3.jpg")
par(mar = c(3,3,3,4), mgp = c(2,1,0))
sim1 = hom_pois_sim(rpois(1,lambda_hat*10^2*900), 0,300,0,300)
plot(sim1, pch=19, cex = 0.4, xlab = 'x', ylab = 'y', col = '#3949AB')
dev.off()


### 2d ###

lambda_int = matrix((1-df_prob$alpha)*lambda_hat, nrow = 30)


sim_cond <- function(){
  sim = list(x = NULL, y = NULL)
  
  for (i in (1:30)){
    for (j in (1:30)){
      grid_sim = hom_pois_sim(rpois(1,lambda_int[i,j]*10^2),(i-1)*10, i*10,(j-1)*10, j*10)
      sim$x = append(sim$x,grid_sim$x)
      sim$y = append(sim$y, grid_sim$y)
    }
  }
  
  return(sim)
}


sim_cond1 = sim_cond()

sim_cond2 = sim_cond()

sim_cond3 = sim_cond()


jpeg("sim_cond1.jpg")
par(mar = c(3,3,3,4), mgp = c(2,1,0))
plot(sim_cond1, pch=19, cex = 0.4, xlab = 'x', ylab = 'y', col = '#E65100')
dev.off()

jpeg("sim_cond2.jpg")
par(mar = c(3,3,3,4), mgp = c(2,1,0))
plot(sim_cond2, pch=19, cex = 0.4, xlab = 'x', ylab = 'y', col = '#E65100')
dev.off()

jpeg("sim_cond3.jpg")
par(mar = c(3,3,3,4), mgp = c(2,1,0))
plot(sim_cond3, pch=19, cex = 0.4, xlab = 'x', ylab = 'y', col = '#E65100')
dev.off()




# --------------- 500 realizations -------------

sim_500 = matrix(NA, nrow = 30, ncol = 30)
sim_500_sd = matrix(NA, nrow = 30, ncol = 30)


for (i in (1:30)){
  for (j in (1:30)){
    p_sim = rpois(500, lambda = lambda_hat*10^2)
    sim_500[i,j] = mean(p_sim)
    sim_500_sd[i,j] = sd(p_sim)
  }
}

  
sim_cond_500 = matrix(NA, nrow = 30, ncol = 30)
sim_cond_500_sd = matrix(NA, nrow = 30, ncol = 30)

for (i in (1:30)){
  for (j in (1:30)){
   p_sim = rpois(500, lambda = lambda_int[i,j]*10^2)
   sim_cond_500[i,j] = mean(p_sim) + matrix(df_obs$N_obs,nrow = 30, ncol = 30)[i,j]
   sim_cond_500_sd[i,j] = sd(p_sim)
  }
}


x = unique(df_obs$x)
y = unique(df_obs$y)
z = sim_500

jpeg("priori_plot.jpg", height = 550, width = 550)
par(cex = 1.5, mar=c(4.5,4.5,3.5,3.5),mgp = c(2, 1, 0))

image.plot(x,y,z, 
           col = viridis(8),
           legend.args = list( text = "Expected obs",
                               cex = 1.3,
                               side = 3,
                               line = 0.6))
dev.off()

z = sim_500_sd

jpeg("priori_sd.jpg", height = 550, width = 550)
par(cex = 1.5, mar=c(4.5,4.5,3.5,3.5),mgp = c(2, 1, 0))

image.plot(x,y,z, 
           col = inferno(8),
           legend.args = list( text = "std",
                               cex = 1.3,
                               side = 3,
                               line = 0.6))
dev.off()

z = sim_cond_500

jpeg("posterior_plot.jpg", height = 550, width = 550)
par(cex = 1.5,mar=c(4.5,4.5,3.5,3.5),mgp = c(2,1,0))

image.plot(x,y,z, 
           col = viridis(8),
           legend.args = list( text = "Expected obs",
                               cex = 1.3,
                               side = 3,
                               line = 0.6))
dev.off()

z = sim_cond_500_sd

jpeg("posterior_sd.jpg", height = 550, width = 550)
par(cex = 1.5,mar=c(4.5,4.5,3.5,3.5),mgp = c(2,1,0))

image.plot(x,y,z, 
           col = inferno(8),
           legend.args = list( text = "std",
                               cex = 1.3,
                               side = 3,
                               line = 0.6))
dev.off()



##### -------------- Problem 3 -------------- #####


# We first plot the data again to try to guess parameters.

plot(df_redwood, main = 'Plot of redwood.dat', xlab = 'x', ylab = 'y')

sim_NS <- function(lambda_p, lambda_d, sigma, xl = 0, xu = 1, yl = -1, yu = 0){
  # Extend window before simulating
  ext_window_x = c(xl - 2*sigma,xu + 2*sigma)
  ext_window_y = c(yl - 2*sigma,yu + 2*sigma)
  ext_window_area = diff(ext_window_x)*diff(ext_window_y)
  
  parents = hom_pois_sim(rpois(1,lambda_p*ext_window_area),xl,xu,yl,yu)
  
  x_p = parents$x
  y_p = parents$y
  
  cov_mat = sigma^2*diag(1,2)
  
  x_d = 1
  y_d = 1
  for (i in (1:length(x_p))){
    n_d = rpois(1, lambda_d)
    if(n_d != 0){
      mu = c(x_p[i],y_p[i])
      d_points = mvrnorm(n_d,mu,cov_mat)
      
      if(n_d == 1){
        x_d = append(x_d, d_points[1])
        y_d = append(y_d, d_points[2])
      }
      
      else{
        x_d = append(x_d, d_points[,1])
        y_d = append(y_d, d_points[,2])
      }
      
    }
    
  }
  
  x_d[-1]
  y_d[-1]
  
  remove = which(x_d>xu | x_d<xl | y_d>yu|y_d< yl)
  x_d = x_d[-remove]
  y_d = y_d[-remove]
  
  return(d_points = list(x=x_d, y=y_d))
}



## ------------------ Guestimate 1 -----------------

# Since around 40% of the points fall into a circle of radius sigma, 
# we should have the "radius" such that the circle around "the center of 
# the cluster contains 0.4*lambda_c points as sigma. This 0.4*lambda_c is about 2.75.
# Guessing 0.05 seems fair, based on how I have "selected"/guessed my clusters.

sigma = 0.05

# guess parents 
lambda_p = 9

# guess daughters
guessed_daughter_count = sum(c(5,9,8,13,6,7,6,2,6))
lambda_d = guessed_daughter_count/lambda_p

jpeg("first_prob3.jpg")
plot(sim_NS(lambda_p,lambda_d,sigma), xlab = 'x', ylab= 'y', main = "A realization from the first guestimated Neyman-Scott process")
dev.off()


## ------------------ Guestimate 2 -----------------

# 

sigma = 0.05

# guess parents 
lambda_p = 11

# guess daughters
lambda_d = 6

jpeg("second_prob3.jpg")
plot(sim_NS(lambda_p,lambda_d,sigma), xlab = 'x', ylab= 'y',  main = "A realization from the last guestimated Neyman-Scott process")
dev.off()

# ---------------- Construct empirical L-functions -----------

# Construct a matrix to store the L-function values for each realization
L_3_sims = matrix(NA, nrow = 100, ncol = 70)
L_3_sims_upper = rep(NA,70)
L_3_sims_lower = rep(NA,70)

# Iterate through 100 realizations
for (i in (1:100)){
  sim_L = sim_NS(lambda_p, lambda_d, sigma)
  L_NS_temp = Kfn(sim_L, 1, 100)
  L_3_sims[i,] = L_NS_temp$y
}

# Find upper and lower prediction interval
for(i in 1:70){
  L_3_sims_upper[i] = quantile(L_3_sims[,i], 0.95)
  L_3_sims_lower[i] = quantile(L_3_sims[,i], 0.05)
}

#jpeg("3firstguess.jpg",width = 350, height = 350)

plot(L_redwood$x, L_redwood$y, col = "blue", type = "l", lwd = 2, xlab = "r", ylab = "L(r)", main = "Plot of empirical L-function with empirical prediction interval from guestimated model.")
lines(L_redwood$x, L_3_sims_upper, type = "l", col = "darkgreen", lwd = 2, lty = 2)
lines(L_redwood$x, L_3_sims_lower, type = "l", col = "darkgreen", lwd = 2, lty = 2)
legend('topleft', cex = 0.7, c("Upper/lower empirical 90% prediction interval from guestimated model", "Empirical L-function from the data"), lty=c(2,1), col = c("darkgreen", "blue"), lwd = c(1.5,2))

#dev.off()










##### -------------- Problem 4 -------------- #####

# -------- Guestimating parameters -------------

# We plot the data to try to guess.
plot(df_cells, main = 'Plot of cells.dat', xlab = 'x', ylab = 'y')

# We know the true number of points.
num_strauss = length(df_cells$x) # = 42

### GUESTIMATE 1
# It looks like points at least stay outside of balls of radius 0.1 from each other, so we guestimate r_0=0.1.

#r_strauss = 0.1

# If we assume r_0 = r_strauss = 0.1, we see that few points violate this. We guess $beta = 1$.

#beta_strauss = 1
#c_strauss = exp(-beta_strauss)


# ### GUESTIMATE 2
# r_strauss = 0.08
# beta_strauss = 3
# c_strauss = exp(-beta_strauss)


### GUESTIMATE 3
r_strauss = 0.08
beta_strauss = 10
c_strauss = exp(-beta_strauss)

# -------- Construct realizations -------------

# Make an empty matrix to fill with values for L-func and pred.ints.
L_emp_strauss = matrix(NA, nrow = 100, ncol= 70)
L_strauss_upper = rep(NA,70)
L_strauss_lower = rep(NA,70)

ppregion()

x= Strauss(9,0.5,3.5)

plot(x)
# Simulate Strauss to get L-function.
for (i in (1:100)){
  strauss_sim = Strauss(n = num_strauss, c = c_strauss, r = r_strauss)
  L_strauss_temp = Kfn(strauss_sim, 1, 100)
  L_emp_strauss[i,] = L_strauss_temp$y
}

for (j in (1:70)){
  L_strauss_upper[j] = quantile(L_emp_strauss[,j], 0.95)
  L_strauss_lower[j] = quantile(L_emp_strauss[,j], 0.05)
}

#jpeg("4firstguess_strauss.jpg",width = 450, height = 450)
plot(L_cells$x, L_cells$y, col = "blue", lwd = 1, ylim = c(0,0.7))
lines(L_cells$x, L_strauss_upper, type = "l", col = "darkgreen", lwd = 1)
lines(L_cells$x, L_strauss_lower, type = "l", col = "darkgreen", lwd = 1)
#dev.off()

# --------- Generate three realizations from the guestimated model -----------

ppregion()
par(mfrow = c(2,2), cex.main = 1.5)

# Plot cell data
plot(df_cells, main = 'Plot of cells.dat', xlab = 'x', ylab = 'y', xlim=c(0,1), ylim = c(0,1))

# Make and plot three realizations from the last parameters
for (i in (1:3)){
  strauss_sim = Strauss(n = num_strauss, c = c_strauss, r = r_strauss)
  plot(strauss_sim, main = "Plot of realization from guestimated Strauss model", xlab = 'x', ylab = 'y', xlim=c(0,1), ylim = c(0,1))
}









