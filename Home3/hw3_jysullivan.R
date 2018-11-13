# Assignment: write function that take information of natural mortality,
# fecundity-at-age, weight-at-age, selectivity, R0, steepness, the form of the
# stock-recruitment relationships, and gamma
library(knitr)
library(ggplot2)
library(dplyr)
library(kableExtra)
library(data.table)
# library(cowplot)

# Stupid function to add a main title to base r paneled plots (https://stackoverflow.com/questions/14660372/common-main-title-of-a-figure-panel-compiled-with-parmfrow)
line2user <- function(line, side) {
  lh <- par('cin')[2] * par('cex') * par('lheight')
  x_off <- diff(grconvertX(0:1, 'inches', 'user'))
  y_off <- diff(grconvertY(0:1, 'inches', 'user'))
  switch(side,
         `1` = par('usr')[3] - line * y_off * lh,
         `2` = par('usr')[1] - line * x_off * lh,
         `3` = par('usr')[4] + line * y_off * lh,
         `4` = par('usr')[2] + line * x_off * lh,
         stop("side must be 1, 2, 3, or 4", call.=FALSE))
}

# Function to create Z matrix
get_Z <-
  function(age, # vector of ages
           M = 0.15, # natural mortality
           sel, # vector selectivity-at-age (sex-specific)
           F_vals) { # vector of fishing mortality values
    
    Z <- matrix(ncol = length(unique(age)), 
                nrow = length(F_vals))
    
    for(f in 1:length(F_vals)){
      for(a in 1:length(unique(age))) {
        Z[f,a] <- M + F_vals[f] * sel[a] }}
    return(Z)
  }

# Function to get numbers-at-age matrix: N is numbers-at-age relative to the
# number of animals in age-0
get_N <- 
  function(age, # vector of ages
           Z) { # matrix of total mortality-at-age
    
    # N matrix of same dimensions as Z matrix
    N <- matrix(ncol = length(unique(age)), 
                nrow = nrow(Z))       
    a_plus <- max(age) # plus group
    
    # Loop over all input values of mortality
    for(f in 1:length(Z[,1])) {
      
      # Initialize N matrix
      N[f,1] <- 0.5  
      
      for(a in 2:a_plus) {
        N[f,a] <- N[f,a-1] * exp(-Z[f,a-1]) }
      
      # Plus group
      N[f,a_plus+1] <- (N[f,a_plus] * exp(-Z[f,a_plus])) / (1 - exp(-Z[f,a_plus+1]))
    }
    return(N)
  }

# Function to get yield-per-recruit
get_ypr <- 
  function(age, # vector of ages
           wt, # vector of weight-at-age
           M = 0.15, # instantaneous natural mortality
           Z, # matrix of instantaneous total mortality-at-age
           N) { # matrix of numbers-at-age
    
    # Y matrix of same dimensions as Z matrix
    Y <- matrix(ncol = length(unique(age)), 
                nrow = nrow(Z)) 
    
    # Loop over all input values of mortality
    for(f in 1:length(Z[,1])) {
      for(a in 1:length(Z[1,])) {
        Y[f,a] <- wt[a] * (Z[f,a]-M) / Z[f,a] * N[f,a] * (1 - exp(-Z[f,a])) }}
    
    # Sum across ages for each F value
    Y <- rowSums(Y) 
    return(Y)
  }

# Function to get spawner-per-recruit
get_spr <-
  function(age, # vector of ages
           fec, # vector of fecundity-at-age
           N) { # matrix of numbers-at-age
    
    # spr matrix of same dimensions as N matrix
    spr <- matrix(ncol = length(unique(age)), 
                  nrow = nrow(N))
    
    # Loop over all input values of fishing mortality
    for(f in 1:length(N[,1])) {
      for(a in 1:length(unique(age))){
        spr[f,a] <- fec[a] * N[f,a] }}
    
    # Sum across ages for each F value
    spr <- rowSums(spr) 
    return(spr)
  }

# Function to derive stock-recruitment parameters 
get_sr_params <- 
  function(rec_mod, # recruitment model
           spr0, # spawner biomass-per-recruit (unfished)
           h, # steepness
           r0, # unfished recruitment
           gamma) {# pella-tom parameter
    
    # Beverton-Holt model
    if(rec_mod == 1) {
      # alpha <- (5*h/0.8)/spr0
      # beta <- log(5*h)/(0.8*spr0*r0) }
      alpha <- spr0*(1-h)/(4*h)
      beta <- (5*h-1)/(4*h*r0) }
    
    # Ricker model
    if(rec_mod == 2) {
      alpha <- exp(log(5*h)/0.8)/spr0
      beta <- log(5*h)/(0.8*spr0*r0) }
    
    # Pella-Tomlinson model
    if(rec_mod == 3) {
      alpha <- 1/spr0
      beta <- (5*h-1)/(1-0.2^gamma) }
    
    recruit_params <- list(alpha, beta)
    return(recruit_params)
  }

# Function to get recruitment conditioned on alpha/beta and spr
get_recruit <-
  function(rec_mod, # recruitment model
           spr, # spawner biomass-per-recruit (fished)
           s0, # unfished spawning biomass
           alpha, # sr parameter
           beta, # sr parameter
           gamma) { # pella-tomlinson 3rd parameter
    # Beverton-Holt model (Eqn 9)
    if(rec_mod == 1) {
      r <- (spr-alpha)/(beta*spr) }
    
    # Ricker model (Eqn 10a, put in terms of spr)
    if(rec_mod == 2) {
      r <- log(alpha*spr)/(beta*spr) }
    
    # Pella-Tomlinson model (Eqn 10b, put in terms of spr)
    if(rec_mod == 3) {
      r <-  (s0/spr) * (1 - (1-alpha*spr)/(alpha*beta*spr))^(1/gamma)}
    
    return(r)
  }

# Inputs: natural mortality, fecundity-at-age, weight-at-age, 
# and selectivity

dat <- read.table("HOME3/HOME3.txt", header = TRUE)
colnames(dat) <- c("age", "fec", "wt_f", "wt_m", "sel_f", "sel_m")

FF <- 0.01 #seq(0, 1, 0.001) # vector of fishing mortality
h <- 0.5 # steepness parameter
r0 <- 1 # unfished recruitment
gamma <- 1 # third parameter for pella-tomlinson model
MOD <- 1
sprf0 <- 1
s0 <- 1
M <- 0.15

yield_f <- function(FF, h, MOD) {
  
  # Get matrices of sex-specific total instantaneous mortalities 
  Zf <- get_Z(age = dat$age, sel = dat$sel_f, F_vals = FF)
  Zm <- get_Z(age = dat$age, sel = dat$sel_m, F_vals = FF)
  
  # Get sex-specific numbers-at-age matrices
  Nf <- get_N(Z = Zf, age = dat$age)
  Nm <- get_N(Z = Zm, age = dat$age)
  
  # Get yield vectors different values of F
  ypr_f <- get_ypr(age = dat$age, wt = dat$wt_f, Z = Zf, N = Nf)   
  ypr_m <- get_ypr(age = dat$age, wt = dat$wt_m, Z = Zm, N = Nm)   
  ypr <- ypr_f + ypr_m
  
  # Get spawner biomass-per-recruit across different values of F
  spr <- get_spr(age = dat$age, fec = dat$fec, N = Nf)
  ### AEP - spr0 should be global
  
  # Get stock-recruitment parameters
  ab <- get_sr_params(rec_mod = MOD, spr0 = sprf0, h = h, r0 = r0, gamma = gamma)
  
  # Get recruitment values
  r <- get_recruit(rec_mod = MOD, spr = spr, s0 = s0, alpha = ab[[1]], beta = ab[[2]], gamma = gamma)
  
  out <- NULL
  out$s0 <- s0
  out$sprf0 <- sprf0
  out$r <- r
  out$spr <- spr # Yield
  out$s <- (spr*r)
  out$s_ratio <- (spr*r)/s0 # Express S(F) as a proportion of S(0)
  out$y <- ypr*r
  out$ypr <- ypr
  return(out)
}

sprf0 <- yield_f(FF = 0, h = h, MOD = MOD)$spr 
s0 <-  r0*sprf0 # unfished ssb 

FF <- seq(0, 1, 0.0001)

# Function to get derivative of Y(F) with respect to F (use uniroot() to solve for dY(F)/df=0)
dfx.dx <- function(FF, yield_f, delta, h, MOD) {
  y1 <- yield_f(FF-delta/2, h = h, MOD = MOD)$y
  y2 <- yield_f(FF+delta/2, h = h, MOD = MOD)$y
  approx.gradient <- (y2-y1)/delta
  return(approx.gradient)
}

# Get F crash (where F is minimized and S ~ 0 using the bisection method - MK helped me with this.
bisect <- function(Fmin = 0, Fmax = 1){
  for(b in 1:10000){
    F_tst <- (Fmin + Fmax)/2 # update
    s_tmp <- yield_f(FF = F_tst, h = h, MOD = MOD)$s
    if(round(s_tmp,4) == 0 & (Fmax - Fmin) > 0.0002){ return(F_tst)
    } else if(round(s_tmp,4) > 0) { Fmin <- F_tst 
    } else if(round(s_tmp,4) < 0) { Fmax <- F_tst }
  }
  print('max iter')
}

MOD_NAMES <- c("Beverton-Holt", "Ricker", "Pella-Tomlinson")

# Get biological reference points (brp's)
brps <- data.frame(Model = NA, Fmsy = NA, MSY = NA, Fcrash = NA)

# Loop over the 3 S-R models
for(MOD in 1:3){
  # uniroot will only take one value at a time, it will pick values on the interval and feed them as FF's to dfx.dx
  Fmsy <- uniroot(f = dfx.dx, interval = c(0.01,1), tol = 0.000001, 
                  yield_f = yield_f, delta = 0.001, h = h, MOD = MOD)$root[1]
  msy <- yield_f(FF = Fmsy, h = h, MOD = MOD)$y
  Fcrash <- bisect()
  brps[MOD,] <- c(MOD_NAMES[MOD], round(Fmsy, 4), round(msy, 4), round(Fcrash, 4))
}


for(MOD in 1:3){
  
  out <- yield_f(FF, h, MOD = MOD)
  refs <- brps[MOD,]
  
  par(mar=c(4,4,3,3), mfrow = c(2,2))
  
  plot(out$s_ratio, out$r, type = "l", xlab = "S(F)", ylab = "R(F)",
       xlim = c(0, max(out$s_ratio)), ylim = c(0, max(out$r)*1.1))
  abline(0,1, col = "grey", lty = 2)
  
  plot(out$s_ratio, out$y, type = "l", xlab = "S(F)", ylab = "Y(F)",
       xlim = c(0, max(out$s_ratio)), ylim = c(0, max(out$y)*1.1))
  
  text(line2user(line=mean(par('mar')[c(2, 4)]), side=2),
       line2user(line=2, side=3), MOD_NAMES[MOD], xpd=NA, cex=2, font=2)
  
  plot(FF, out$y, type = "l", xlab = "F", ylab = "Y(F)",
       xlim = c(0, as.numeric(refs$Fcrash)*1.2), ylim = c(0, max(out$y)*1.1))
  abline(h = refs$MSY, col = "grey", lty = 2)
  abline(v = refs$Fmsy, col = "grey", lty = 2)
  abline(v = refs$Fcrash, col = "grey", lty = 2)
  
  plot(FF, out$s_ratio, type = "l", xlab = "F", ylab = "S(F)",
       xlim = c(0, as.numeric(refs$Fcrash)*1.2), ylim = c(0, max(out$s_ratio)*1.1))
  abline(v = refs$Fcrash, col = "grey", lty = 2)
  
}

# Question 3b
kable(brps) %>% kable_styling("striped", full_width = F)

# # Compare the ratio S(Fmsy)/S(0) across different values of h for each SR model
h_vals <- seq(0.25, 0.95, 0.05)
FF <- seq(0, 1, 0.01)

out <- list()

for(j in 1:3){
  
  res <- list()
  
  for(i in 1:length(h_vals)) {
    
    h <- h_vals[i]   
    
    tmp <- yield_f(FF, h = h, MOD = j)
    
    Fmsy <- uniroot(f = dfx.dx, interval = c(0.01,1), tol = 0.000001, 
                    yield_f = yield_f, h = h, MOD = j, delta = 0.001)$root[1]
    
    s_Fmsy <- yield_f(FF = Fmsy, h = h, MOD = j)$s
    
    s0 <- yield_f(FF = Fmsy, h = h, MOD = j)$s0
    
    ratio <- rep(s_Fmsy/s0, length(FF))
    
    res[[i]] <- data.frame(Model = rep(MOD_NAMES[j], length(FF)),
                           h = rep(h, length(FF)),
                           ratio = rep(ratio, length(FF)),
                           s = tmp$s_ratio,
                           r = tmp$r)
  }
  
  out[[j]] <- do.call(rbind, res)
}

out <- do.call("rbind", out)

ggplot(out, aes(s, r, size = h, group = h)) +
  geom_line(colour = "grey50") +
  scale_size_continuous(range = c(.1,1), breaks = seq(0.25,0.95, 0.1),
                        labels = paste0(seq(0.25, 0.95, 0.1))) +
  scale_x_continuous(limits = c(0, max(out$s))) +
  scale_y_continuous(limits = c(0, max(out$r))) +
  facet_wrap(~Model) +
  labs(x = "S(F)", y = "R(F)")

out %>%
  distinct(Model, h, ratio) %>%
  mutate(ratio = round(ratio, 4)) %>%
  select(Model, h, `S(Fmsy)/S(0)` = ratio) %>%
  dcast(h ~ Model) %>%
  kable() %>%
  kable_styling("striped", full_width = F) %>%
  add_header_above(c(" " = 1, "S(Fmsy)/S(0)" = 3))

# # The target depletion for many west coast groundfish, including widow rockfish,
# # is 0.4*S(0). Comment on the implications of these results for management.
# # every time the ratio S(Fmsy)/S(0) exceeds 0.4, the target depletion exceeded
# # Fmsy.
# 
# # Percentage of times each model was on target or resulted in overfishing:
out %>%
  group_by(Model) %>%
  summarise(`On or below target` = paste0(round(length(ratio[ratio <= 0.4])/length(ratio) * 100, 1), "%"),
            `Overfishing` = paste0(round(length(ratio[ratio > 0.4])/length(ratio) * 100, 1), "%")) %>%
  kable() %>%
  kable_styling("striped", full_width = F)
