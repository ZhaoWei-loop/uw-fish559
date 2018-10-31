# Assignment: write function that take information of natural mortality,
# fecundity-at-age, weight-at-age, selectivity, R0, steepness, the form of the
# stock-recruitment relationships, and gamma
library(knitr)
library(ggplot2)
library(dplyr)
library(kableExtra)
library(data.table)
# library(cowplot)

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
      N[f,1] <- 0.5 #/2 
      
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

dat <- read.table("Home3/HOME3.txt", header = TRUE)
colnames(dat) <- c("age", "fec", "wt_f", "wt_m", "sel_f", "sel_m")

F_vals <- seq(0, 1, 0.001) # vector of fishing mortality
h <- 0.5 # steepness parameter
r0 <- 1 # unfished recruitment
gamma <- 1 # third parameter for pella-tomlinson model

# Get matrices of sex-specific total instantaneous mortalities 
Zf <- get_Z(age = dat$age, sel = dat$sel_f, F_vals = F_vals)
Zm <- get_Z(age = dat$age, sel = dat$sel_m, F_vals = F_vals)

# Get sex-specific numbers-at-age matrices
Nf <- get_N(Z = Zf, age = dat$age)
Nm <- get_N(Z = Zm, age = dat$age)

# Get yield vectors different values of F
ypr_f <- get_ypr(age = dat$age, wt = dat$wt_f, Z = Zf, N = Nf)   
ypr_m <- get_ypr(age = dat$age, wt = dat$wt_m, Z = Zm, N = Nm)   
ypr <- ypr_f + ypr_m

# Get spawner biomass-per-recruit across different values of F
spr <- get_spr(age = dat$age, fec = dat$fec, N = Nf)
spr0 <- spr[1] # F=0
s0 <-  r0*spr0 # unfished ssb 

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

# Get plots function
get_results <- function() {
  
  # Get stock-recruitment parameters
  ab <- get_sr_params(rec_mod = MOD, spr0 = spr0, h = h, r0 = r0, gamma = gamma)
  
  # Get recruitment values
  r <- get_recruit(rec_mod = MOD, spr = spr, s0 = s0, alpha = ab[[1]], beta = ab[[2]], gamma = gamma)
  
  # Yield
  y <- ypr*r 
  
  par(mar=c(4,4,3,3), mfrow = c(2,2))
  
  plot(spr*r, r, type = "l", xlab = "S(F)", ylab = "R(F)",
       xlim = c(0, max(spr*r)), ylim = c(0, max(r)*1.1))
  abline(0,1, col = "grey", lty = 2)
  
  plot(spr*r, y, type = "l", xlab = "S(F)", ylab = "Y(F)",
       xlim = c(0, max(spr*r)), ylim = c(0, max(y)*1.1))
  
  text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
       line2user(line=2, side=3), MOD_NAME, xpd=NA, cex=2, font=2)
  
  # MSY, Fmsy and Fcrash
  imax <- match(max(y), y)
  icrash <- match(max(y[y < 0]), y)
  out <- data.frame(Model = MOD_NAME,
                    MSY = y[imax],
                    Fmsy = F_vals[imax],
                    Fcrash = F_vals[icrash])
  
  plot(F_vals, y, type = "l", xlab = "F", ylab = "Y(F)",
       xlim = c(0, F_vals[icrash]*1.2), ylim = c(0, max(y)*1.1))
  abline(h = out$MSY, col = "grey", lty = 2)
  abline(v = out$Fmsy, col = "grey", lty = 2)
  abline(v = out$Fcrash, col = "grey", lty = 2)
  
  plot(F_vals, spr*r,type = "l", xlab = "F", ylab = "S(F)",
       xlim = c(0, F_vals[icrash]*1.2), ylim = c(0, max(spr*r)*1.1))
  abline(v = out$Fcrash, col = "grey", lty = 2)
  
  return(out)
}

# Question 3a
MOD <- 1; MOD_NAME <- "Beverton-Holt stock-recruitment"
out <- get_results()

MOD <- 2; MOD_NAME <- "Ricker stock-recruitment"
out2 <- get_results()

MOD <- 3; MOD_NAME <- "Pella-Tomlinson stock-recruitment"
out3 <- get_results()

# Question 3b
out <- rbind(out,out2,out3)
out$MSY <- round(out$MSY, 2)
kable(out)  %>%
  kable_styling("striped", full_width = F) 

# Compare the ratio S(Fmsy)/S(0) across different values of h for each SR model
h_vals <- seq(0.25, 0.95, 0.05)

get_ratio <- function() {
  out <- list() 
  for(i in 1:length(h_vals)) {
    # Get stock-recruitment parameters
    ab <- get_sr_params(rec_mod = MOD, spr0 = spr0, h = h_vals[i], r0 = r0, gamma = gamma)
    # Get recruitment values
    r <- get_recruit(rec_mod = MOD, spr = spr, s0 = s0, alpha = ab[[1]], beta = ab[[2]], gamma = gamma)
    # Yield
    y <- ypr*r 
    # Spawning stock biomass
    s <- spr*r
    
    target <- round(0.4*s0,2) # target depletion according to assignment
    itarget <- match(target, round(s,2)) # associated S value
    Ftarget <- F_vals[itarget] # associated F value
    
    # Fmsy and S_Fmsy and 
    imax <- match(max(y), y)
    F_msy <- F_vals[imax]
    s_Fmsy <- s[imax]

    out[[i]] <- data.frame(mod = rep(MOD_NAME, length(s)),
                           h = rep(h_vals[i], length(s)),
                           r = r,
                           s = s,
                           s_target = rep(target, length(s)),
                           Ftarget = rep(Ftarget, length(s)),
                           Fmsy = rep(F_msy, length(s)),
                           s_Fmsy = rep(s_Fmsy, length(s)),
                           ratio = rep(s_Fmsy/s0, length(s)))
  }
  return(out)
}

par(mar=c(4,4,3,3), 
    mfrow = c(1,3))

MOD <- 1; MOD_NAME <- "Beverton-Holt stock-recruitment"
out <- get_ratio() 
out <- do.call("rbind", out)

MOD <- 2; MOD_NAME <- "Ricker stock-recruitment"
out2 <- get_ratio()
out2 <- do.call("rbind", out2)

MOD <- 3; MOD_NAME <- "Pella-Tomlinson stock-recruitment"
out3 <- get_ratio()
out3 <- do.call("rbind", out3)

out <- rbind(out,out2,out3)

ggplot(out, aes(s, r, size = h, group = h)) +
  geom_line(colour = "grey50") +
  scale_size_continuous(range = c(.1,1), breaks = seq(0.25,0.95, 0.1),
                        labels = paste0(seq(0.25, 0.95, 0.1))) +
  scale_x_continuous(limits = c(0, max(out$s))) +
  scale_y_continuous(limits = c(0, max(out$r))) +
  facet_wrap(~mod) +
  labs(x = "S(F)", y = "R(F)")

out %>% 
  distinct(mod, h, ratio) %>% 
  mutate(ratio = round(ratio, 2)) %>% 
  select(mod, h, `S(Fmsy)/S(0)` = ratio) %>% 
  dcast(h ~ mod) %>% 
  kable() %>% 
  kable_styling("striped", full_width = F) %>% 
  add_header_above(c(" " = 1, "S(Fmsy)/S(0)" = 3))

# The target depletion for many west coast groundfish, including widow rockfish,
# is 0.4*S(0). Comment on the implications of these results for management.
# every time the ratio S(Fmsy)/S(0) exceeds 0.4, the target depletion exceeded
# Fmsy.

# Percentage of times each model was on target or resulted in overfishing:
out %>% 
  group_by(Model = mod) %>% 
  summarise(`On or below target` = paste0(round(length(ratio[ratio <= 0.4])/length(ratio) * 100, 1), "%"),
            `Overfishing` = paste0(round(length(ratio[ratio > 0.4])/length(ratio) * 100, 1), "%")) %>% 
  kable() %>% 
  kable_styling("striped", full_width = F)
  