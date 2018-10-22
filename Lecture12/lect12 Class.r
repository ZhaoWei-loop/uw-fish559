Lect12<-function()
{
 set.seed(128)
 ParVals <- GenData(Nsim=500)
 sumpars(ParVals)

 set.seed(128)
 ParVals <- GenData(Nsim=50)
 sumpars(ParVals)
 
 ProjRec1 <- DoProject(ParVals$r,ParVals$K,ParVals$z,ParVals$N1,ParVals$Like,100,0,CVabund=0.2,DoPlot=T,RecY0=NULL)
 TestS <- rep(c(1,2),50)
 ProjRec2 <- DoProject(ParVals$r,ParVals$K,ParVals$z,ParVals$N1,ParVals$Like,100,1,TheCat=TestS,CVabund=0.2,DoPlot=T,RecY0=ProjRec1$recYear)
 TestS <- c(rep(1,10),rep(2,90))
 ProjRec3 <- DoProject(ParVals$r,ParVals$K,ParVals$z,ParVals$N1,ParVals$Like,100,1,TheCat=TestS,CVabund=0.2,DoPlot=T,RecY0=ProjRec1$recYear)
 ProjRec4 <- DoProject(ParVals$r,ParVals$K,ParVals$z,ParVals$N1,ParVals$Like,100,2,TheCat=TestS,CVabund=0.2,DoPlot=T,RecY0=ProjRec1$recYear)
 NULL
    
}

NmsyKz <- function(z,NmsyK) { 1-(z+1)*NmsyK^z }

# =========================================================================================================================================================================

GenData<-function(Nsim=100)
{
 print(Nsim)    
 
# Create the vectors of parameters
 NmsyK <- runif(Nsim,0.5,0.7)
 K <- runif(Nsim,1000,2000)
 N1 <- rnorm(Nsim,500,50)
 MSYNmsy <- runif(Nsim,0.01,0.05)
 z <- rep(0,Nsim)

 # Solve for z given NMsy/K
 for (I in 1:Nsim)
  {
    xx <- uniroot(NmsyKz,NmsyK=NmsyK[I],lower=1,upper=100)
    z[I] <- xx$root
    print(c(z[I],NmsyK[I],NmsyKz(z[I],NmsyK[I])))
  }

 # Compute r
 r <- MSYNmsy/(1-NmsyK^z)

 Ncurr <- N1
 for (Iyear in 1:14)
   Ncurr <- Ncurr + r*Ncurr*(1-(Ncurr/K)^z) - 1
 Like <- exp(-1*(Ncurr-625)^2/(2*75*75))
 TotLike <- sum(Like)
 Like <- Like / sum(Like)*100

 # Output
 Params <- NULL
 Params$r <- r
 Params$K <- K
 Params$z <- z
 Params$N1 <- N1
 Params$Like <- Like
 Params$NmsyK <- NmsyK
 Params$MSYNmsy <- MSYNmsy
 return(Params)
}

# =========================================================================================================================================================================

sumpars<-function(ParVals)
 {


 Npost = 1000
 Nsim <- length(ParVals$N1)
 Post <- matrix(0,ncol=4,nrow=Npost)

 }
# =========================================================================================================================================================================
DoProject <- function(r,K,z,N1,Like,NprojYear=100,Ctype=0,RecDef=0.8,TheCat=NULL,Rmax=0.02,Fr=0.5,CVabund=0.2,DoPlot=F,RecY0)
{
 # Extract Nsim
 Nsim <- length(r)

 Simout <- matrix(0,ncol=Nsim,nrow=(16+NprojYear))
 Catches <- rep(1,15+NprojYear)
 recYear <- rep(0,Nsim)

 # Set up the numbers and project to the start of the first "real" year
 Simout[1,] <- N1
 for (Iyear in 1:15)
   Simout[Iyear+1,] <- Simout[Iyear,] + r*Simout[Iyear,]*(1-(Simout[Iyear,]/K)^z) - 1
   
 # Project forward (simulation by simulation)
 for (Isim in 1:Nsim)
  {
   found <- F
   for (Iyear in 1:NprojYear)
     {
      IyrAct <- Iyear+15    
    
      # Check for recovery
      if (Simout[IyrAct,Isim] > RecDef*K[Isim] && found == F) 
          { found <- T; recYear[Isim] <- IyrAct }
    
      # Catch strategy
      if (Ctype==0) Catch <- 0  
      if (Ctype==1) Catch <- TheCat[Iyear]
      if (Ctype==2)
        {
         # PBR happens every 5th year
         if (Iyear %% 5 == 1)
           {
            Nval <- Simout[IyrAct,Isim]*exp(rnorm(1,0,CVabund))
            Nminq <- qnorm(0.2,0,1)
            Nmin <- Nval*exp(Nminq*CVabund)
            Catch <- 0.5*Rmax*Nmin*Fr
            if (Catch < 0) Catch = 0
           }
         }
    
      # Store the catch
      Catches[IyrAct] <- Catch 
    
      # Project forward
      Simout[IyrAct+1,Isim] <- Simout[IyrAct,Isim] + r[Isim]*Simout[IyrAct,Isim]*(1-(Simout[IyrAct,Isim]/K[Isim])^z[Isim]) - Catch
      if (Simout[IyrAct+1,Isim] < 1) Simout[IyrAct+1,Isim] <- 1
      }  
      Simout[,Isim] <- Simout[,Isim]/K[Isim] 
     if (found == F)
      recYear[Isim] <- 15+NprojYear
    }

 # Results
 # print(recYear)
 if (Ctype!=0)
  {
   DiffProb <- 0
   DiffProb <- sum(Like * (recYear-RecY0)/(RecY0-15))
   print(c(paste("Percentage difference in recovery probability ",format(DiffProb,digits=5) )))
  }

  
 if (DoPlot)
  {
   par(mfrow=c(2,2))
   
   # Histogram opf recovery
   hist(recYear,xlab="Recovery Year")

   # catches (last simulation)
   Years <-c(1:(15+NprojYear))
   plot(Years,Catches,xlab="Year",ylab="Catch",type='l',lty=1,lwd=2) 

   # Time-trajectories of depletion
   Y <- rep(1.0,15+NprojYear)
   maxy <- 1.05
   for (Iyear in 1:(15+NprojYear))
    Y[Iyear] <- quantile(Simout[Iyear,],prob=0.95)
   plot(Years,Y,xlab="Year",ylab="Propulation Depletion",type='n',ylim=c(0,maxy))
   lines(Years,Y,lty=2,lwd=2)
   for (Iyear in 1:(15+NprojYear))
    Y[Iyear] <- quantile(Simout[Iyear,],prob=0.5)
   lines(Years,Y,lty=1,lwd=2)
   for (Iyear in 1:(15+NprojYear))
    Y[Iyear] <- quantile(Simout[Iyear,],prob=0.05)
   lines(Years,Y,lty=2,lwd=2)
   abline(RecDef,0,lty=3,lwd=3)
  }
 NULL

 results <- NULL
 results$recYear <- recYear
 return(results)
 NULL
}
# =========================================================================================================================================================================

Lect12()
