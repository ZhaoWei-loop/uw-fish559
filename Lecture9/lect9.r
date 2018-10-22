setwd("Lecture9")

library(stats4)
library(MASS)
data(faithful)
waiting <- faithful[,2]



Lect8b<-function()
{
 Case1() # Optim   
# Case2() # Simplex
# Case3() # Powell
#  Case5() # mle

}

# =====================================================================================================================================================================

Lect8a <-function()
{

 # Specify the function to minimize
 # ================================
 mix.obj <- function(p,x) {
    e <- p[1]*dnorm((x-p[2])/p[3])/p[3]+(1-p[1])*dnorm((x-p[4])/p[5])/p[5]
    return(- sum(log(e)))
 }

 # Vector function 
 tstfn <- function(x,p,u1,s1,u2,s2) p*dnorm(x,u1,s1)+(1-p)*dnorm(x,u2,s2)

 # Plot the histogram and add a line
 par(mfrow=c(2,2))  
 truehist(waiting,xlim=c(35,110),h=5,ymax=0.07) 
 wait.dns <- density(waiting,200,width=10.24)
 lines(wait.dns$x,tstfn(wait.dns$x,0.3,50,5,80,5),lty=1,lwd=3,col=2)
 lines(wait.dns$x,tstfn(wait.dns$x,0.3,50,5,80,5),lty=1,lwd=3,col=2)
 lines(wait.dns$x,tstfn(wait.dns$x,0.3,50,5,80,5),lty=1,lwd=3,col=2)

 # Fit the modle using optim
 wait.init <- c(0.3,50,5,80,5)
 mix.n1 <- optim(wait.init,mix.obj,method="L-BFGS-B",lower=c(0,-Inf,0,-Inf,0),upper=c(1,rep(Inf,4)),x=waiting)
 print(mix.n1)
 truehist(waiting,xlim=c(35,110),h=5,ymax=0.07) 
 wait.dns <- density(waiting,200,width=10.24)
 lines(wait.dns$x,tstfn(wait.dns$x,mix.n1$par[1],mix.n1$par[2],mix.n1$par[3],mix.n1$par[4],mix.n1$par[5]),lty=1,lwd=3,col=2)


 # Add in the derivartives
 # =======================
 lmix2a <- deriv( ~-log(p*dnorm((x-u1)/s1)/s1 + (1-p)*dnorm((x-u2)/s2)/s2), c("p","u1","s1","u2","s2"), function(x,p,u1,s1,u2,s2) NULL)
  
 # gradient function
 mix.gr <- function(p,x) {
  u1 <- p[2]; s1 <- p[3]; u2 <- p[4]; s2 <- p[5]; p <- p[1]
  e <- lmix2a(x,p,u1,s1,u2,s2)
  print(e)
  print(rep(1,length(x)) %*% attr(e,"gradient"))
  rep(1,length(x)) %*% attr(e,"gradient")
  
   }
 
 #supplying the gradient function - computes a gradient for each function call
 mix.n11 <-optim(wait.init,mix.obj,gr=mix.gr,method="L-BFGS-B",lower=c(0,-Inf,0,-Inf,0),
          upper=c(1,Inf,Inf,Inf,Inf),x=waiting)
 print(mix.n11)

}

# =====================================================================================================================================================================

Sign<-function(a,b)
{
 if (b > 0)
  return(abs(a))
 else
  return(-1*abs(a))
}

# =====================================================================================================================================================================

Fit<-function(Xinit,Ndim,FUNK,doGraph=F)
{
 # Set up some basics
 Tol <- 0.00000001
 GRD <- 1.2

 # Set up storage
 P <- matrix(0,Ndim+1,Ndim)
 XO <- rep(0,Ndim)
 Y <- rep(0,Ndim)

 # Set up the initial simplex
 PVecs <- NULL
 PVals <- NULL
 for (I in 1:(Ndim+1))
  {
    P[I,] <- Xinit
    if (I > 1) P[I,I-1] <- GRD*P[I,I-1]
    XO <- P[I,] 
    Y[I] <- FUNK(XO)
    PVals <- c(PVals,Y[I])
    PVecs <- c(PVecs,XO)
  }
 print(PVecs)
 
 ZZ <- Amoeba(P,Y,Ndim,Tol,FUNK,PVals,PVecs,T)
 print(ZZ)
 print(ZZ$Zfin)
 print(ZZ$Sol)

 if (doGraph)
  {
   par(mfrow=c(2,2))
   X <- seq(1,ZZ$Nfunc,1)
   Y <- ZZ$TraceY
   plot(X,log(Y),xlab="Function call",ylab="log(Function value)",type='b',pch=16,lwd=2)
   X <- NULL; Y <-NULL
   for (I in 1:ZZ$Nfunc)
    {
     X <- c(X,ZZ$TraceX[(I-1)*2+1])
     Y <- c(Y,ZZ$TraceX[(I-1)*2+2])
    }
    print(X)
    print(Y)
    plot(exp(X),exp(Y),xlim=c(500,5000),ylim=c(0.2,0.4),type='n',lwd=2,xlab="K",ylab="r")
   lines(exp(X),exp(Y),col=5,lty=1,lwd=2,type='b',pch=16)
 }
 return(ZZ)
}

# ==============================================================================================================================================================================

Amoeba<-function(P,Y,Ndim,Ftol,FUNK,PVals,PVecs,doPlot)
{
 # Fixed parameters (expansion, contraction, extra expansion)
 ALPHA <- 1; BETA <- 0.5; GAMMA <- 2; Mpnts <- Ndim+1   
    
 # Temporary storage    
 PR <- rep(0,Ndim)
 PRR <- rep(0,Ndim)
 Pbar <- rep(0,Ndim)
    
 print("In")    
 print(Ndim)
 print(Ftol)
 print(P)
 print(Y)
 print("here")
 
 if (doPlot) par(mfrow=c(3,4))
 for (Iter in 1:500)
  {

   # Find the highest (worst), next-highest, and lowest (best)
   ILO <- 1;
   if (Y[1] > Y[2]) 
    { IHI <- 1; INHI <- 2 }
   else
    { IHI <- 2; INHI <- 1 }     
   for (I in 1:Mpnts)
   {
     if (Y[I] < Y[ILO]) ILO <- I
     if (Y[I] > Y[IHI])
      {
        INHI <- IHI; IHI <- I 
      }
     else   
      if ((Y[I] > Y[INHI]) & (I != IHI)) INHI <- I
   }  

    # Check for convergence
    Rtol <- 2.0*abs(Y[IHI]-Y[ILO])/(abs(Y[IHI]) + abs(Y[ILO]))
    if (Rtol < Ftol) break;
    
    # Begin a new iteration (set up the average of the NON-highest points
        Pbar <- rep(0,Ndim)
    for (I in 1:Mpnts)
     if (I != IHI)
      {
       for (J in 1:Ndim)    
         Pbar[J] <- Pbar[J] + P[I,J]
      }     
    # Reflect the simplest through the highest point and compute the function
    for (J in 1:Ndim)
     {
      Pbar[J] <- Pbar[J] / Ndim
      PR[J] <-  (1.0 + ALPHA)*Pbar[J] - ALPHA*P[IHI,J]
     }  
    YPR <- FUNK(PR) 
        PVecs <- c(PVecs,PR)
        PVals <- c(PVals,YPR)
    
    # Is the new point better than the last
    if (YPR <= Y[ILO])
     {
    # Improved so keep trying in this direction
      for (J in 1:Ndim)
       PRR[J] <- GAMMA*PR[J] + (1-GAMMA)*Pbar[J]
      YPRR <- FUNK(PRR)
          PVecs <- c(PVecs,PRR)
          PVals <- c(PVals,YPRR)

      # Keep PRR or PR depending on which leads to a lower function value   
      if (YPRR < Y[ILO])
       { P[IHI,] <- PRR; Y[IHI] <- YPRR }
      else
       { P[IHI,] <- PR; Y[IHI] <- YPR }
     }   
    else
      # Whoops, at least as bad as the second worst point
     if (YPR >= Y[INHI]) 
      {

       # Replace the worst point if this point is worse than the worst
       if (YPR < Y[IHI])
        {P[IHI,] <- PR; Y[IHI] <- YPR }

       # contract - this wasn't the way forward
        for (J in 1:Ndim)
         PRR[J] <- BETA*P[IHI,J] + (1 - BETA)*Pbar[J]
       YPRR <- FUNK(PRR) 
       PVecs <- c(PVecs,PRR)
       PVals <- c(PVals,YPRR)

       # contraction gave an improvement
       if (YPRR < Y[IHI])
        { P[IHI,] <- PRR; Y[IHI] <- YPRR }
       else
        {
         # Nothing seemed to help so contract towards the lowest point
         for (I in 1:Mpnts)   
           if (I != ILO)
           {
             for (J in 1:Ndim)
              { PR[J] <- 0.5*(P[I,J]+P[ILO,J]); 
                P[I,J] <- PR[J] }
             Y[I] <- FUNK(PR)
            PVecs <- c(PVecs,PR)
            PVals <- c(PVals,Y[I])
           }
        }
      }
   else
     # Not best but better than worst so replce
    {P[IHI,] <- PR; Y[IHI] <- YPR }
   ItOut <- Iter
  }
    
 Xout <- P[ILO,]
 Sout <- FUNK(Xout)

 # Return
 ZZ <- NULL
 ZZ$Zfin <- Xout
 ZZ$Sol <- Sout
 ZZ$TraceY <- PVals
 ZZ$TraceX <- PVecs
 ZZ$Dim <- Ndim
 ZZ$NIter <- ItOut
 ZZ$Nfunc <- length(PVals)
 return(ZZ)
}

# =====================================================================================================================================================================

Fit2<-function(Xinit,Ndim,FUNK,doGraph)
{
 # Tolerances, etc,
 Ftol <- 0.00000001
    
 # Set up a matrix of initial directions    
 XI <- matrix(0,Ndim,Ndim)  
 for (I in 1:Ndim) XI[I,I] <- 1
  
 ZZ <- Powell(Xinit,XI,Ndim,FUNK,Ftol)
 print(c("final",ZZ))

 if (doGraph)
  {
   par(mfrow=c(2,2))
   Nfunc <- length(ZZ$TraceY)
   X <- seq(1,Nfunc,1)
   Y <- ZZ$TraceY
   plot(X,log(Y),xlab="Function call",ylab="log(Function value)",type='b',pch=16,lwd=2)
   X <- NULL; Y <-NULL
   Npnt <- length(ZZ$TraceX)/2
   for (I in 1:Npnt)
    {
     X <- c(X,ZZ$TraceX[(I-1)*2+1])
     Y <- c(Y,ZZ$TraceX[(I-1)*2+2])
    }
    print(X)
    print(Y)
    plot(exp(X),exp(Y),xlim=c(500,5000),ylim=c(0.2,0.4),type='n',lwd=2,xlab="K",ylab="r")
   lines(exp(X),exp(Y),col=5,lty=1,lwd=2,type='b',pch=16)
 }
 return(ZZ)

}

# =====================================================================================================================================================================

Powell<-function(X,XI,N,FUNK,Ftol)
{
 # Declarations
 XIT <- rep(0,N)    
    
 # Prleiminaries 
 PT <- rep(0,N)
 Fret <- FUNK(X)
 PT <- X
 FVecs <- X
 FVals <- Fret

 for (Iter in 1:10)
  {
    FP <- Fret
    IBig <- 0
    Del  <- 0
    
    # Do a line minimization in each direction
    for (I in 1:N)
     {   
      XIT <- XI[,I]
      FPTT <- Fret
      ZZ<- LinMin(X,XIT,N,FUNK)
     FVals <- c(FVals,ZZ$FVals)
     FVecs <- c(FVecs,ZZ$X)
     Fret <- ZZ$Fret; X <- ZZ$X; XIT <- ZZ$XIT

      # Record the difference and store the largest different
      if (abs(FPTT-Fret) > Del)
       { Del <- abs(FPTT-Fret); IBIG <- I }
    
     }
    
    # Any progres, else quit
    if (2.0*abs(FP-Fret) <= Ftol*(abs(FP)+abs(Fret))) break
    
    # construnct the extrapolated point and the average direction moved
    PTT <- 2*X - PT
    XIT <- X - PT
    PT <- X
    
    FPTT <- FUNK(PTT)
    FVals <- c(FVals,FPTT)
    FVecs <- c(FVecs,PTT)
    if (FPTT < FP)
     {
      TT <- 2*(FP-2*Fret+FPTT)*(FP-Fret-Del)^2-Del*(FP-FPTT)^2
      print(TT)
      if (TT < 0)
       {
         ZZ <- LinMin(X,XIT,N,FUNK)
        FVals <- c(FVals,ZZ$FVals)
        FVecs <- c(FVecs,ZZ$X)
         Fret <- ZZ$Fret; X <- ZZ$X; XIT <- ZZ$XIT
         XI[,IBIG] <- XIT
       }    
     }  
   }
    
 ZZ <- NULL
 ZZ$Zfin <- X
 ZZ$Sol <- FUNK(X)
 ZZ$TraceY <- FVals
 ZZ$TraceX <- FVecs
 return(ZZ)     
}

# =====================================================================================================================================================================

LinMin <- function(X,XI,N,FUNK)
{
 TOL <- 1.0E-4  

 Mout <- MnBrak(0,1,F1DIM,FUNK,X,XI)
 Bout <- Brent(Mout$AX,Mout$BX,Mout$CX,F1DIM,FUNK,TOL,X,XI)
 Xmin <- Bout$X
 Fret <- Bout$FX
 XI <- Xmin * XI
 X <- X + XI
 
 ZZ <- NULL
 ZZ$Fret <- FUNK(X)
 ZZ$XIT <- XI
 ZZ$X <- X
 ZZ$FVals <- c(Mout$FVals,Bout$FVals)
 ZZ$FVecs <- Bout$FVecs
 return(ZZ)
}

# =====================================================================================================================================================================

F1DIM <- function(X,Pcom,XIcom,FUNK)
{
    
 # Form the vector  
 XT <- Pcom + X*XIcom   
 ZZ <- FUNK(XT)
 return(ZZ)

}

# =====================================================================================================================================================================

Brent <- function(AX,BX,CX,F1DIM,FUNK,TOL,Pcom,XIcom)
{
 # Defaults
 CGold <- 0.3819660; Zeps <- 1.0E-10    

 # Initial function value counter
 FVals <- NULL
 FVecs <- NULL
    
 # Initializations
 A <- min(AX,CX)
 B <- max(AX,CX)
 V <- BX; W <- V; X<- V; E <- 0
 FX <- F1DIM(X,Pcom,XIcom,FUNK)
 FV <- FX; FW <- FX;
 FVecs <- c(FVecs,X)
 FVals <- c(FVals,FX)

 for (Iter in 1:100)
  {
    XM <- 0.5*(A+B)
    TOL1 <- TOL*abs(X)+Zeps
    TOL2 <- 2.0*TOL1
    
    # Test for completeness
    if (abs(X-XM) <= (TOL2-0.5*(B-A))) break
        
    # Construct a trial parabolic fit   
    Ok <- T
    if (abs(E) > TOL1)  
     {
      R <- (X-W)*(FX-FW); Q <- (X-V)*(FX-FW); P <- (X-V)*Q-(X-W)*R
      Q <- 2.0*(Q-R)
      if (Q > 0) P <- -1*P
      Q <- abs(Q)
      Etemp <- E; E <- D
      Ok <- T
      if (abs(P) < abs(0.5*Q*Etemp) && (P > Q*(A-X)) && P <= Q*(B-X))
       {
         D <- P/Q; U <- X +D
         if (U-A < TOL2 || B-U < TOL2) D <- Sign(TOL1,XM-X)
         Ok <- F
       }    
     }  

   # Gold section step
   if (Ok)
    {
      if (X >= XM) 
       E <- A - X 
      else 
       E <- B - X
      D <- CGold*E
    }

   # Proceed either from gold step or parabolic fit
   if (abs(D) >= TOL1) 
    U <- X+D
   else
    U <- X + Sign(TOL1,D)
   FU <- F1DIM(U,Pcom,XIcom,FUNK)
   FVecs <- c(FVecs,U)
   FVals <- c(FVals,FU)
   if (FU <= FX)
    {
      if (U >= X)
       A <- X
      else
       B <- X
      V <- W; FV <- FW; W <- X; FW <- FX; X <- U; FX <- FU
    
    }
   else
    {
      if (U < X)
       A <- U
      else
       B <- U
      if (FU <= FW || W == X)
       { V <- W; FV <- FW; W <- U; FW <- FU }
      else
       {
         if (FU <= FV || V == X || V == W)
          { V <- U; FV <- FU }
       }
    }
    
  } # for (Iter) 
    
 ZZ <- NULL
 ZZ$X <- X
 ZZ$FX <- FX
 ZZ$FVals <- FVals
 ZZ$FVecs <- FVecs
 return(ZZ)
}

# =====================================================================================================================================================================

MnBrak<-function(AX,BX,F1DIM,FUNK,Pcom,XIcom)
{
 # various fixed parameters
 Gold <- 1.618034; Glimit <- 100; Tiny <- 1.0E-20   

 # switch A and B so that we can go DONWHILL in the direction from A to B; 
 # compute a first guess for C
 FA <- F1DIM(AX,Pcom,XIcom,FUNK)    
 FB <- F1DIM(BX,Pcom,XIcom,FUNK)    
 FVals <- c(FA,FB)

 if (FB > FA)
  { Dum <- AX; AX <- BX; BX <- Dum; Dum <- FB; FB <- FA; FA <- Dum }
 CX <- BX + Gold*(BX-AX)
 FC <- F1DIM(CX,Pcom,XIcom,FUNK)
 FVals <- c(FVals,FC)

 while(FB >= FC)
  {
   R <- (BX-AX)*(FB-FC)
   Q <- (BX-CX)*(FB-FA)
    U <- BX - ((BX-CX)*Q-(BX-AX)*R)/(2.0*Sign(max(abs(Q-R),Tiny),Q-R))
    Ulim <- BX + Glimit*(CX-BX)
      
   # try various possibilities
    if ((BX-U)*(U-CX) > 0) 
     {
      FU <- F1DIM(U,Pcom,XIcom,FUNK)
     FVals <- c(FVals,FU)
      if (FU < FC)
       { AX <- BX; FA <- FB; BX <- U; FB <- FU; break }
      else
        if (FU > FB)
         { CX <- U; FC <- FU; break}
      U <- CX + Gold*(CX-BX)
      FU <- F1DIM(U,Pcom,XIcom,FUNK)
      FVals <- c(FVals,FU)
     }
    else
      {
       # parabolic fit between C and its allowed limit
       if ((CX-U)*(U-Ulim) > 0)
         {
          FU <- F1DIM(U,Pcom,XIcom,FUNK)
         FVals <- c(FVals,FU)
          if (FU < FC)
            { BX <- CX; CX <- U; U <- CX + Gold*(CX-BX); FB <- FC; FC <- FU; FU <- F1DIM(U,Pcom,XIcom,FUNK); FVals <- c(FVals,FU) }
        }
        else
         {
          if ((U-Ulim)*(Ulim-CX) >=0)
            { U <- Ulim; FU <- F1DIM(U,Pcom,XIcom,FUNK);     FVals <- c(FVals,FU) }
          else
            { U <- CX + Gold*(CX-BX);   FU <- F1DIM(U,Pcom,XIcom,FUNK); FVals <- c(FVals,FU) }  
         }
     } # if ((BX-U)*(U-CX) > 0) 

     # Elimiate oldest point and continue
     AX <- BX; BX <- CX; CX <- U; FA <- FB; FB <- FC; FC <- FU
    
  } # while

 ZZ<-NULL
 ZZ$AX <- AX
 ZZ$BX <- BX
 ZZ$CX <- CX
 ZZ$FA <- FA
 ZZ$FB <- FB
 ZZ$FC <- FC
 ZZ$FVals <- FVals
 return(ZZ)

}
# =====================================================================================================================================================================
# =====================================================================================================================================================================

f1 <- function(x)
{
 K <- exp(x[1])
 r <- exp(x[2])

 # Temp storage
 CpuePred <- rep(0,co$Nyear)
 Biomass <- rep(0,co$Nyear+1)

 # Do projections
 Biomass[1] <- K
 for (Year in 1:co$Nyear)
  {
    Biomass[Year+1] <- Biomass[Year] + r*Biomass[Year]*(1.0-Biomass[Year]/K) - co$Catch[Year]
    if (Biomass[Year+1] < 0.01) Biomass[Year+1] <- 0.01
  }

 # Calculate the ML estimate of q
 qbar <- 0; npar <- 0
 for (Year in 1:co$Nyear)
  if (co$CPUE[Year] > 0)
   { npar <- npar + 1; qbar <- qbar + log(co$CPUE[Year]/Biomass[Year]) }
 qbar <- exp(qbar / npar)

 # Find the SS
 SS <- 0
 for (Year in 1:co$Nyear)
  if (co$CPUE[Year] > 0)
   {
     CpuePred[Year] <- qbar*Biomass[Year]
     SS <- SS + log(CpuePred[Year]/co$CPUE[Year])^2
   }
 # print(SS)

 # Diagnostic stuff
 if (co$Final)
  {
    par(mfrow=c(2,2))
    Years <- seq(1917,1992)
    ymax <- max(co$CPUE,CpuePred)*1.1
    plot(Years[co$CPUE > 0],co$CPUE[co$CPUE > 0],xlab="Year",ylab="CPUE",type='b',pch=16,ylim=c(0,ymax))
    lines(Years[co$CPUE > 0],CpuePred[co$CPUE > 0],lty=2,lwd=2)
    ymax <- max(Biomass)*1.1

    plot(c(Years,1993),Biomass,xlab="Year",ylab="Biomass",type='l',pch=16,lty=1,lwd=2,ylim=c(0,ymax))
    FinalPar <- NULL
    FinalPar$Biomass <- Biomass
    FinalPar$CPUE <- CpuePred
    FinalPar$Pars <- c(r,K,qbar)
   assign("FinalPar",FinalPar,pos=1)
  }

 return(SS)
}

# ========================================================================================================================================

f2 <- function(LogK,logr,Nyear,Catch,CPUE,Final)
{
 K <- exp(LogK)
 r <- exp(logr)

 # Temp storage
 CpuePred <- rep(0,Nyear)
 Biomass <- rep(0,Nyear+1)

 # Do projections
 Biomass[1] <- K
 for (Year in 1:Nyear)
  {
    Biomass[Year+1] <- Biomass[Year] + r*Biomass[Year]*(1.0-Biomass[Year]/K) - Catch[Year]
    if (Biomass[Year+1] < 0.01) Biomass[Year+1] <- 0.01
  }

 # Calculate the ML estimate of q
 qbar <- 0; npar <- 0
 for (Year in 1:Nyear)
  if (CPUE[Year] > 0)
   { npar <- npar + 1; qbar <- qbar + log(CPUE[Year]/Biomass[Year]) }
 qbar <- exp(qbar / npar)

 # Find the SS
 SS <- 0
 for (Year in 1:Nyear)
  if (CPUE[Year] > 0)
   {
     CpuePred[Year] <- qbar*Biomass[Year]
     SS <- SS + log(CpuePred[Year]/CPUE[Year])^2
   }
 # print(SS)

 # Diagnostic stuff
 if (Final)
  {
    par(mfrow=c(2,2))
    Years <- seq(1917,1992)
    ymax <- max(CPUE,CpuePred)*1.1
    plot(Years[CPUE > 0],CPUE[CPUE > 0],xlab="Year",ylab="CPUE",type='b',pch=16,ylim=c(0,ymax))
    lines(Years[CPUE > 0],CpuePred[CPUE > 0],lty=2,lwd=2)
    ymax <- max(Biomass)*1.1

    plot(c(Years,1993),Biomass,xlab="Year",ylab="Biomass",type='l',pch=16,lty=1,lwd=2,ylim=c(0,ymax))
    FinalPar <- NULL
    FinalPar$Biomass <- Biomass
    FinalPar$CPUE <- CpuePred
    FinalPar$Pars <- c(r,K,qbar)
   assign("FinalPar",FinalPar,pos=1)
  }

 return(SS)
}

# ========================================================================================================================================

Case1 <- function()
{
 Nyear <- 76    
 FileN <- "Lect8.txt"
 TheData <- scan(FileN,what=list(NULL,Catch=0,CPUE=0),n=Nyear*3)
 print(TheData) 
    
 co <- list()
 co$Final <- T
 co$Nyear <- Nyear
 co$Catch <- TheData$Catch
 co$CPUEObs <- TheData$CPUE
 assign("co",co,pos=1)
 
 xinit <- c(log(1800),log(0.4))
 result<- optim(xinit,f1,method="L-BFGS-B") 
 print(result)
 co$Final <- T
 assign("co",co,pos=1)
 xfin <- result$par
 f1(xfin)   
 NULL
 print(FinalPar)
}

# ========================================================================================================================================

Case2 <- function()
{
 Nyear <- 76    
 FileN <- "Lect8.txt"
 TheData <- scan(FileN,what=list(NULL,Catch=0,CPUE=0),n=Nyear*3)
 print(TheData) 
    
 co <- list()
 co$Final <- T
 co$Nyear <- Nyear
 co$Catch <- TheData$Catch
 co$CPUEObs <- TheData$CPUE
 assign("co",co,pos=1)
 
 xinit <- c(log(1800),log(0.4))
 result<-Fit(xinit,2,f1,T)
 print(result)
 co$Final <- T

 assign("co",co,pos=1)
 xfin <- result$Zfin
 f1(xfin)   
 NULL
 print(FinalPar)
}
# ------------------------------------------------------------------------------------------------------------
Case3 <- function()
{
 Nyear <- 76    
 FileN <- "Lect8.txt"
 TheData <- scan(FileN,what=list(NULL,Catch=0,CPUE=0),n=Nyear*3)
 print(TheData) 
    
 co <- NULL
 co$Final <- F
 co$Nyear <- Nyear
 co$Catch <- TheData$Catch
 co$CPUEObs <- TheData$CPUE
 assign("co",co,pos=1)
 
 xinit <- c(log(1800),log(0.4))
 result<-Fit2(xinit,2,f1,T)
 print(result)
 co$Final <- T

 assign("co",co,pos=1)
 xfin <- result$Zfin
 f1(xfin)   
 NULL
 print(FinalPar)
}
# ------------------------------------------------------------------------------------------------------------
Case4<-function()
{
 Nyear <- 76    
 FileN <- "Lect8.txt"
 TheData <- scan(FileN,what=list(NULL,Catch=0,CPUE=0),n=Nyear*3)
 print(TheData) 
    
 co <- NULL
 co$Final <- F
 co$Nyear <- Nyear
 co$Catch <- TheData$Catch
 co$CPUEObs <- TheData$CPUE
 assign("co",co,pos=1)

 Npnt <- 50
 Y <- seq(0.1,0.6,length=Npnt)
 X <- seq(1000,3000,length=Npnt)
 print(X)
 print(Y)   
 Z <- matrix(0,Npnt,Npnt)
 for (I in 1:Npnt)
  for (J in 1:Npnt)
   {
     v <- c(log(X[I]),log(Y[J]))
     Z[I,J] <- log(f1(v))
     print(Z[I,J])
   }
 par(mfrow=c(1,1))
 contour(X,Y,Z)
 
}

# ------------------------------------------------------------------------------------------------------------
Case5 <- function()
{
 Nyear <- 76    
 FileN <- "Lect8.txt"
 TheData <- scan(FileN,what=list(NULL,Catch=0,CPUE=0),n=Nyear*3)
 
 result<- mle(f2,start=list(LogK=log(1800),logr=log(0.4)),method="L-BFGS-B",fixed=list(Nyear=Nyear,Catch=TheData$Catch,CPUE=TheData$CPUE,Final=F)) 
 print(summary(result))
 print(coef(result))
 logK <- coef(result)[1]
 logr <- coef(result)[2]

 f2(logK,logr,Nyear=Nyear,Catch=TheData$Catch,CPUE=TheData$CPUE,Final=T)   
 print(FinalPar)
 NULL
}
# ------------------------------------------------------------------------------------------------------------

Case1()
#Lect8a()
Lect8b()
Case2()
