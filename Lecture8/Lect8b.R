# ==============================================================================================================================================================================

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
   plot(X,Y,xlab="Function call",ylab="Function value",type='b',pch=16,lwd=2)
   X <- NULL; Y <-NULL
   for (I in 1:ZZ$Nfunc)
    {
     X <- c(X,ZZ$TraceX[(I-1)*2+1])
     Y <- c(Y,ZZ$TraceX[(I-1)*2+2])
    }
    plot(X,Y,xlim=c(-2,3),ylim=c(-2,4),type='n',lwd=2)
   lines(c(0,0),c(1,1),lwd=10,pch=16,csi=10)
   lines(X,Y,col=5,lty=1,lwd=2,type='b',pch=16)
   TT <- seq(0,2*pi,length=100)
   for (I in 1:3)
     {
      XX <- I*sin(TT)
      YY <- 1+I*cos(TT)
      lines(XX,YY,lty=1,lwd=2)
     } 
 }

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
    # Plot the simplex?
    if (doPlot)
     {
      XX <- c(P[1,1],P[2,1],P[3,1],P[1,1])
      YY <- c(P[1,2],P[2,2],P[3,2],P[1,2])
      plot(XX,YY,xlim=c(-1,3),ylim=c(-1,3),type='n',lwd=2)
      polygon(XX,YY,col=5)
      lines(c(0,0),c(1,1),lwd=10,pch=16,csi=10)
      TT <- seq(0,2*pi,length=100)
      for (I in 1:3)
       {
        XX <- I*sin(TT)
        YY <- 1+I*cos(TT)
        lines(XX,YY,lty=1,lwd=2)
       } 
     }
    
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
    YPR = FUNK(PR) 
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
              { PR[J] = 0.5*(P[I,J]+P[ILO,J]); 
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

# ==============================================================================================================================================================================
g <- function(x) { return(1+x[1]*x[1]+(x[2]-1)*(x[2]-1)) }

test <- function()
{
 X<- c(2,2)     
 Fit(X,2,g,T)   

}

test()
