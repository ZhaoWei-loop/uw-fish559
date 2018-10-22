Sign<-function(a,b)
{
 if (b > 0)
  return(abs(a))
 else
  return(-1*abs(a))
}

# =====================================================================================================================================================================

Fit<-function(Xinit,Ndim,FUNK,doGraph)
{
 # Tolerances, etc,
 Ftol = 0.00000001
    
 # Set up a matrix of initial directions    
 XI <- matrix(0,Ndim,Ndim)  
 for (I in 1:Ndim) XI[I,I] = 1
  
 ZZ <- Powell(Xinit,XI,Ndim,FUNK,Ftol)
 print(c("final",ZZ))

 if (doGraph)
  {
   par(mfrow=c(2,2))
   Nfunc <- length(ZZ$TraceY)
   X <- seq(1,Nfunc,1)
   Y <- ZZ$TraceY
   plot(X,Y,xlab="Function call",ylab="Function value",type='b',pch=16,lwd=2)
   X <- NULL; Y <-NULL
   Npnt <- length(ZZ$TraceX)/2
   for (I in 1:Npnt)
    {
     X <- c(X,ZZ$TraceX[(I-1)*2+1])
     Y <- c(Y,ZZ$TraceX[(I-1)*2+2])
    }
    plot(X,Y,xlim=c(-2,3),ylim=c(-2,4),type='n',lwd=2)
   lines(c(0,0),c(1,1),lwd=10,pch=16)
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

 for (Iter in 1:200)
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
       { Del = abs(FPTT-Fret); IBIG = I }
    
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
      T = 2*(FP-2*Fret+FPTT)*(FP-Fret-Del)^2-Del*(FP-FPTT)^2
      if (T < 0)
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
 ZZ$X <- X
 ZZ$Fret <- FUNK(X)
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

 print(c("start",FV))
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
      R <- (X-W)*(FX-FW); Q = (X-V)*(FX-FW); P = (X-V)*Q-(X-W)*R
      Q <- 2.0*(Q-R)
      if (Q > 0) P = -1*P
      Q <- abs(Q)
      Etemp <- E; E <- D
      Ok <- T
      if (abs(P) < abs(0.5*Q*Etemp) && (P > Q*(A-X)) && P <= Q*(B-X))
       {
         D = P/Q; U <- X +D
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
       A = X
      else
       B = X
      V = W; FV <- FW; W <= X; FW <- FX; X <- U; FX <- FU
    
    }
   else
    {
      if (U < X)
       A = U
      else
       B = U
      if (FU <= FW || W == X)
       { V <- W; FV <- FW; W <- U; FW <- FU }
      else
       {
         if (FU <= FV || V == X || V == W)
          { V = U; FV = FU }
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
    print("continue")
   R <- (BX-AX)*(FB-FC)
   Q <- (BX-CX)*(FB-FA)
    U <- BX - ((BX-CX)*Q-(BX-AX)*R)/(2.0*Sign(max(abs(Q-R),Tiny),Q-R))
    Ulim <- BX + Glimit*(CX-BX)
      
   # try various possibilities
    if ((BX-U)*(U-CX) > 0) 
     {
      print("here1")
      FU <- F1DIM(U,Pcom,XIcom,FUNK)
     FVals <- c(FVals,FU)
      if (FU < FC)
       { AX <- BX; FA <- FB; BX <- U; FB <- FU; print("beaking1"); break }
      else
        if (FU > FB)
         { CX <- U; FC <- FU; print("beaking2"); break}
      print("no break")
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
          print("here3")    
          if ((U-Ulim)*(Ulim-CX) >=0)
            { U <- Ulim; FU <- F1DIM(U,Pcom,XIcom,FUNK);    FVecs <- c(FVecs,U); FVals <- c(FVals,FU) }
          else
            { U <- CX + Gold*(CX-BX);   FU <- F1DIM(U,Pcom,XIcom,FUNK); FVals <- c(FVals,FU) }  
         }
     } # if ((BX-U)*(U-CX) > 0) 

     # Elimiate oldest point and continue
     print("elim")
     AX <- BX; BX <- CX; CX <- U; FA <- FB; FB <- FC; FC <- FU
 } # while
  print("done")

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

g <- function(x) { return(1+x[1]*x[1]+(x[2]-1)*(x[2]-1)) }

#g <- function(x) { return(1+x*x) }

test <- function()
{
 X<- c(2,2)     
 Fit(X,2,g,T)   
 
}

test()
