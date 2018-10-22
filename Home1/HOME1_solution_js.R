Homework1 <-function()
{
  library(stats4)
  
  #Ex1()  
  
  Ex2()  
}

# =======================================================================================================================================================================
setwd("Home1")

Ex1<-function()
{
  Streams <- read.table("HOME1A.TXT",sep="") 
  Streams <- matrix(Streams[,3],ncol=3,byrow=T)
  
  low <- -10
  hi <- 10
  Step <- 0.025
  Nstep <- (hi-low)/Step+1
  x <- rep(0,Nstep)
  Weights <- c(1,rep(c(4,2),(Nstep)/2-1),4,1)/3.0*Step
  b <- low + Step*seq(from=0,to=(Nstep-1))

    Verbose <- F
  
  MinFunc <- function(Beta,SigmaB,Sigma,Streams)
  {
    
    #  Find dimensions   
    NStream <- length(Streams[,1])
    Nobs <- length(Streams[1,])
    
    #  Random effect component (note that this is N(0,1) variable)
    yy <- dnorm(b,0,1) 
    
    #  Negative LogLikelihood
    NLL <- 0
    
    #  Work across streams
    for (Istream in 1:NStream)
    {
      
      #    Vec across values of bi
      Vec <- rep(1,length=Nstep)
      for (Iobs in 1:Nobs)
      {
        xx <- dnorm(Streams[Istream,Iobs],(Beta+b*SigmaB),Sigma)
        Vec <- Vec*xx
      }
      Vec <- Vec*yy
      Integral1 <- sum(Vec*Weights)
      NLL <- NLL - log(Integral1)
    }
    if (Verbose) cat(NLL,Vars,"\n")
    return(NLL)   
  }
  
  Vars <- c(70,10,6)
  Test <- MinFunc(Vars[1],Vars[2],Vars[3],Streams)
  print(Test)
  opt <- mle(MinFunc,start=list(Beta=70,SigmaB=10,Sigma=6),fixed=list(Streams=Streams),method="BFGS")
  print("Solution")
  print(summary(opt))
  MLE <- -1*logLik(opt)[1]
  print(MLE)
  coef(opt)
  #AAA
  
  # Check solution
  #Verbose <- T
  #Vars <- coef(opt)
  #MinFunc(Vars[1],Vars[2],Vars[3],Streams)
  
  # Do profile
  Betas <- seq(from=55,to=95,by=2)
  Betas <- c(Betas,seq(from=63,to=65,by=0.01),seq(from=91,to=93,by=0.01))
  Betas <- sort(Betas)
  Verbose <- F
  Likes <- NULL
  for (BetaFix in Betas)
  {
    opt1 <- mle(MinFunc,start=list(SigmaB=10,Sigma=6),fixed=list(Beta=BetaFix,Streams=Streams),method="Nelder-Mead")
    #   cat(BetaFix,-1*logLik(opt1)[1],"\n")
    Likes <- c(Likes,-1*logLik(opt1)[1])
  }
  par(mfrow=c(2,2))
  plot(Betas,Likes-MLE,xlab="Beta",ylab="-Log Likelihood",lty=1,type='l')
  abline(h=1.92) 
  print(cbind(Betas,Likes-MLE)) 
  
  # Fits streams (use ML not REML!)
  library(nlme)
  TheData <- scan("D:\\Courses\\FISH 559_18\\HOMES\\AHOME1\\HOME1A.TXT",what=list(Stream=0,Tst=0,Density=0),n=3*18) 
  Streams <- groupedData(Density~1 | Stream,data=as.data.frame(TheData))
  lm1 <- lme(fixed = Density ~ 1,data=Streams,random = ~ 1 | Stream,method="ML") 
  print(summary(lm1))
  
  NULL
}

# =======================================================================================================================================================================

Ex2<-function()
{
  
  # Read in the basic data and transform the data as needed
  TheData<-scan("Home1b.txt",what=list(SSB=0,Rec=0,Spec2=0,SBPR=0),skip=1,n=4*294)   
  TheData$Spec2 <- factor(TheData$Spec2)
  TheData$REC <- log(TheData$Rec*TheData$SBPR/TheData$SSB)
  print(TheData[1:5])
  
  # Set up a grouped data class
  library(nlme)
  SBPR2 <- groupedData(REC~SSB|Spec2,data=as.data.frame(TheData),FUN=mean)
  
  # Fit the fixed effects model
  print("\nFixed effects models\n")
  fm1 <- lm(REC~SSB*Spec2-1-SSB,data=SBPR2)
  
  print(" ")
  print("Mixed effects models")
  fm2 <- lme(REC~SSB:Spec2,data=SBPR2,random=~1|Spec2,method="REML")
  #fm2 <- lme(REC~SSB:Spec2,data=SBPR2,random=~1|Spec2,method="ML")
  print(summary(fm2))
  #AAA
  
  library(lattice)
  plot(fm2,REC~fitted(.)|Spec2,pch=20,ylab="Observed",abline=c(a=0,b=1))
  plot(fm2,form=resid(.,type="p")~fitted(.)|Spec2,abline=0,pch=20,ylab="Residual")
  #print(residuals(fm2))
  
  par(mfrow=c(2,3))
  
  # QQ plot
  qqnorm(fm2$residuals/0.8122,ylab="Quantiles of Residuals")
  qqline(fm2$residuals/0.8122)
  
  plot(residuals(fm2),ylab="Residuals")
  
  hist(residuals(fm2),xlab="Residuals",main="")
  
  #homogeneity of withing group variance
  boxplot(split(fm2$residuals/0.8122,SBPR2$Spec2),ylab="Residual",xlab="Species",csi=0.2)
  abline(0,0,lwd=3)
  
  #normality of the between-group residuals
  print(fm2$coefficients$random$Spec2)
  re<-fm2$coefficients$random$Spec2/0.6184
  qqnorm(re,ylab="Quantiles of random effects")
  qqline(re)
  
  hist(fm2$coefficients$random$Spec2/0.6184,xlab="random effects",main="")
  
  NULL
}

Homework1()