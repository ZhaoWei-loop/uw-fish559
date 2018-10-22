Lect3<-function()
{
 Case1()    
}
Case1<-function()
{
 par(mfrow=c(2,3))
 xlow <- -2
 xhi <- 3

 PX1 <- c(0,0,1,0)
 PY1 <- c(1,3,2,1)
 plot(PX1,PY1,type='b',lwd=2,xlab="parameter 1",ylab="parameter 2",xlim=c(xlow,xhi),ylim=c(-1,5))
 PX2 <- c(0,0,-1,0)
 PY2 <- c(1,3,2,1)
 plot(PX2,PY2,type='b',lwd=2,xlab="parameter 1",ylab="parameter 2",xlim=c(xlow,xhi),ylim=c(-1,5))
 PX3 <- c(0,0,-2,0)
 PY3 <- c(1,3,2,1)
 plot(PX3,PY3,type='b',lwd=2,xlab="parameter 1",ylab="parameter 2",xlim=c(xlow,xhi),ylim=c(-1,5))
 PX4 <- c(0,0,-0.5,0)
 PY4 <- c(1,3,2,1)
 plot(PX4,PY4,type='b',lwd=2,xlab="parameter 1",ylab="parameter 2",xlim=c(xlow,xhi),ylim=c(-1,5))
 PX5 <- c(0,0,0.5,0)
 PY5 <- c(1,2,1.5,1)
 plot(PX5,PY5,type='b',lwd=2,xlab="parameter 1",ylab="parameter 2",xlim=c(xlow,xhi),ylim=c(-1,5))
 }
Lect3()
