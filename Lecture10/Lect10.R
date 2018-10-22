x <- seq(-1,1,0.1)
D(expression(tan(x)),"x")
D(D(expression(tan(x)),"x"),"x")

dtanx.dx <- function(x)
 {
 gradient <- 1/cos(x)^2
 return(gradient)
 }
cbind(x,dtanx.dx(x))

dtanx2.dx <- D(expression(tan(x)),"x")
eval(dtanx2.dx)

dtanx3.dx <- deriv(expression(tan(x)),"x")
eval(dtanx3.dx)

dtanxy.dxy <- deriv(expression(tan(x*y)),c("x","y"),hessian=TRUE)
y <- seq(-1,1,0.1)
eval(dtanxy.dxy)

dfx.dx <- function(f,x,delta)
{
 y1 <- f(x-delta/2)
 y2 <- f(x+delta/2)
 approx.gradient <- (y2-y1)/delta
 return(approx.gradient)
}

dfx.dx( tan,x,0.001)

indef.tanx.sx <- function(x)
{
 answer <- -log(cos(x))
 return(answer)
}

tan2 <- function(x)
 {
  y <- tan(x)
  print(x)
  print(y)
  return(y)
  
 }

integrate(tan2, 0, 0.5)

library(MASS)
area(tan,0,0.5)


library(cubature)

tan.xy <- function(xy)
 {
  answer <- tan(xy[1] * xy[2])
  return(answer)
 }

adaptIntegrate(tan.xy,lowerLimit=c(0,0),upperLimit=c(0.5,0.5))
