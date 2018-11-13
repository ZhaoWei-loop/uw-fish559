Lecture19 <- function()
 {
  set.seed(1999)
 # Case1()
  #Case2()
  Case3()
  #Case4()
 }

# ===============================================================

Case1 <- function()
 {
  Ndata <- 10
  Vector <- rbeta(Ndata,5,6)*10
  Vector <- sort(Vector)
  par(mfrow=c(2,2))
  plot(1:Ndata,Vector,type="l",xlab="i",ylab="x")
  abline(h=5,col="red",lwd=4)
  
  #x <- rnorm(1e7)
  #print(system.time(x1 <- sort(x, method = "shell")))
  #print(system.time(x2 <- sort(x, method = "quick")))
  print(Vector[which(Vector==max(Vector[Vector < 5]))])
  
  Vector <- sort(Vector)  
  Lo <- 1
  Hi <- length(Vector)
  cat(Lo,Hi,"\n")
  while (round(Lo) < round(Hi)-1)
   {
    Test <- (Lo+Hi)%/%2
    if (Vector[Test] < 5)
     Lo <- Test
    else 
     Hi <- Test
    cat(Lo,Hi,"\n")
   }  
  print(Vector[Test])
  AAA
 }

# ===============================================================

Case2 <- function()
 {
  Ndata <- 10
  Mat <- matrix(NA,nrow=Ndata,ncol=2) 
  Mat[,1] <- rnorm(Ndata,0,1)
  Mat[,2] <- rnorm(Ndata,0,1)
  Index <- sort.int(Mat[,1],index.return=T)$ix
  Mat <- Mat[Index,]
 }  
  

# ===============================================================

Case3 <- function()
{

 #xx <- sample(1:100000,size=2000,rep=T)
 set.seed(18019)
 xx <- sample(1:20,size=15,rep=F)
 #print(xx)


 BubbleSort <- function(Vec)
  {
   # This is a really dumb way to sort a vector
   n <- length(Vec)
   swapped <- T
   #print(Vec)
   while (swapped)
   {
    swapped <- F  
    for (II in 2:n)  
      if (Vec[II-1] > Vec[II]) 
      { AA <- Vec[II]; Vec[II] <- Vec[II-1]; Vec[II-1] <- AA; swapped <- T }
    #print(Vec)
   }  
   return(Vec)    
  }

 InsertionSort <- function(Vec)
  {
   #print(Vec)
   for (II in 1:length(Vec))
    {
     AA <- Vec[II] 
     JJ <- II
     while (JJ > 1 && Vec[JJ-1] > AA) { Vec[JJ] <- Vec[JJ-1]; JJ <- JJ - 1}
     Vec[JJ] <- AA;
     #print(Vec)
    }  
   return(Vec)
  }  

 Shell <- function(Vec)
  {
   n <- length(Vec) 
   inc <- 1
   while (inc <= n) inc <- inc + 3
   #print(Vec)
   while (inc > 1)
   {
    inc <- inc %/% 3
    for (II in (inc+1):n)
    {
      AA <- Vec[II]
      JJ <- II
      while (JJ > inc && Vec[JJ-inc] > AA) { Vec[JJ] <- Vec[JJ-inc]; JJ <- JJ-inc}
      Vec[JJ] <- AA      
    }  
    #print(Vec)
   }  
   return(Vec)
  }  

 QuickSort <- function(Vec)
  {
   n <- length(xx) 
   ST1 <- rep(NA,n)
   Left <- rep(NA,n)
   Right <- rep(NA,n)
   StkLen <- 1
   Left[1] <- 1
   Right[1] <- n
   
   while(StkLen > 0)
   {
    print(Vec)
    Mid <- Vec[Left[StkLen]]
    Lefts <- Left[StkLen]
    Rights <- Right[StkLen]
    LS <- Left[StkLen]
    RS <- Right[StkLen]
    #cat(Mid,StkLen,LS,RS,"\n")
    
    # Do a one-level sort
    for (IC in c((Left[StkLen]+1):Right[StkLen]))
    {
      # Check whether the current is less than the middle
      if (Vec[IC] > Mid) 
      { ST1[Rights] <- Vec[IC]; Rights <- Rights - 1; }
      else
      { ST1[Lefts] <- Vec[IC]; Lefts <- Lefts + 1; }
    }  
    
    # Store the middle value
    ST1[Lefts] <- Vec[Left[StkLen]]
    
    # Replace the data
    for (IC in c(Left[StkLen]:Right[StkLen])) Vec[IC] <- ST1[IC]  
    StkLen <- StkLen - 1
    
    # update right points
    if ((Lefts-LS) > 1)
    {
      StkLen <- StkLen + 1
      Left[StkLen] <- LS
      Right[StkLen] <- Lefts - 1
    }
    
    # update right points
    if ((RS-Rights) > 1)
    {
      StkLen <- StkLen + 1
      Left[StkLen] <- Rights + 1
      Right[StkLen] <- RS
    }
    
   }  
   print(Vec)
   return(Vec)
 }   

 print(system.time(yy1 <- BubbleSort(xx)))
 print(system.time(yy2 <- InsertionSort(xx)))
 print(system.time(yy3 <- Shell(xx)))
 print(system.time(yy4 <- QuickSort(xx)))
 print(system.time(yy5 <- sort(xx)))
 print(sum(yy1-yy2))
 print(sum(yy1-yy3))
 print(sum(yy1-yy4))
 print(sum(yy1-yy5))

 #QuickSort(xx)
 
}

Case4 <- function()
 {
  lo <- 0
  hi <- 1
 }

Lecture19()