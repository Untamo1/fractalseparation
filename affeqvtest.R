set.seed(12345)

w <- 30
h <- 90
n <- h*w
simn <- n

ts1 <- arima.sim(model=list(ar=c(0.25,0.6)),n=simn) + rt(n,9)
ts2 <- arima.sim(model=list(ar=c(0.15,0.4)),n=simn) + rt(n,9)
ts3 <- arima.sim(model=list(ar=c(0.1,0.2)) ,n=simn) + rt(n,9)

ts1c <- -4.2*ts1  +rt(n,10)
ts2c <-  2.0*ts2  +rbeta(n,0.5,0.5)
ts3c <-  0.9*ts3  +runif(n,0,1) 

Z<- cbind(complex(re=ts1,im=ts1c),complex(re=ts2,im=ts2c),complex(re=ts3,im=ts3c))



##################Real
ini.cov <- cov1(Z)
ini.acov <- acov(Z,4)
ini.2cov <- cov2(Z)

A.re <- c(1,-3,2,-4,9,0,-2,1,1)
A.im <- c(9,-2,1,1.5,-20,-5,6,2,1)
A <- matrix(complex(re=A.re,im=A.im),ncol=3)





trans <- Z %*% t(A) # or tcrossprod



trans.cov <- cov1(trans)
trans.acov <- acov(trans,4)
trans.2cov <-  cov2(trans)


tmpcov <- A %*% ini.cov %*% t(Conj(A))
tmpacov <- A %*% ini.acov %*% t(Conj(A))
tmp2cov <- A %*% ini.2cov %*% t(Conj(A))

trans.cov -tmpcov

norm(tmpcov - trans.cov,type="F")
norm(tmp2cov - trans.2cov,type="F")
norm(tmpacov - trans.acov,type="F")

FOBIori <- NFOBI(Z)
names(FOBIori)
GF <- FOBIori$Gamma

GF %*% cov1(Z) %*% t(Conj(GF))
GF %*% cov2(Z) %*% t(Conj(GF))

cov1(FOBIori$Data)
cov2(FOBIori$Data)

Ftrans <- NFOBI(trans)$Gamma

MD_fun(GF,diag(3))
MD_fun(Ftrans,A)


AMUSEori <- NAMUSE(Z,5)
GA <- AMUSEori$Gamma

GA %*% cov1(Z) %*% t(Conj(GA))
GA %*% acov(Z,5) %*% t(Conj(GA))

cov1(AMUSEori$Data)
acov(AMUSEori$Data,5)

GA2 <- NAMUSE(trans,5)$Gamma

MD_fun(GA,diag(3))
MD_fun(GA2,A)

SOBIori <- NSOBI(Z,1)
GS <- SOBIori$Gamma
GS2 <- NSOBI(trans,1)$Gamma

GS %*% cov1(Z) %*% t(Conj(GS))
GS %*% acov(Z,1) %*% t(Conj(GS))

cov1(SOBIori$Data)
acov(SOBIori$Data,1)

MD_fun(GS,diag(3)) - MD_fun(GS2,A)




library(JADE)
library(clue)
