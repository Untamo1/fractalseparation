imagecheck <- function(F1,F2,F3,he,wi){
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  
  for(i in 1:he){
    for(j in 1:wi){
      if( sum(F1[i,j,]==c(0,0,1))==3){
        sum1 <- sum1+1
        }
      if( sum(F2[i,j,]==c(0,0,1))==3){
        sum2 <- sum2+1
        }
      if( sum(F3[i,j,]==c(0,0,1))==3){
        sum3 <- sum3+1
        }
      }
    }
  cat(sprintf("First image has %.f pixels on the north pole \n", sum1 ))
  cat(sprintf("Second image has %.f pixels on the north pole \n", sum2 ))
  cat(sprintf("Third image has %.f pixels on the north pole \n", sum3 ))
}



correct_pixels <- function(X,he,wi){

  
  sh <- max(0,floor(dim(X)[1]/2-he/2)) 
  sw <- max(0,floor(dim(X)[2]/2-wi/2)) 
  
  X.scale.tmp <- X[sh:(sh+he),sw:(sw+wi),]
  X.scale <- X.scale.tmp[1:he,1:wi,]
  
  return(X.scale)
  
}


generateunitaryomega <- function(theta,alpha,beta,gamma){

  
  I <- (cos(gamma) + 2*cos(theta)*cos(alpha))/3
  J <- (cos(gamma) + 2*cos(theta)*cos(alpha+2*pi/3))/3
  K <- (cos(gamma) + 2*cos(theta)*cos(alpha+4*pi/3))/3
  R <- (sin(gamma) + 2*sin(theta)*cos(beta))/3
  G <- (sin(gamma) + 2*sin(theta)*cos(beta+ 2*pi/3))/3
  B <- (sin(gamma) + 2*sin(theta)*cos(beta+4*pi/3))/3
  
  tmp1 <- matrix(c(I,J,K,K,I,J,J,K,I),nrow=3,byrow=T)
  tmp2 <- matrix( c( complex(im=R),complex(im=B),complex(im=G),
                     complex(im=B),complex(im=G),complex(im=R),
                     complex(im=G),complex(im=R),complex(im=B)),nrow=3,
                  byrow=T)
  
  OMEGA <- tmp1 + tmp2
  
  return(OMEGA)
  
}


plot_jpeg = function(jpg, add=FALSE)
{
  res = dim(jpg)[1:2] # get the resolution
  if (!add) # initialize an empty plot area if add==FALSE
    plot(0,0,xlim=c(0,res[1]),ylim=c(0,res[2]),asp=1,type='n',xaxs='i'
         ,yaxs='i',xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  rasterImage(jpg,0,0,res[1],res[2])
}


#Project the pixels to the surface of the RBG cube
cor_pix <- function(A){
  d1 <- dim(A)[1]
  d2 <- dim(A)[2]
  d3 <- dim(A)[3]
  
  for( i in 1:d1){
    for( j in 1:d2){
      tmp_max <- max(A[i,j,])
      tmp_min <- min(A[i,j,])
      diff1 <- 1 - tmp_max
      
      if( diff1 <= tmp_min){
        A[i,j,match(tmp_max,A[i,j,])] <- 1
      }
      else{
        A[i,j,match(tmp_min,A[i,j,])] <- 0
        
        
      }
      
    }
    
  }
  return(A)
  
}


cube_to_sphere <- function(A){
  D <- A
  A2 <- A*2 - 1
  n <- dim(A)[1]
  p <- dim(A)[2]
  for(i in 1:n){
    for(j in 1:p){
      x <- A2[i,j,1]
      y <- A2[i,j,2]
      z <- A2[i,j,3]
      D[i,j,1] <- x*sqrt(1- y^2/2- z^2/2 + (y^2*z^2)/3)
      D[i,j,2] <- y*sqrt(1- z^2/2- x^2/2 + (x^2*z^2)/3)
      D[i,j,3] <- z*sqrt(1- x^2/2- y^2/2 + (x^2*y^2)/3)
    }
  }
  return(D)
}


sphere_to_cube <- function(A){
 
  n <- dim(A)[1]
  p <- dim(A)[2]
  D <- array(dim = c(n,p,3))
  for(i in 1:n){
    for(j in 1:p){
      x <- A[i,j,1]
      y <- A[i,j,2]
      z <- A[i,j,3]
      
      fx = abs(x)
      fy = abs(y)
      fz = abs(z)
      if (fy >= fx && fy >= fz) {
        a2 = x * x * 2.0
        b2 = z * z * 2.0
        inner = -a2 + b2 -3
        innersqrt = -sqrt((inner * inner) - 12.0 * a2)
        
        if(x == 0.0 || x == -0.0) { 
          position.x = 0.0
        }else {
          position.x = sqrt(innersqrt + a2 - b2 + 3.0) * 1/sqrt(2)
        }
        
        if(z == 0.0 || z == -0.0) {
          position.z = 0.0
        }else {
          position.z = sqrt(innersqrt - a2 + b2 + 3.0) * 1/sqrt(2)
        }
        
        if(position.x > 1.0) position.x = 1.0
        if(position.z > 1.0) position.z = 1.0
        
        if(x < 0) position.x = -position.x
        if(z < 0) position.z = -position.z
        
        if (y > 0) {
          #top face
          position.y = 1.0;
        }else {
          # bottom face
          position.y = -1.0;
        }
      }else if (fx >= fy && fx >= fz) {
        a2 = y * y * 2.0;
        b2 = z * z * 2.0;
        inner = -a2 + b2 -3;
        innersqrt = -sqrt((inner * inner) - 12.0 * a2);
        
        if(y == 0.0 || y == -0.0) { 
          position.y = 0.0; 
        }else {
          position.y = sqrt(innersqrt + a2 - b2 + 3.0) * 1/sqrt(2)
        }
        
        if(z == 0.0 || z == -0.0) {
          position.z = 0.0;
        }else {
          position.z = sqrt(innersqrt - a2 + b2 + 3.0) * 1/sqrt(2)
        }
        
        if(position.y > 1.0) position.y = 1.0;
        if(position.z > 1.0) position.z = 1.0;
        
        if(y < 0) position.y = -position.y;
        if(z < 0) position.z = -position.z;
        
        if (x > 0) {
          # right face
          position.x = 1.0;
        }else {
          # left face
          position.x = -1.0;
        }
      }else {
        a2 = x * x * 2.0;
        b2 = y * y * 2.0;
        inner = -a2 + b2 -3;
        innersqrt = -sqrt((inner * inner) - 12.0 * a2);
        
        if(x == 0.0 || x == -0.0) { 
          position.x = 0.0; 
        }else {
          position.x = sqrt(innersqrt + a2 - b2 + 3.0) * 1/sqrt(2)
        }
        
        if(y == 0.0 || y == -0.0) {
          position.y = 0.0;
        }else {
          position.y = sqrt(innersqrt - a2 + b2 + 3.0) * 1/sqrt(2)
        }
        
        if(position.x > 1.0) position.x = 1.0;
        if(position.y > 1.0) position.y = 1.0;
        
        if(x < 0) position.x = -position.x;
        if(y < 0) position.y = -position.y;
        
        if (z > 0) {
          # front face
          position.z = 1.0;
        }else {
          # back face
          position.z = -1.0;
        }
      }
      D[i,j,1] <- position.x
      D[i,j,2] <- position.y
      D[i,j,3] <- position.z
      
    }
  }
  D.scaled <- 1/2*(1+D)
  return(D.scaled)
}





comptopix <- function(A,he,wi){
  height <- he
  width  <- wi
  TM <- array(dim=c(he,wi,3))
  for(i in 1:height){

    for(j in 1:width){

      a <- Re(A[i,j])
      b <- Im(A[i,j])
      TM[i,j,1] <- (2*a)/(a^2+b^2+1)
      TM[i,j,2] <- (2*b)/(a^2+b^2+1)
      TM[i,j,3] <- (a^2+b^2-1)/(a^2+b^2+1)
    }
  }
  return(TM)
}


pixtocomp <- function(D){
  #Scale and shift the location to fit the transformation
  A <- D
  height <- dim(A)[1]
  width  <- dim(A)[2]
  TM <- matrix(NA,nrow=height,ncol=width)
  for(i in 1:height){
    for(j in 1:width){
      repart <- A[i,j,1]/(1-A[i,j,3])
      impart <- A[i,j,2]/(1-A[i,j,3])
      TM[i,j] <- complex(re=repart,im=impart)
    }
  }
  return(TM)
}





