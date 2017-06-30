
# F1,F2,F3 are the images
# OMEGA = mixing matrix
# taus = set of taus for AMUSE/SOBI (separately for AMUSE and use them together for SOBI)
# name = name of the printed images
# checkimage = TRUE/FALSE, if true checks whether any pixels are on the north-pole of the stereographic projection
ImageSep <- function(F1,F2,F3,OMEGA,taus,name,checkimage){
  he <- min(dim(F1)[1],dim(F2)[1],dim(F3)[1])
  wi <- min(dim(F1)[2],dim(F2)[2],dim(F3)[2])

  # Project the colors to the surface of the RBG cube
  F1.corrected <- correct_pixels(F1,he,wi)
  F2.corrected <- correct_pixels(F2,he,wi)
  F3.corrected <- correct_pixels(F3,he,wi)
  
  # Plot the corrected figures
  nameoriginal <- paste("results/original",name,".pdf",sep="")
  
  pdf(nameoriginal,width=10,height=10,paper="special")
  par(mfrow = c(1,3),mar = c(0.1,0.1,0.1,0.1))
  plot_jpeg(F1.corrected)
  plot_jpeg(F2.corrected)
  plot_jpeg(F3.corrected)
  dev.off()
  
  # Transform the cube surface to sphere and
  # stereohraphic transformation from the sphere to the complex plane
  
  # Check wheter there are points on the north-pole
  if(checkimage==TRUE){
    imagecheck(F1.corrected,F2.corrected,F3.corrected,he,wi)
  }
  
  # Transform the images to the complex plane
  F1.complex <- pixtocomp(cube_to_sphere(F1.corrected))
  F2.complex <- pixtocomp(cube_to_sphere(F2.corrected))
  F3.complex <- pixtocomp(cube_to_sphere(F3.corrected))

  # Vectorize the complex valued images
  Z <- cbind(as.vector(F1.complex),as.vector(F2.complex),as.vector(F3.complex))
  
  # Mix the images
  MIXED <- Z %*% t(OMEGA)
  
  # Reverse vectorize the images
  M1 <- matrix(MIXED[,1],he,wi)
  M2 <- matrix(MIXED[,2],he,wi)
  M3 <- matrix(MIXED[,3],he,wi)
  
  #Transform from complex plane to sphere and from sphere to cube surface
  M1.RBG <- sphere_to_cube(comptopix(M1,he,wi))
  M2.RBG <- sphere_to_cube(comptopix(M2,he,wi))
  M3.RBG <- sphere_to_cube(comptopix(M3,he,wi))

  namemixed <- paste("results/mixed",name,".pdf",sep="")
  
  pdf(namemixed,width=10,height=10,paper="special")
  par(mfrow = c(1,3),mar = c(0.1,0.1,0.1,0.1))
  plot_jpeg(M1.RBG)
  plot_jpeg(M2.RBG)
  plot_jpeg(M3.RBG)
  dev.off()
  
  FOB <- NFOBI(MIXED)
  
  FOBI.EST <- FOB$Data
  FOBI.GAMMA <- FOB$Gamma
  
  
  cat(sprintf("The MD index for FOBI is %.8f \n",  MD_fun(FOBI.GAMMA,OMEGA)))
  # cat(sprintf("The MD index for SOBI is %.8f \n",  MD_fun(t(SOBI.GAMMA),OMEGA)))
  
  

  EST1.FOBI <- matrix(FOBI.EST[,1],he,wi)
  EST2.FOBI <- matrix(FOBI.EST[,2],he,wi)
  EST3.FOBI <- matrix(FOBI.EST[,3],he,wi)
  

  
  EST1.FOBI.RBG <- sphere_to_cube(comptopix(EST1.FOBI,he,wi))
  EST2.FOBI.RBG <- sphere_to_cube(comptopix(EST2.FOBI,he,wi))
  EST3.FOBI.RBG <- sphere_to_cube(comptopix(EST3.FOBI,he,wi))

  
  namefobi <- paste("results/fobiunmixed",name,".pdf",sep="")
  # namesobi <- paste("sobiunmixed",name,".jpeg",sep="")
  pdf(namefobi,width=10,height=10,paper="special")
  par(mfrow = c(1,3),mar = c(0.1,0.1,0.1,0.1))
  plot_jpeg(EST1.FOBI.RBG)
  plot_jpeg(EST2.FOBI.RBG)
  plot_jpeg(EST3.FOBI.RBG)
  dev.off()
  
  # jpeg(namesobi,width=480,height=150)
  # par(mfrow = c(1,3),mar = c(0.1,0.1,0.1,0.1))
  # plot_jpeg(EST1.SOBI.RBG)
  # plot_jpeg(EST2.SOBI.RBG)
  # plot_jpeg(EST3.SOBI.RBG)
  # dev.off()
  
  
  n <- length(taus)
  
  for(i in 1:n){
    AMU <- NAMUSE(MIXED,taus[i])
    AMUSE.EST <- AMU$Data
    EST.AMUSE.MAT <- AMU$Gamma
    
    
    cat(sprintf("The MD index for AMUSE is %.8f, with tau %.d \n",  MD_fun(EST.AMUSE.MAT,OMEGA),taus[i]))
 
 
    EST1.AMUSE <- matrix(AMUSE.EST[,1],he,wi)
    EST2.AMUSE <- matrix(AMUSE.EST[,2],he,wi)
    EST3.AMUSE <- matrix(AMUSE.EST[,3],he,wi)

    
    EST1.AMUSE.RBG <- sphere_to_cube(comptopix(EST1.AMUSE,he,wi))
    EST2.AMUSE.RBG <- sphere_to_cube(comptopix(EST2.AMUSE,he,wi))
    EST3.AMUSE.RBG <- sphere_to_cube(comptopix(EST3.AMUSE,he,wi))
    
    picname <- paste("results/amuseunmixed_tau",taus[i],name,".pdf",sep="")
    
    pdf(picname,width=10,height=10,paper="special")
    par(mfrow = c(1,3),mar = c(0.1,0.1,0.1,0.1))
    plot_jpeg(EST1.AMUSE.RBG)
    plot_jpeg(EST2.AMUSE.RBG)
    plot_jpeg(EST3.AMUSE.RBG)
    dev.off()
  }
}

