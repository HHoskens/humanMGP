#' Human MGP
#' 
#' Generates associated MGP effects from "runHumanMGP" output 
#' 
#' @param obj runHumanMGP object
#' @param comp Which PLS component to visualize (default = 1)
#' @param type Type of plot (options: "max", "min", "overlay", "heatmap"; default = "none")
#' @param sdy Standard deviation of the yscores. sdy=3 will show the effects of -3sd vs +3sd (default = 3)
#' @param lambda Regularization parameter of the TPS (only valid for sparse landmarks; default = 1e-8)
#' @param ncores Number of CPU cores to be used for parallellization (default = 1)
#' 
#' @return 
#' \item{predmax }{Predicted MGP effect at +sdy standard deviations from the mean}
#' \item{predmin }{Predicted MGP effect at -sdy standard deviations from the mean}
#'
#' @author Hanne Hoskens
#' 
#' @export
plotMGPEffect <- function(obj, comp, type, sdy, lambda, ncores){
  
  #  Required
  if (missing(obj)) {
    stop("Please provide runHumanMGP object")
  }
  
  # Optional
  if (missing(comp)) { comp = 1 }
  if (missing(type)) { type = "none" }
  if (missing(sdy)) { sdy = 3 }
  if (missing(lambda)) { lambda = 1e-08 }
  if (missing(ncores)) { ncores = 1 }
  
  max_comp = dim(obj$GeneLoadings)[2]
  if (length(comp)>1){
    stop("Please only select one PLS component at a time")
  }
  if (comp > max_comp) {
    cat("\033[33m", paste("Maximum number of components is ", max_comp, ". Plotting PLS1 instead.", sep = ""), "\033[0m", "\n")
    comp = 1
  }
  

  mntpath = '/mnt/BHServer4/'
  
  
  # Sparse
  if (obj$ShapeType == "sparse") { 
    mesh = Morpho::file2mesh(paste(mntpath,"FaceBase_3/Data/Images/Atlas/Dense_2k_ears/dense_2k_ears_atlas.ply",sep=""))
    ind_sparse = read.csv(paste(mntpath,'FaceBase_3/Data/Images/Atlas/Dense_2k_ears/sparse_65_ind.txt',sep=""),header=F);
    sparse_lm = mesh$vb[1:3,as.matrix(ind_sparse)]
    
    nlm = dim(mesh$vb)[2]
    nsplm = dim(ind_sparse)[1]
    
    
    avg = colMeans(obj$PLS$y)
    
    plsEffects = plsCoVar(obj$PLS, i=comp, sdy=sdy)
    predminlm = matrix((avg + plsEffects$y[1,]),3,nsplm)
    predmaxlm = matrix((avg + plsEffects$y[2,]),3,nsplm)
    
    # Transform sparse to dense
    out = list()
    out$predmin = tps3d(mesh,t(sparse_lm),t(predminlm),lambda=lambda,threads=ncores)
    out$predmax = tps3d(mesh,t(sparse_lm),t(predmaxlm),lambda=lambda,threads=ncores)
  
    
  # Dense
  } else if (obj$ShapeType == "dense") { 
    mesh = Morpho::file2mesh(paste(mntpath,"FaceBase_3/Data/Images/Atlas/Dense_5k/dense_5k_atlas.ply",sep=""))
    nlm = dim(mesh$vb)[2]
    load(paste(mntpath,"FaceBase_3/Analysis/HumanMGP/Data/",toupper(obj$ShapeType),"PCA_",toupper(obj$Cohort),".RData",sep=""))
    
    avg = row2array3d(pheno.avg, Nlandmarks = nlm)[,,1]
    
    plsEffects = plsCoVar(obj$PLS, i=comp, sdy=sdy)
    predminlm = avg + row2array3d(t(pca.eigvec %*% plsEffects$y[1,]), Nlandmarks = nlm)[,,1]
    predmaxlm = avg + row2array3d(t(pca.eigvec %*% plsEffects$y[2,]), Nlandmarks = nlm)[,,1]
    
    # Save mesh
    out = list()
    out$predmin = mesh; out$predmin$vb[1:3,] = t(predminlm)
    out$predmax = mesh; out$predmax$vb[1:3,] = t(predmaxlm)
  }
  
  
  
  # Generate shape viewer
  if (type == "none"){
    # skip
    
  } else if (type == "max"){
    open3d(zoom = .65); 
    plot3d(out$predmax, col = adjustcolor("darkgrey", .3), alpha = .9, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso", add = F)
    view3d(theta=1, phi=5, fov=45)
    
  } else if (type == "min"){
    open3d(zoom = .65); 
    plot3d(out$predmin, col = adjustcolor("darkgrey", .3), alpha = .9, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso", add = F)
    view3d(theta=1, phi=5, fov=45)
    
  } else if (type == "overlay"){
    open3d(zoom = .65); 
    plot3d(out$predmin, col = "darkgrey", alpha = .7, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso", add = T)
    plot3d(out$predmax, col = "lightgrey", alpha = .7, specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "", main = "", aspect = "iso", add = T, type = "wire")
    view3d(theta=98, phi=-5, fov=45)
    
  } else if (type == "heatmap"){
    avg = mesh; avg$vb = (out$predmin$vb + out$predmax$vb)/2
    
    tmpdist = meshDist(avg, out$predmax, plot = F)
    
    open3d(zoom = .65)
    meshDist(avg, out$predmax, rampcolors = coolwarm(5), plot = T, steps = 100, scaleramp = T, from = -max(abs(tmpdist$dists)), to = max(abs(tmpdist$dists)))
    view3d(theta=1, phi=5, fov=45)
  }
  
  
  return(out)
  
}
