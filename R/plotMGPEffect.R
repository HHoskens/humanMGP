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
  
  
  #nlandmarks = 
  nsplandmarks = 65
  
  
  
  nlm = dim(obj$PLS$y)[2]

  # Sparse
  if (nlm == 3*nsplandmarks) { 
    mesh = Morpho::file2mesh("/mnt/BHServer4/FaceBase_3/Data/Images/Atlas/Dense_2k_ears/dense_2k_ears_atlas.ply")
    ind_sparse = read.csv('/mnt/BHServer4/FaceBase_3/Data/Images/Atlas/Dense_2k_ears/sparse_65_ind.txt',header=F);
    sparse_lm = mesh$vb[1:3,as.matrix(ind_sparse)]
    
    avg = colMeans(obj$PLS$y)
    
    plsEffects = plsCoVar(obj$PLS, i=comp, sdy=sdy)
    predminlm = matrix((avg + plsEffects$y[1,]),3,nsplandmarks)
    predmaxlm = matrix((avg + plsEffects$y[2,]),3,nsplandmarks)
    
    # Transform sparse to dense
    out = list()
    out$predmin = tps3d(mesh,t(sparse_lm),t(predminlm),lambda=lambda,threads=ncores)
    out$predmax = tps3d(mesh,t(sparse_lm),t(predmaxlm),lambda=lambda,threads=ncores)
  }
  # Dense
  
  
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
