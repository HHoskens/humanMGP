MGPobj = setClass("MGPobj")

plotMGPEffect = function(object){
  UseMethod("plotMGPEffect")
}

plotMGPEffect.MGPobj <- function(object, comp, sdy, lambda, ncores){
  
  if (missing(comp)) { comp = 1 }
  if (missing(sdy)) { sdy = 3 }
  if (missing(lambda)) { lambda = 1e-08 }
  if (missing(ncores)) { ncores = 0 }
  
  nsplandmarks = 65
  
  nlm = dim(object$PLS$y)[2]

  # Sparse
  if (nlm == 3*nsplandmarks) { 
    mesh = Morpho::file2mesh("/mnt/BHServer4/FaceBase_3/Data/Images/Atlas/Dense_2k_ears/dense_2k_ears_atlas.ply")
    ind_sparse = read.csv('/mnt/BHServer4/FaceBase_3/Data/Images/Atlas/Dense_2k_ears/sparse_65_ind.txt',header=F);
    sparse_lm = mesh$vb[1:3,as.matrix(ind_sparse)]
    
    avg = colMeans(object$PLS$y)
    
    plsEffects = plsCoVar(object$PLS, i=comp, sdy=sdy)
    predminlm = matrix((avg + plsEffects$y[1,]),3,nsplandmarks)
    predmaxlm = matrix((avg + plsEffects$y[2,]),3,nsplandmarks)
    
    # Transform sparse to dense
    predmin = tps3d(mesh,t(sparse_lm),t(predminlm),lambda=lambda,threads=ncores)
    predmax = tps3d(mesh,t(sparse_lm),t(predmaxlm),lambda=lambda,threads=ncores)
  }
  
  return(predmin,predmax)
  
}
