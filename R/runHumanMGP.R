#' Human MGP
#' 
#' Performs a two-block PLS on sets of genetic and shape variables based on "pls2B" from the Morpho package.
#' The goal is to find pairs of latent variables that show maximal covariation with a set of genetic markers.
#' 
#' Note: mount Storage4 as "/mnt/BHServer4/"
#' 
#' @param GOterm GO term (name or ID) or gene list (names or ensembl ID) of interest
#' @param cohort Cohort to be analyzed ("3DFN" or "TANZ")
#' @param lm Sparse or dense landmarking scheme to be used ("sparse" or "dense"; default = sparse)
#' @param covs List of covariates to standardize for (options: "none", "age", "age2", "age3", "sex", "height", "weight", "csize")
#' @param window Search window (in kb) around transcription start and end sites (default = 0kb)
#' @param ncomp Number of PLS components to consider (default = 1)
#' @param npc Number of PCs used to build gene composite score (1 or "all"; default = 1)
#' @param signif Logical value indicating whether or not to perform significance testing (default = F)
#' @param nperm Number of permutations for significance testing (default = 99)
#' @param ncores Number of CPU cores to be used for parallellization (default = 1)
#' 
#' @return 
#' \item{GOterm }{Search term}
#' \item{cohort }{Cohort}
#' \item{PLS }{Summary of pls2B output}
#' \item{GeneLoadings }{Matrix (nxncomp) containing the relative loadings per gene composite score}
#' \item{ShapeEffectMax }{Matrix (3xnLMxncomp) containing the shape effects associated with each PLS component}
#' \item{ShapeEffectMin }{Matrix (3xnLMxncomp) containing the "inverse" shape effects associated with each PLS component}
#' \item{ShapeAverage }{Matrix (3xnLM) containing the mean shape of the specified cohort}
#' \item{ShapeType }{sparse or dense landmarking scheme; required for plotting functions}
#' \item{SV }{Singular values (covariance) per pair of PLS axes; plus significance}
#' \item{R2 }{Shape variance explained per PLS component; plus significance}
#' \item{PD }{Procrustes Distance between mean shape and MGP effect; plus significance}
#'
#' @author Hanne Hoskens
#' 
#' @export
runHumanMGP <- function(GOterm,cohort,lm,covs,window,ncomp,npc,signif,nperm,ncores){
  
  ## SPECIFY INPUT VARIABLES ####
  #  Required
  if (missing(GOterm)) {
    stop("Please provide GO term or gene list")
  }
  
  if (missing(cohort)) {
    stop("Please specify cohort (3DFN or TANZ)")
  }
  cohort = toupper (cohort)
  if (cohort != "3DFN" & cohort != "TANZ") {
    stop("Please specify valid cohort (3DFN or TANZ)")
  }

  #  Optional
  if (missing(lm)) { lm = "sparse" }
  if (missing(covs)) { covs = "none" }
  if (missing(window)) { window = 0e3 }
  if (missing(ncomp)) { ncomp = 1 }
  if (missing(npc)) { npc = 1 }
  if (missing(signif)) { signif = F }
  if (missing(nperm)) { nperm = 99 }
  if (missing(ncores)) { ncores = 1 }

  

  
  ## LOAD DATA ####
  cat("\033[32m", "STEP 1: LOADING PHENOTYPES", "\033[0m", "\n")
  mntpath = "/mnt/BHServer4/"
  loadpath1 = paste(mntpath,"FaceBase_3/Analysis/HumanMGP/Data/",sep="")
  loadpath2 = paste(mntpath,"FaceBase_3/Data/Genetics/",sep="")
  
  # Phenotypes + covariates
  # sparse: pheno.id, pheno.coeff, pheno.avg, meta.cov
  # dense: pheno.id, pheno.coeff, pheno.avg, meta.cov, pca.eigvec, pca.eigstd
  load(paste(loadpath1,toupper(lm),"DATA_",cohort,".RData",sep=""))

  #if (lm == "sparse"){
    #tmp = read.csv(paste(loadpath1,"LM_",cohort,"_",lm,".csv",sep=""), header = F, sep = ",")
    #pheno.id = as.character(tmp[,1])
    #pheno.coeff = as.matrix(tmp[,2:dim(tmp)[2]]); colnames(tmp) <- NULL
    #pheno.avg = colMeans(pheno.coeff); pheno.avg = as.matrix(t(pheno.avg))
    #meta.cov = read.csv(paste(loadpath1,"COV_",cohort,"_",lm,".csv",sep=""), header = T, sep = ",")
  #} else if (lm == "dense"){
    #tmp = read.csv(paste(loadpath1,"PC_",cohort,"_",tolower(lm),".csv",sep=""), header = F, sep = ",")
    #pheno.id = as.character(tmp[,1])
    #pheno.coeff = as.matrix(tmp[,2:dim(tmp)[2]]); colnames(tmp) <- NULL
    #pca.eigvec = read.csv(paste(loadpath1,"EIGVEC_",cohort,"_",tolower(lm),".csv",sep=""), header = F, sep = ","); pca.eigvec = as.matrix(pca.eigvec)
    #pca.eigstd = read.csv(paste(loadpath1,"EIGVAL_",cohort,"_",tolower(lm),".csv",sep = ""), header = F, sep = ","); pca.eigstd = as.matrix(pca.eigstd)
    #pheno.avg = read.csv(paste(loadpath1,"AVG_",cohort,"_",tolower(lm),".csv",sep = ""), header = F, sep = ","); pheno.avg = as.matrix(pheno.avg)
    #meta.cov = read.csv(paste(loadpath1,"COV_",cohort,"_",lm,".csv",sep=""), header = T, sep = ",")
  #}

  nlandmarks = 5629
  if (lm == "sparse"){ nsplandmarks = 65; nlandmarks = 2565 }
  
  # Genotypes
  if (cohort == "3DFN"){ 
    bfile = paste(loadpath2,"Pitt/QC/Prune/Marazita_imputed_qc_prune_rmrel",sep="") 
  } else if (cohort == "TANZ"){ 
    bfile = paste(loadpath2,"Tanz/QC/Prune/Spritz_imputed_qc_prune_rmrel",sep="") 
  }
  geno.full.bim = readBIM(bfile)
  geno.full.fam =  readFAM(bfile)
  GRCh = 37
  
  
  
  
  ## COVARIATE ADJUSTMENT ####
  cat("\033[32m", "STEP 2: COVARIATE STANDARDIZATION", "\033[0m", "\n")
  #  Only keep individuals that have covariate info
  sort_idx1 = pheno.id %in% meta.cov$IID
  pheno.id = pheno.id[sort_idx1]
  pheno.coeff = pheno.coeff[sort_idx1,]
  #  Make sure shape data and meta data are in the same order
  sort_idx2 = match(pheno.id,meta.cov$IID)
  meta.cov = meta.cov[sort_idx2,]
  
  # Standardization
  if ("none" %in% covs){
    # No standardization
    pheno.coeff.adj = pheno.coeff
    pheno.avg.adj = pheno.avg
    
  } else {
    cov_idx = match(tolower(covs),tolower(colnames(meta.cov)))
    cov_keep = meta.cov[,cov_idx,drop=F]
    
    # Remove individuals with missing data
    keep_id = rowSums(is.na(cov_keep)) < 1
    cov_keep = cov_keep[keep_id,,drop=F]
    pheno.coeff = pheno.coeff[keep_id,,drop=F]
    pheno.id = pheno.id[keep_id]
    
    # Make geomorph dataframe
    df = geomorph.data.frame(cov_keep)
    df$Shape = pheno.coeff
    
    # Fit model
    f1 = as.formula(paste(names(df)[length(df)],(paste(names(df)[1:(length(df)-1)], collapse = " + ")), sep = " ~ "))
    mod = procD.lm(f1, data = df, RRPP = T, iter = 99, Parallel = ncores)

    # Get residuals
    #pheno.coeff.adj = mod$residuals + matrix(data=pheno.avg,nrow=nrow(pheno.coeff),ncol=ncol(pheno.coeff),byrow = T)
    pheno.coeff.adj = mod$residuals + matrix(data=colMeans(df$Shape),nrow=nrow(pheno.coeff),ncol=ncol(pheno.coeff),byrow = T)
    pheno.avg.adj = colMeans(pheno.coeff.adj)
  }

  
  
  

  ## GENE/GO SEARCH TERM ####
  cat("\033[32m", "STEP 3: MATCH GENES TO PROVIDED SEARCH TERM", "\033[0m", "\n")
  # List all IDs of GO terms available in org.Hs object, and reduce to biological processes
  GO_ID = toTable(org.Hs.egGO)
  GO_ID = unique(GO_ID$go_id[GO_ID$Ontology == "BP"])
  # Match GO IDs with process names
  GO_all = toTable(GOTERM)
  GO_BP = GO_all[GO_all$Ontology == "BP",]
  GO_BP = GO_BP[match(GO_ID,GO_BP$go_id),]
  
  # List all gene names available in org.Hs object
  GENE_all = toTable(org.Hs.egSYMBOL)
  # Match ensembl IDs with gene names
  ENSG_all = toTable(org.Hs.egENSEMBL)
  ind = match(ENSG_all$gene_id,GENE_all$gene_id)
  GENE_all$ensembl_id = ""; GENE_all$ensembl_id[ind] = ENSG_all$ensembl_id
  
  # Biological process (NAME / ID)
  if (sum(tolower(GOterm) %in% tolower(GO_BP$Term) | tolower(GOterm) %in% tolower(GO_BP$go_id))) {
    ind1 = match(tolower(GOterm),tolower(GO_BP$Term))
    ind2 = match(tolower(GOterm),tolower(GO_BP$go_id))
    ind = c(ind1,ind2); ind = ind[!is.na(ind)]
    
    go_id = GO_BP$go_id[ind]
    go2symbol = unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = go_id, columns = c("ENSEMBL", "SYMBOL"), keytype = "GO")[,-2:-3]))
    
  # Gene (NAME / ID)
  } else if (sum(toupper(GOterm) %in% toupper(GENE_all$symbol) | toupper(GOterm) %in% toupper(GENE_all$ensembl_id))) {   
    ind1 = match(toupper(GOterm),toupper(GENE_all$symbol))
    ind2 = match(toupper(GOterm),toupper(GENE_all$ensembl_id))
    ind = ind1; ind[is.na(ind1)] = ind2[is.na(ind1)]
    
    go2symbol = unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = GENE_all$symbol[ind[!is.na(ind)]], columns = c("ENSEMBL", "SYMBOL"), keytype = "SYMBOL")))
    
    # Non coding [to be added] 
    #gene.list2 <- chdr.list[grepl("TCONS_", chdr.list$gene_id2, fixed = TRUE),]
    #db2 <- TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts
    #symbol2info2 <- transcriptsBy(db2, by=c("gene"), use.names=FALSE)[gene.list2]
    #symbol2info2 <- symbol2info2@unlistData
    
    # List genes that were not found in list
    if (sum(is.na(ind))>0) { cat("\033[33m", paste("No match was found for ", gsub(" ",", ",paste(GOterm[is.na(ind)],collapse = " ")), sep=""), "\033[0m", "\n") }
    
    
  # No match is found  
  } else {
    stop("No match was found for GO term")
  }
  
  
  
  
  ## GENES TX START-END SITE ####   
  cat("\033[32m", "STEP 4A: MAP SNPs TO GENE LIST", "\033[0m", "\n")
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=GRCh)
  symbol2info = getBM(attributes = c('ensembl_gene_id','chromosome_name','transcript_start', 'transcript_end', 'hgnc_symbol'), 
                      mart = ensembl, 
                      filters = 'ensembl_gene_id', 
                      values = go2symbol$ENSEMBL)
  
  # Replace empty gene names
  tmp.ind = symbol2info$hgnc_symbol==""
  tmp.symbol = go2symbol$SYMBOL[match(symbol2info$ensembl_gene_id[tmp.ind],go2symbol$ENSEMBL)]
  symbol2info$hgnc_symbol[tmp.ind] = tmp.symbol 
    
  # remove weird chr names
  chr = symbol2info$chromosome_name
  chr = gsub("HSCHR","",chr);  chr = sapply(strsplit(chr,"_"),"[[",1)
  symbol2info$chromosome_name = chr
  if (length(grep(symbol2info$chromosome_name, pattern = "HG")) > 0) { symbol2info = symbol2info[-grep(symbol2info$chromosome_name, pattern = "HG"),] }
  
  # Rename sex chromosomes
  symbol2info$chromosome_name = gsub("X","23",symbol2info$chromosome_name)
  symbol2info$chromosome_name = gsub("Y","24",symbol2info$chromosome_name)
  
  # Skip if no match is found
  if (dim(symbol2info)[1] == 0) { stop("No genes could be mapped to GO search term") }
  
  # Select maximum transcript size per gene
  n = length(unique(symbol2info$hgnc_symbol))
  
  genvar.name = rep(NA, n)
  genvar.start = rep(NA,  n)
  genvar.end = rep(NA,  n)
  chr.name = rep(NA,  n)
  for(i in 1:n){
    ind = symbol2info$hgnc_symbol == unique(symbol2info$hgnc_symbol)[i]

    genvar.name[i] = unique(symbol2info$hgnc_symbol[ind])
    genvar.start[i] = min(symbol2info$transcript_start[ind])
    genvar.end[i] = max(symbol2info$transcript_end[ind])
    chr.name[i] = unique(symbol2info$chromosome_name[ind])
  }
  seq.info = data.frame(mgi_symbol = genvar.name, chromosome_name = chr.name, start_position = genvar.start, end_position = genvar.end)
  
  
  

  ## FIND MARKERS WITHIN WINDOW ####
  seq.indexes = matrix(nrow = 0, ncol = 3)
  for (i in 1 : dim(seq.info)[1]) {
    tmp.markers = geno.full.bim[which(geno.full.bim$chr == as.numeric(seq.info$chromosome_name[i])),]
    tmp.dist1 = (tmp.markers$bps) > (as.numeric(seq.info$start_position[i]) - window)
    tmp.dist2 = (tmp.markers$bps) < (as.numeric(seq.info$end_position[i]) + window)
    ind = which(tmp.dist1 & tmp.dist2)
    
    if (length(ind) > 0) { seq.indexes = rbind(seq.indexes,as.matrix(cbind(seq.info$mgi_symbol[i],tmp.markers$vid[ind],tmp.markers$bps[ind]))) }
  }
  # Make dataframe
  seq.indexes = data.frame(gene = seq.indexes[,1], snp = seq.indexes[,2], pos = seq.indexes[,3])

  # Skip if no markers are found within window
  if (dim(seq.indexes)[1] == 0) { stop("No SNPs could be mapped to GO search term") }
  
  
  
  
  ## REDUCE GENOTYPES IN PLINK #### 
  plink_path = paste(mntpath,"FaceBase_3/Analysis/HumanMGP/Plink/plink",sep="")
  tmp_path = paste(path.expand("~"),"tmp_mgp/",sep = "/") 
  if (!dir.exists(tmp_path)){ dir.create(tmp_path) }
  
  # Temporarily export SNPs to keep
  rnd_i = sample(1:1e6,1)
  tmp = seq.indexes$snp
  write.csv(tmp, file = paste(tmp_path,"tmp_snp_list_",rnd_i,".txt",sep=""),row.names = F, quote = F)

  # Run PLINK
  # Extract SNPs + remove related individuals
  system2(plink_path, args = paste(" --bfile ", bfile, " --threads ", ncores, " --extract ", tmp_path, "tmp_snp_list_", rnd_i, ".txt --make-bed --out ", tmp_path, "tmp_geno_", rnd_i, " > NUL", sep=""))

  # Import geno files
  geno.bim = readBIM(paste(tmp_path,"tmp_geno_",rnd_i,sep=""))
  geno.fam = readFAM(paste(tmp_path,"tmp_geno_",rnd_i,sep=""))
  geno.bed = readBED(paste(tmp_path,"tmp_geno_",rnd_i,sep="")); geno.bed = as.data.frame(geno.bed)
  
  # Remove tmp files
  unlink(paste(tmp_path,"tmp_snp_list_",rnd_i,".txt",sep=""))
  unlink(paste(tmp_path,"tmp_geno_",rnd_i,".bed",sep=""))
  unlink(paste(tmp_path,"tmp_geno_",rnd_i,".bim",sep=""))
  unlink(paste(tmp_path,"tmp_geno_",rnd_i,".fam",sep=""))
  unlink(paste(tmp_path,"tmp_geno_",rnd_i,".log",sep=""))
  
  
  
  
  ## CHECK ORDER OF GENO AND PHENO ####
  sort.ind1 = seq.indexes$snp %in% geno.bim$vid
  seq.indexes = seq.indexes[sort.ind1,]
  
  sort.ind2 = match(seq.indexes$snp,geno.bim$vid)
  geno.bim = geno.bim[sort.ind2,,drop=F]
  geno.bed = geno.bed[,sort.ind2,drop=F]
  
  sort.ind3 = pheno.id %in% geno.fam$iid
  pheno.id = pheno.id[sort.ind3]
  pheno.coeff.adj = pheno.coeff.adj[sort.ind3,]

  sort.ind4 = match(pheno.id,geno.fam$iid)
  geno.fam = geno.fam[sort.ind4,,drop=F]
  geno.bed = geno.bed[sort.ind4,,drop=F]
  
  
  
  
  ## REPLACE MISSING GENO WITH MAJOR GENOTYPE #### 
  geno.0 = colSums(geno.bed==0,na.rm = T)
  geno.1 = colSums(geno.bed==1,na.rm = T)
  geno.2 = colSums(geno.bed==2,na.rm = T)
  tmp.geno = rbind(geno.0,geno.1,geno.2)
 
  max.ind = apply(tmp.geno, 2, function(x) which.max(x))
  max.val = gsub("geno.","",rownames(tmp.geno)[max.ind])
  max.val = matrix(data=as.numeric(max.val),nrow=dim(geno.bed)[1],ncol=dim(geno.bed)[2],byrow = T)
  
  geno.bed[is.na(geno.bed)] = max.val[is.na(geno.bed)]
  

  

  ## GET GENE COMPOSITE SCORE ####
  cat("\033[32m", "STEP 4B: GET GENE COMPOSITE SCORE", "\033[0m", "\n")
  
  #  Do PCA per gene
  gene.names <- unique(seq.indexes$gene)
  n_gen <- length(gene.names)
  
  gene.score = matrix(NA,dim(pheno.coeff.adj)[1],1e4)
  gene.id = matrix(NA,1,1e4)
  gene.n = matrix(NA,1,1e4)
  count = 0
  for (i in 1:n_gen) {
    ind <- seq.indexes$gene %in% gene.names[i]
    
    # PCA
    tmp.pca <- fast.prcomp(geno.bed[,ind])

    pc_axis = 1
    n = 1
    if (npc == "all"){
      # Parallel analysis (only needed if more than 1 SNP is defined within gene)
      if (sum(ind)>1){
        tmp.pca.pa = paran(geno.bed[,ind], quietly = T, status = F, all = T, iterations = 1000, centile = 95)
        pc_axis = which(tmp.pca.pa$Ev>tmp.pca.pa$RndEv)
        n = length(pc_axis)
        if (n==0){ pc_axis = 1; n = 1 }
      }
    }
  
    gene.score[,(count+1):(count+n)] = as.matrix(tmp.pca$x[,pc_axis])
    tmp.gene.id = as.data.frame(cbind(rep(gene.names[i],1,n),as.character(1:n)))
    gene.id[(count+1):(count+n)] = apply(tmp.gene.id[ , 1:2 ], 1, paste, collapse = "_" )
    gene.n[(count+1):(count+n)] = matrix(n,1,n)
    count = count + n
  }
  gene.score = gene.score[,1:count]
  gene.id = gene.id[1:count]
  gene.n = gene.n[1:count]
  
  
  
  

  ## TWO BLOCK PLS ####
  cat("\033[32m", "STEP 5: TWO-BLOCK PLS", "\033[0m", "\n")
  
  X = as.matrix(gene.score)
  colnames(X) = gene.id
  Y = as.matrix(pheno.coeff.adj)
  
  # Two-block PLS
  mgp.pls <- pls2B(X, Y, rounds = nperm, mc.cores = ncores, useCor = T, cv = F, same.config = F)
  
  # Reduce number of components to be computed if specified number (ncomp) exceeds max possible components 
  max_comp = min(dim(X)[2],dim(Y)[2])
  if (ncomp > max_comp) { 
    ncomp = max_comp
    cat("\033[33m", paste("Maximum number of PLS components is ", ncomp, sep=""), "\033[0m", "\n")
  }
  
  
  
  
  ## SIGNIFICANCE TESTING ####
  cat("\033[32m", "STEP 6: SIGNIFICANCE TESTING", "\033[0m", "\n")
  
  # SV
  mgp.sv = mgp.pls$CoVar$`singular value`[1:ncomp]
    
  # R2
  mgp.r2 = matrix(data=NA, nrow=1, ncol = ncomp)
  mgp.pls$ylm$coefficients <- as.matrix(mgp.pls$ylm$coefficients)
  for (i in 1:ncomp){
    # Adapted from Morpho::predictPLSfromData()
    xs = mgp.pls$Xscores[,i,drop=F] 
    yest <- t(t(mgp.pls$ylm$coefficients[i, , drop = FALSE]) %*% t(xs)) 
      
    ypred <- t(mgp.pls$svd$v %*% t(yest))
    ypred <- sweep(ypred, 2, -mgp.pls$ycenter)
    mgp.pred <- as.matrix(ypred)
      
    ymean = matrix(colMeans(Y),nrow = dim(Y)[1],ncol=dim(Y)[2],byrow = T)
      
    ess = sum((mgp.pred - ymean)^2)
    tss = sum((Y - ymean)^2)
    rss = sum((Y - mgp.pred)^2)
    
    mgp.r2[i] = ess/tss
  }

  # ProcDist
  mgp.pdist = matrix(data=NA, nrow=1, ncol = ncomp)
  for (i in 1:ncomp){
    mgp.predmax = mgp.pls$svd$v[,i] * apply(Y, 2, sd)
    mgp.pdist[i] <- sqrt(sum((mgp.predmax - colMeans(Y))^2))
  }
    
  if (signif==T){
    registerDoParallel(cores=ncores)
    
    # Permutation
    out = foreach(j = 1:nperm, .combine = 'rbind', .packages = c("Morpho")) %dopar% {
      Yperm = Y[sample(1:nrow(Y), size = nrow(Y)),]
      perm.pls <- pls2B(X, Yperm, rounds = 1, mc.cores = ncores, useCor = T, cv = F, same.config = F)
      
      # SV
      forperm.sv = perm.pls$CoVar$`singular value`[1:ncomp]
      
      forperm.r2 = matrix(data=NA, nrow=1, ncol = ncomp);
      forperm.pdist = matrix(data=NA, nrow=1, ncol = ncomp);
      for (i in 1:ncomp){
        # R2
        perm.pls$ylm$coefficients <- as.matrix(perm.pls$ylm$coefficients)
        xs = perm.pls$Xscores[,i,drop=F]
        yest <- t(t(perm.pls$ylm$coefficients[i, , drop = FALSE]) %*% t(xs)) 
        
        ypred <- t(perm.pls$svd$v %*% t(yest))
        ypred <- sweep(ypred, 2, -perm.pls$ycenter)
        mgp.pred <- as.matrix(ypred)
        
        ymean = matrix(colMeans(Y),nrow = dim(Y)[1],ncol=dim(Y)[2],byrow = T)
        
        ess = sum((mgp.pred - ymean)^2)
        tss = sum((Y - ymean)^2)
        rss = sum((Y - mgp.pred)^2)
        
        forperm.r2[i] = ess/tss
        
        # ProcDist
        perm.predmax = perm.pls$svd$v[,i] * apply(Yperm, 2, sd)
        forperm.pdist[i] <- sqrt(sum((perm.predmax - colMeans(Y))^2))
      }
      return(list(forperm.sv,forperm.r2,forperm.pdist))
    }
    stopImplicitCluster()
    
    perm.sv = matrix(unlist(out[,1],use.names = F),nrow = nperm, ncol = ncomp, byrow = T)
    perm.r2 = matrix(unlist(out[,2],use.names = F),nrow = nperm, ncol = ncomp, byrow = T)
    perm.pdist =matrix(unlist(out[,3],use.names = F),nrow = nperm, ncol = ncomp, byrow = T)
    
    mgp.sv.p <- (colSums(matrix(data=mgp.sv,nrow = nperm, ncol = ncomp, byrow = T) <= perm.sv) + 1)/(nperm + 1)
    mgp.r2.p <- (colSums(matrix(data=mgp.r2,nrow = nperm, ncol = ncomp, byrow = T) <= perm.r2) + 1)/(nperm + 1)
    mgp.pdist.p <- (colSums(matrix(data=mgp.pdist,nrow = nperm, ncol = ncomp, byrow = T) <= perm.pdist) + 1)/(nperm + 1)
    
  }
  
  
  
  
  ## EXPORT RESULTS ####
  cat("\033[32m", "STEP 7: SAVE RESULTS", "\033[0m", "\n")
  
  out <- list()
  out$GOterm = GOterm
  out$Cohort = cohort

  # PLS summary
  out$PLS = mgp.pls

  # Geno
  XL = as.matrix(mgp.pls$svd$u[,1:ncomp])
  rownames(XL) = gene.id
  colnames(XL) = paste("PLS",1:ncomp,sep="")
  out$GeneLoadings = XL
  
  # Pheno
  predmin = array(data=NA, dim = c(3,nlandmarks,ncomp))
  predmax = array(data=NA, dim = c(3,nlandmarks,ncomp))
  predavg = array(data=NA, dim = c(3,nlandmarks))
  for (i in 1:ncomp){
    plsEffects = plsCoVar(mgp.pls, i=i, sdy=6)
    if (lm == "sparse"){
      mesh = Morpho::file2mesh(paste(mntpath,"FaceBase_3/Data/Images/Atlas/Dense_2k_ears/dense_2k_ears_atlas.ply",sep=""))
      ind_sparse = read.csv(paste(mntpath,"FaceBase_3/Data/Images/Atlas/Dense_2k_ears/sparse_65_ind.txt",sep=""),header=F);
      sparse_lm = mesh$vb[1:3,as.matrix(ind_sparse)]
      
      predminlm = matrix((pheno.avg.adj + plsEffects$y[1,]),3,nsplandmarks)
      predmaxlm = matrix((pheno.avg.adj + plsEffects$y[2,]),3,nsplandmarks)
      predmin[,,i] = tps3d(mesh,t(sparse_lm),t(predminlm),threads=ncores)$vb[1:3,]
      predmax[,,i] = tps3d(mesh,t(sparse_lm),t(predmaxlm),threads=ncores)$vb[1:3,]
      if (i==1){ predavg = tps3d(mesh,t(sparse_lm),t(matrix(pheno.avg.adj,3,nsplandmarks)),threads=ncores)$vb[1:3,] }
      
    } else if (lm == "dense"){
      mesh = Morpho::file2mesh(paste(mntpath,"FaceBase_3/Data/Images/Atlas/Dense_5k/dense_5k_atlas.ply",sep=""))
      avg = t(row2array3d(pheno.avg,Nlandmarks = nlandmarks)[,,1])
      predmin[,,i] = avg + t(row2array3d(t(pca.eigvec %*% plsEffects$y[1,]), Nlandmarks = nlandmarks)[,,1])
      predmax[,,i] = avg + t(row2array3d(t(pca.eigvec %*% plsEffects$y[2,]), Nlandmarks = nlandmarks)[,,1])
      if (i==1) { predavg = avg }
    }
  }
  dimnames(predmin)[[1]] <- c("x","y","z");   dimnames(predmin)[[2]] <- paste("LM",1:nlandmarks, sep = "");   dimnames(predmin)[[3]] <- paste("PLS",1:ncomp, sep = "")
  dimnames(predmax)[[1]] <- c("x","y","z");   dimnames(predmax)[[2]] <- paste("LM",1:nlandmarks, sep = "");   dimnames(predmax)[[3]] <- paste("PLS",1:ncomp, sep = "")
  dimnames(predavg)[[1]] <- c("x","y","z");   dimnames(predavg)[[2]] <- paste("LM",1:nlandmarks, sep = "")
  
  out$ShapeEffectMax = predmax
  out$ShapeEffectMin = predmin
  out$ShapeAverage = predavg
  out$ShapeType = lm
  
  # Measures of effect size
  mgp.sv = matrix(mgp.sv,nrow=1,ncol=ncomp,byrow=T)
  mgp.sv.perc = mgp.pls$CoVar$`% total covar.`[1:ncomp]
  colnames(mgp.sv) <- paste("PLS",1:ncomp, sep = "")
  colnames(mgp.r2) <- paste("PLS",1:ncomp, sep = "")
  colnames(mgp.pdist) <- paste("PLS",1:ncomp, sep = "")
  
  if (signif == T){
    mgp.sv = rbind(mgp.sv,mgp.sv.perc,mgp.sv.p); rownames(mgp.sv) = c("SV","%Expl Cov","pval (permute shapes)");  
    mgp.r2 = rbind((100*mgp.r2),mgp.r2.p); rownames(mgp.r2) = c("R2","pval (permute shapes)");
    mgp.pdist = rbind(mgp.pdist,mgp.pdist.p); rownames(mgp.pdist) = c("PD","pval (permute shapes)"); 
  }
  out$SV = mgp.sv
  out$R2 = mgp.r2
  out$PD = mgp.pdist

  
  cat("\033[32m", "FINISHED", "\033[0m", "\n")
  
  return(out)

}
