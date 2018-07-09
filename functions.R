require(parallel, quietly = T)
require(reshape2, quietly = T)

#_______________________________________________________________________________________________

double_standardize <- function(X, niter = 10, use = "double"){
  if(use == "double"){
    for(i in 1:niter){
      XT <- t(scale(X))
      X <- t(scale(XT))
    }
  } else if(use == "single") {
    X <- scale(X)
  } else {
    X <- X
  }
  return(X[1:nrow(X), 1:ncol(X)])
}

make_random_array <- function(ngenes, nsamples, random.seed = NULL){
  if(!is.null(random.seed)) set.seed(random.seed)
  X0 <- t(mvrnorm(n = nsamples,
                  mu = rep(0, ngenes),
                  Sigma = diag(rep(1, ngenes))))
  rownames(X0) <- paste0("gene_", 1:ngenes)
  colnames(X0) <- paste0("sample_", 1:nsamples)
  return(X0)
}

#_______________________________________________________________________________________________

simulate_data <- function(ngenes,
                          nsamples,
                          modules,
                          batches,
                          treatments,
                          means.genes,
                          means.samples,
                          means.treatments.list, 
                          means.batches.list,
                          means.modules.list,
                          means.modules_treatments.list,
                          cor.intramodules,
                          cor.intrabatch,
                          background.genecor = 0,
                          background.samplecor = 0,
                          sd.genes,
                          sd.samples,
                          w.SX0S,
                          w.g,
                          w.s,
                          w.batch,
                          w.treatment,
                          w.module,
                          w.dex,
                          double_std_iter = 300,
                          random.seed = NULL){
  
  require(reshape2, quietly = T)
  require(MASS, quietly = T)
  
  # generate random matrix ~ N(0x0,1x1)
  Xraw <-make_random_array(ngenes, nsamples, random.seed = random.seed)
  X0 <- double_standardize(Xraw, niter = 300)
  #X0 <- Xraw
  
  #__________________________________________________________________________________________________________________
  # housekeeping
  
  # gather any genes not assigned to a module into an unclustered module
  module.leftovers <- which(!1:ngenes %in% unlist(modules))
  if(length(module.leftovers > 0)) modules[["unclustered"]] <- which(!1:ngenes %in% unlist(modules))
  
  # gather any samples not assigned to a batch into an unclustered batch
  batch.leftovers <- which(!1:nsamples %in% unlist(batches))
  if(length(batch.leftovers) > 0) batches[["unclustered"]] <- which(!1:nsamples %in% unlist(batches))
  
  #__________________________________________________________________________________________________________________
  # mean effects
  
  # genewise and samplewise mean expression
  means.outergenes <- outer(means.genes, rep(1, length(means.samples)))
  if(sd(unlist(means.outergenes)) > 0){
    z.outergenes <- (means.outergenes-mean(means.outergenes))/sd(unlist(means.outergenes))
  } else {
    z.outergenes <- matrix(0, nrow = ngenes, ncol = nsamples)
  }
  
  # genewise and samplewise mean expression
  means.outersamples <- outer(rep(1, length(means.genes)), means.samples)
  if(sd(unlist(means.outersamples)) > 0){
    z.outersamples <- (means.outersamples-mean(means.outersamples))/sd(unlist(means.outersamples))
  } else {
    z.outersamples <- matrix(0, nrow = ngenes, ncol = nsamples)
  }
  
  # treatment group means
  means.treatment <- matrix(0, nrow = ngenes, ncol = nsamples)
  for(ii in 1:length(means.treatments.list)){
    this_treatment_name <- names(means.treatments.list)[ii]
    this_treatment <- treatments[[this_treatment_name]]
    means.treatment[,this_treatment] <- means.treatment[,this_treatment] + means.treatments.list[[this_treatment_name]]
  }
  if(sd(means.treatment) > 0){
    z.treatment <- (means.treatment - mean(means.treatment))/sd(means.treatment)
  } else {
    z.treatment <- matrix(0, nrow = ngenes, ncol = nsamples)
  }
  
  
  # batch means
  means.batch <- matrix(0, nrow = ngenes, ncol = nsamples)
  for(ii in 1:length(means.batches.list)){
    this_batch_name <- names(means.batches.list)[ii]
    this_batch <- batches[[this_batch_name]]
    means.batch[,this_batch] <- means.batch[,this_batch] + means.batches.list[[this_batch_name]]
  }
  if(sd(means.batch) > 0){
    z.batch <- (means.batch - mean(means.batch))/sd(means.batch)
  } else {
    z.batch <- matrix(0, nrow = ngenes, ncol = nsamples)
  }
  
  
  # adjust module means
  means.module <- matrix(0, nrow = ngenes, ncol = nsamples)
  for(ii in 1:length(means.modules.list)){
    this_module_name <- names(means.modules.list)[ii]
    this_module <- modules[[this_module_name]]
    means.module[this_module,] <- means.module[this_module,] + means.modules.list[[this_module_name]]
  }
  if(sd(means.module) > 0){
    z.module <- (means.module - mean(means.module))/sd(means.module)
  } else {
    z.module <- matrix(0, nrow = ngenes, ncol = nsamples)
  }
  
  
  # adjust means.modules_groups (DEX)
  means.module_treatment <- matrix(0, nrow = ngenes, ncol = nsamples)
  for(ii in 1:length(means.modules_treatments.list)){
    this_module_name <- means.modules_treatments.list[[ii]][["module"]]
    this_treatment_name <-  means.modules_treatments.list[[ii]][["treatment"]]
    this_shift <- means.modules_treatments.list[[ii]][["shift"]]
    this_module <- modules[[as.character(this_module_name)]]
    this_treatment <- treatments[[as.character(this_treatment_name)]]
    means.module_treatment[this_module, this_treatment] <- means.module_treatment[this_module, this_treatment] + this_shift
  }
  if(sd(means.module_treatment) > 0){
    z.module_treatment <- (means.module_treatment - mean(means.module_treatment))/sd(means.module_treatment)
  } else {
    z.module_treatment <- matrix(0, nrow = ngenes, ncol = nsamples)
  }
  
  
  #__________________________________________________________________________________________________________________
  # covariance effects
  
  # target covariance matrices
  Sigma.genes.diag <- diag(sd.genes^2) # diagonal
  Sigma.samples.diag <- diag(sd.samples^2) # diagonal
  
  cov.genes <- outer(sd.genes, sd.genes) # off-diagonal
  cov.samples <- outer(sd.samples, sd.samples) # off-diagonal
  
  if(background.genecor != 0){
    Sigma.genes <- Sigma.genes.diag + background.genecor*cov.genes
  } else {
    Sigma.genes <- Sigma.genes.diag
  }
  for(ii in names(cor.intramodules)){
    this_module <- modules[[ii]]
    rho <- cor.intramodules[[ii]]
    edgecor <- matrix(rho, nrow = length(this_module), ncol = length(this_module)) - diag(rep(rho, length(this_module)))
    Sigma.genes[this_module, this_module] <- Sigma.genes.diag[this_module, this_module] + edgecor*cov.genes[this_module, this_module]
  }
  
  if(background.samplecor != 0){
    Sigma.samples <- Sigma.samples.diag + background.samplecor*cov.samples
  } else {
    Sigma.samples <- Sigma.samples.diag
  }
  
  for(ii in names(cor.intrabatch)){
    this_batch <- batches[[ii]]
    rho <- cor.intrabatch[[ii]]
    edgecor <- matrix(rho, nrow = length(this_batch), ncol = length(this_batch)) - diag(rep(rho, length(this_batch)))
    Sigma.samples[this_batch, this_batch] <- Sigma.samples.diag[this_batch, this_batch] + edgecor*cov.samples[this_batch, this_batch]
  }
  
  # Cholesky decomposition
  S.genes <- t(chol(Sigma.genes))
  S.samples <- t(chol(Sigma.samples))
  Xgs <- S.genes %*% X0 %*% t(S.samples)
  
  
  #__________________________________________________________________________________________________________________
  # annotations
  
  # gene annotations
  module.gene.df <- data.frame(module = unlist(lapply(1:length(modules),
                                                      function(i) rep(names(modules)[i], length(modules[[i]])))),
                               gene = as.numeric(unlist(modules)),
                               row.names = paste0("gene_",as.numeric(unlist(modules))))
  annotations.genes <- module.gene.df
  
  
  # sample annotations
  batch.sample.df <- data.frame(batch = unlist(lapply(1:length(batches),
                                                      function(i) rep(names(batches)[i], length(batches[[i]])))),
                                sample = as.numeric(unlist(batches)),
                                row.names = paste0("batch_",as.numeric(unlist(batches))))
  
  treatment.sample.df <- data.frame(treatment = unlist(lapply(1:length(treatments),
                                                              function(i) rep(names(treatments)[i], length(treatments[[i]])))),
                                    sample = as.numeric(unlist(treatments)),
                                    row.names = paste0("treatment_",as.numeric(unlist(treatments))))
  
  merged.sample.df <- merge(batch.sample.df, treatment.sample.df)
  annotations.samples <- merged.sample.df[order(merged.sample.df$sample),]
  rownames(annotations.samples) <- colnames(X0)
  
  
  #__________________________________________________________________________________________________________________
  # sum contributions
  
  rownames(Xgs) <- rownames(z.batch) <- rownames(z.treatment) <- rownames(z.module) <- rownames(z.module_treatment) <- rownames(z.outergenes) <- rownames(z.outersamples) <- rownames(X0)
  colnames(Xgs) <- colnames(z.batch) <- colnames(z.treatment) <- colnames(z.module) <- colnames(z.module_treatment) <- colnames(z.outergenes) <- colnames(z.outersamples) <- colnames(X0)
  Xfinal <- Xgs*w.SX0S + z.outergenes*w.g + z.outersamples*w.s + z.batch*w.batch + z.treatment*w.treatment + z.module*w.module + z.module_treatment*w.dex
  
  ret <- list(data = Xfinal,
              
              Xraw = Xraw,
              X0 = X0,
              
              means.genes = means.genes,
              means.samples = means.samples,
              means.treatments.list = means.treatments.list,
              means.batches.list = means.batches.list,
              means.modules.list = means.modules.list,
              means.modules_treatments.list = means.modules_treatments.list,
              
              sd.genes = sd.genes,
              sd.samples = sd.samples,
              Sigma.genes = Sigma.genes,
              Sigma.samples = Sigma.samples,
              S.genes = S.genes,
              S.samples = S.samples,
              
              cor.intramodules = cor.intramodules,
              cor.intrabatch = cor.intrabatch,
              background.genecor = background.genecor,
              background.samplecor = background.samplecor,
              
              SX0S = Xgs,
              z.g = z.outergenes,
              z.s = z.outersamples,
              z.batch = z.batch,
              z.treatment = z.treatment,
              z.module = z.module,
              z.dex = z.module_treatment,
              
              w.SX0S = w.SX0S,
              w.g = w.g,
              w.s = w.s,
              w.batch = w.batch,
              w.treatment = w.treatment,
              w.module = w.module,
              w.dex = w.dex,
              
              modules = modules,
              batches = batches,
              treatments = treatments,
              
              annotations.genes = annotations.genes,
              annotations.samples = annotations.samples
  )
  
  return(ret)
}



