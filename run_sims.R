source("functions.R")
ncores <- detectCores() - 1

ngenes = 20
nsamples = 20
cor.gene <- 0.7
cor.sample <- 0.0
background.genecor <- 0
background.samplecor <- 0
efron.iter = 10
nsim <- 10
sims <- mclapply(1:nsim, function(ii){
  
  
  genedraw <- sample(1:ngenes, replace = F)
  modules <- list("module A" = genedraw[1:round(ngenes/4)],
                  "module B" = genedraw[(round(ngenes/4)+1):round(2*ngenes/4)],
                  "module C" = genedraw[(round(2*ngenes/4)+1):round(3*ngenes/4)])
  
  sampledraw <- sample(1:nsamples, replace = F)  
  batches <- list("batch A" = sampledraw[1:round(nsamples/4)],
                  "batch B" = sampledraw[(round(nsamples/4)+1):round(2*nsamples/4)],
                  "batch C" = sampledraw[(round(2*nsamples/4)+1):round(3*nsamples/4)])
  
  means.genes <-rnorm(ngenes, mean = 0, sd = 1) # gaussian distribution for now
  means.samples <- rnorm(nsamples, mean = 0, sd = 1) # gaussian distribution for now
  
  # assign samples to treatment groups
  treatments <- list("case" = 1:round(nsamples/2),
                     "control" = (round(nsamples/2)+1):nsamples)
  
  means.treatments.list <- list("case" = 0, "control" = -1)
  means.batches.list <- list("batch A" = 1, "batch B" = 0, "batch C" = 0, "unclustered" = 0)
  means.modules.list <- list("module A" = 1, "module B" = 0, "module C" = 0, "unclustered" = -1)
  means.modules_treatments.list <- list(data.frame(module = "module A",
                                                   treatment = "case",
                                                   shift = 1),
                                        data.frame(module = "module B",
                                                   treatment = "case",
                                                   shift = -1),
                                        data.frame(module = "module C",
                                                   treatment = "case",
                                                   shift = 0))
  
  sd.genes <- rep(1, ngenes)
  sd.samples <- rep(1, nsamples)
  
  cor.intramodules <- list("module A" = cor.gene,
                           "module B" = cor.gene,
                           "module C" = cor.gene,
                           "unclustered" = background.genecor)
  
  
  cor.intrabatch <- list("batch A" = cor.sample,
                         "batch B" = cor.sample,
                         "batch C" = cor.sample,
                         "unclustered" = background.samplecor)
  
  
  w.SX0S <- 1
  w.g <- 1
  w.s <- 1
  w.batch <- 0
  w.treatment <- 0
  w.module <- 0
  w.dex <- 1
  
  
  dat <- simulate_data(ngenes = ngenes,
                       nsamples = nsamples,
                       modules = modules,
                       batches = batches,
                       treatments = treatments,
                       means.genes = means.genes,
                       means.samples = means.samples,
                       means.treatments.list = means.treatments.list,
                       means.batches.list = means.batches.list,
                       means.modules.list = means.modules.list,
                       means.modules_treatments.list = means.modules_treatments.list,
                       cor.intramodules = cor.intramodules,
                       cor.intrabatch = cor.intrabatch,
                       background.genecor = background.genecor,
                       background.samplecor = background.samplecor,
                       sd.genes = sd.genes,
                       sd.samples = sd.samples,
                       w.SX0S = w.SX0S,
                       w.g = w.g,
                       w.s = w.s,
                       w.batch = w.batch,
                       w.treatment = w.treatment,
                       w.module = w.module,
                       w.dex = w.dex,
                       double_std_iter = efron.iter,
                       random.seed = NULL)
  dat
}, mc.cores = ncores)

sims.RAW <- lapply(sims, function(this_sim) this_sim$Xraw)
sims.X0 <- lapply(sims, function(this_sim) this_sim$X0)
sims.data <- lapply(sims, function(this_sim) this_sim$data)
sims.SX0S <- lapply(sims, function(this_sim) this_sim$SX0S)
sims.g <- lapply(sims, function(this_sim) this_sim$z.g)
sims.s <- lapply(sims, function(this_sim) this_sim$z.s)
sims.dex <- lapply(sims, function(this_sim) this_sim$z.dex)






require(ggbiplot)
w.dex <- sims[[1]]$w.dex

# compare single and double standardization
for(ds in c("none", "double")){
  df.X0 <- do.call(rbind, mclapply(sims.X0, function(this_dat){
    A <- double_standardize(this_dat, efron.iter, use = ds)
    sxeg <- prcomp(t(A), center = TRUE, scale. = TRUE) # samples x eigengenes
    gxea <- prcomp(A, center = TRUE, scale. = TRUE) # genes x eigenarrays
    this_df <- data.frame(rho.reads.eigengene = cor(sxeg$x, colMeans(this_dat)),
                          rho.expr.eigenarray = cor(gxea$x, rowMeans(this_dat)),
                          rho.expr.eigengene = cor(sxeg$rotation, rowMeans(this_dat)),
                          rho.reads.eigenarray = cor(gxea$rotation, colMeans(this_dat)))
    this_df$PC <- substr(rownames(this_df), 3, 6)
    this_df$layer <- "A. Random"
    this_df
  }, mc.cores = ncores))
  
  df.SX0S <- do.call(rbind, mclapply(sims.SX0S, function(this_dat){
    A <- double_standardize(this_dat, efron.iter, use = ds)
    sxeg <- prcomp(t(A), center = TRUE, scale. = TRUE) # samples x eigengenes
    gxea <- prcomp(A, center = TRUE, scale. = TRUE) # genes x eigenarrays
    this_df <- data.frame(rho.reads.eigengene = cor(sxeg$x, colMeans(this_dat)),
                          rho.expr.eigenarray = cor(gxea$x, rowMeans(this_dat)),
                          rho.expr.eigengene = cor(sxeg$rotation, rowMeans(this_dat)),
                          rho.reads.eigenarray = cor(gxea$rotation, colMeans(this_dat)))
    this_df$PC <- substr(rownames(this_df), 3, 6)
    this_df$layer <- "B. Only gene and sample correlations"
    this_df
  }, mc.cores = ncores))
  
  
  df.data <- do.call(rbind, mclapply(sims.data, function(this_dat){
    A <- double_standardize(this_dat, efron.iter, use = ds)
    sxeg <- prcomp(t(A), center = TRUE, scale. = TRUE) # samples x eigengenes
    gxea <- prcomp(A, center = TRUE, scale. = TRUE) # genes x eigenarrays
    this_df <- data.frame(rho.reads.eigengene = cor(sxeg$x, colMeans(this_dat)),
                          rho.expr.eigenarray = cor(gxea$x, rowMeans(this_dat)),
                          rho.expr.eigengene = cor(sxeg$rotation, rowMeans(this_dat)),
                          rho.reads.eigenarray = cor(gxea$rotation, colMeans(this_dat)))
    this_df$PC <- substr(rownames(this_df), 3, 6)
    this_df$layer <- "C. Full Model"
    this_df
  }, mc.cores = ncores))
  
 
  
  gg.df <- melt(rbind(df.X0, df.data, df.SX0S))
  pcmax <- 10
  require(ggsci)
  gg <- ggplot(gg.df) + aes(x = as.numeric(PC), y = value^2, color = layer) + 
    #geom_jitter(width = 0.2, alpha = 0.1) + 
    geom_smooth(method = "loess", span = .2) + 
    facet_grid(~variable, 
               labeller = as_labeller(c("rho.reads.eigengene"="read coverage ~ eigengene", 
                                        "rho.expr.eigenarray"="avg expr ~ eigenarray",
                                        "rho.expr.eigengene" = "avg expr ~ eigengene",
                                        "rho.reads.eigenarray" = "read coverage ~ eigenarray"
               ))) + 
    scale_x_continuous(name ="Eigenvector", limits=c(1,pcmax), breaks = seq(1,pcmax,by=1), labels = seq(1,pcmax,by=1)) +
    scale_y_continuous(name ="Variance Explained", limits=c(0,1), breaks = seq(0,1,by=.1), labels = seq(0,1,by=.1)) +
    scale_color_d3(name = "Generative Model") + theme_bw()
  
  gg
  ggsave(plot=gg, filename=paste0("./figures/correlation_and_meanshift_", ngenes, "_genes_", nsamples, "_samples_", w.dex, "_wdex_", cor.gene, "_genecor_", cor.sample, "_samplecor_", background.samplecor, "_backgroundsamplecor_", ds, "_standardization.pdf"), height=3, width = 13)
  
  
  require(mclust)
  method = "ward.D2"
  rands <- mclapply(sims, function(this_sim){
    clust.gene <- hclust(as.dist(1-cor(double_standardize(t(this_sim$data), niter = 100, use = ds))/2), method = method)
    rand.gene <- adjustedRandIndex(as.numeric(cutree(clust.gene, k=4)[this_sim$annotations.genes$gene]), as.numeric(this_sim$annotations.genes$module))
    
    clust.sample <- hclust(as.dist(1-cor(double_standardize(this_sim$data, niter = 100, use = ds))/2), method = method)
    rand.sample <- adjustedRandIndex(as.numeric(cutree(clust.sample, k=4)[this_sim$annotations.samples$sample]), as.numeric(this_sim$annotations.samples$batch))
    
    this_df <- data.frame(module.randindex = rand.gene,
                          batch.randindex = rand.sample)
    this_df
  }, mc.cores = ncores)
  df.rand <- do.call(rbind, rands)
  
  rands.cor <- mclapply(sims, function(this_sim){
    clust.gene <- hclust(as.dist(1-cor(double_standardize(t(this_sim$SX0S), niter = 100, use = ds))/2), method = method)
    rand.gene <- adjustedRandIndex(as.numeric(cutree(clust.gene, k=4)[this_sim$annotations.genes$gene]), as.numeric(this_sim$annotations.genes$module))
    
    clust.sample <- hclust(as.dist(1-cor(double_standardize(this_sim$SX0S, niter = 100, use = ds))/2), method = method)
    rand.sample <- adjustedRandIndex(as.numeric(cutree(clust.sample, k=4)[this_sim$annotations.samples$sample]), as.numeric(this_sim$annotations.samples$batch))
    
    this_df <- data.frame(module.randindex = rand.gene,
                          batch.randindex = rand.sample)
    this_df
  }, mc.cores = ncores)
  df.rand.cor <- do.call(rbind, rands.cor)
  
  rands.null <- mclapply(sims, function(this_sim){
    clust.gene <- hclust(as.dist(1-cor(double_standardize(t(this_sim$X0), niter = 100, use = ds))/2), method = method)
    rand.gene <- adjustedRandIndex(as.numeric(cutree(clust.gene, k=4)[this_sim$annotations.genes$gene]), as.numeric(this_sim$annotations.genes$module))
    
    clust.sample <- hclust(as.dist(1-cor(double_standardize(this_sim$X0, niter = 100, use = ds))/2), method = method)
    rand.sample <- adjustedRandIndex(as.numeric(cutree(clust.sample, k=4)[this_sim$annotations.samples$sample]), as.numeric(this_sim$annotations.samples$batch))
    
    this_df <- data.frame(module.randindex = rand.gene,
                          batch.randindex = rand.sample)
    this_df
  }, mc.cores = ncores)
  df.rand.null <- do.call(rbind, rands.null)
  
  
  df.rand.null$model <- "A. Null"
  df.rand.cor$model <- "B. Only gene and sample correlations"
  df.rand$model = "C. Full model"
  
  
  gg.rand <- ggplot(rbind(df.rand, df.rand.null, df.rand.cor)) + aes(x = module.randindex, y = batch.randindex, color = model) + 
    stat_ellipse(level = 0.68) +
    stat_ellipse(level = 0.95, linetype = 2) +
    #stat_density_2d(aes(fill = ..level.., color = model), geom = "polygon", colour="white", alpha = .1) +
    geom_point(alpha = 0.05) + 
    #geom_point(x = mean(df.rand$module.randindex, trim = .1), y= mean(df.rand$batch.randindex, trim = .1), size = 5, color = "red") +
    geom_hline(yintercept = 0, lty = 2) + 
    geom_vline(xintercept = 0, lty = 2) + 
    geom_abline(slope = 1, intercept = 0, lty = 2) +
    #stat_ellipse(level = .68)  + 
    #stat_ellipse(level = 0.95) + 
    scale_x_continuous(name ="Gene Module Adjusted Rand Index", limits = c(-0.2,1)) +
    scale_y_continuous(name ="Sample Batch Adjusted Rand Index", limits = c(-0.2,1)) +
    theme_bw() + 
    ggtitle("Adjusted Rand Index") +
    scale_color_d3(name = "Generative Model") +
    coord_fixed()
  gg.rand
  
  ggsave(plot=gg.rand, filename=paste0("./figures/rand_index_", ngenes, "_genes_", nsamples, "_samples_", w.dex, "_wdex_", cor.gene, "_genecor_", cor.sample, "_samplecor_", background.samplecor, "_backgroundsamplecor_", ds, "_standardization.pdf"), height=3.5, width = 6)
 
}





require(pheatmap)
require(RColorBrewer)

dat <- sims.data[[1]]
data.standardized <- double_standardize(dat)
A <- data.standardized
pheatmap(A,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         breaks = c(seq(-max(abs(A)), 0 - 0.01, length.out = 50),
                    seq(0, max(abs(A)), length.out = 50)),
         cluster_rows = T,
         cluster_cols = T,
         scale = "none",
         annotation_row = sims[[1]]$annotations.genes,
         annotation_col = sims[[1]]$annotations.samples,
         show_rownames = F,
         show_colnames = F,
         border_color = NA,
         main = "Expression Levels")

B <- cor(A)
pheatmap(B,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         breaks = c(seq(-max(abs(B)), 0 - 0.01, length.out = 50),
                    seq(0, max(abs(B)), length.out = 50)),
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = F,
         show_colnames = F,
         border_color = NA,
         annotation_row = sims[[1]]$annotations.samples,
         annotation_col = sims[[1]]$annotations.samples,
         main = "Sample Correlation")

C <- cor(t(A))
pheatmap(C,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         breaks = c(seq(-max(abs(C)), 0 - 0.01, length.out = 50),
                    seq(0, max(abs(C)), length.out = 50)),
         cluster_rows = T,
         cluster_cols = T,
         show_rownames = F,
         show_colnames = F,
         border_color = NA,
         annotation_row = sims[[1]]$annotations.genes,
         annotation_col = sims[[1]]$annotations.genes,
         main = "Gene Correlation")


