rm(list=ls())

library(OCNet)
library(rivnet)
library(LaplacesDemon)
library(Rcpp)
library(BayesianTools)
library(eDITH)

if (packageVersion("eDITH") <= "0.3.0"){
  stop("please install eDITH v1.0.0 or higher")
}


# perform simulations with Bayesian method ####
if (!file.exists("results/results_multiplicity.rda")){
  job <- 2 # set to 1 or 2 to partition simulation between 2 cores
  nuFluvial <- 0.5 
  vec_sim <- ((job-1)*30+1):(job*30)
  for (ind_sim in vec_sim){
    set.seed(ind_sim)
    k_mean <- runif(1,-1,0)
    k_cov <- runif(1,0.25,1)
    scaleDistance <- runif(1,5000,20000)
    tau_sol <- runif(1, 1, 24)
    ind_OCN <- sample(50,1)
    omega_k <- runif(1,1,10)  # overdispersion for kicknet
    Cstar <- 5*10^runif(1,0,2)   # eDNA non-detection parameter (5--500, log-spaced)
    sdlog_e <- runif(1,0.1,1) # sd of log(C)
    
    load(paste0("OCN/OCN",ind_OCN,".rda"))
    eval(parse(text=paste0("OCN <- OCN",ind_OCN)))
    
    OCN_S4 <- new("river")
    fieldnames <- names(OCN)
    for (i in 1:length(fieldnames)){slot(OCN_S4, fieldnames[i]) <- OCN[[fieldnames[i]]]}
    OCN <- OCN_S4
    rm(OCN_S4)
    eval(parse(text=paste0("rm(OCN",ind_OCN,")")))
    
    OCN$cellsize <- 100
    OCN$FD$X <- OCN$FD$X*100
    OCN$FD$Y <- OCN$FD$Y*100
    OCN$FD$A <- OCN$FD$A*100^2 # max A: 400 km2
    OCN <- landscape_OCN(OCN)
    OCN <- aggregate_OCN(OCN, thrA=5e5, maxReachLength=1000, equalizeLengths=T) # ~ 600 nodes
    nNodes <- OCN$AG$nNodes
    OCN <- paths_OCN(OCN, level="AG", includeUnconnectedPaths = T)
    
    hd <- data.frame(data=c(10, 1, 1), # max discharge of 10 m3/s
                     type=c("w",  "v", "d"),
                     node=OCN$AG$outlet*c(1,1,1))
    OCN <- hydro_river(hd, OCN)
    
    distMatFluvial <- OCN$AG$downstreamPathLength + t(OCN$AG$downstreamPathLength) +
      OCN$AG$downstreamLengthUnconnected + t(OCN$AG$downstreamLengthUnconnected)
    distMatLand <- as.matrix(dist(cbind(OCN$SC$X,OCN$SC$Y)))
    distMat <- nuFluvial * distMatFluvial + (1-nuFluvial)*distMatLand
    
    covMat <- exp(-distMat/scaleDistance)
    p_sol <- as.numeric(exp(rmvn(1, rep(k_mean,nNodes), Sigma=k_cov*covMat)))
    
    beta <- rep(1, nNodes)
    param <- c(tau_sol, 0, rep(1, nNodes))
    covariates <- data.frame(diag(log(p_sol)))
    names(covariates) <- paste0("beta_",1:nNodes)
    names(param) <-  c("tau", "log_p0", names(covariates))
    
    conc <- eDITH:::run_eDITH_single(param, OCN, covariates, Z.normalize=FALSE) 
    C_sol <- conc$C
    
    # we generate C_obs 3 times, and only pick the first set if multiplicity=1
    C_obs <- rep(C_sol, 3)
    C_obs[runif(length(C_obs)) < exp(-C_obs/Cstar)] <- 0 # non detection probability
    C_obs <- exp(rnorm(length(C_obs),log(C_obs), sdlog_e))
    
    dens_obs <- rnbinom(nNodes, size=p_sol/(omega_k-1), prob=1/omega_k)
    
    nTotSites <- 40
    nSites_e <- 40
    nSites_k <- 0
    ind_str <- 2
    for (multiplicity in c(1,3)){
      fnam <- paste0("results/multiplicity_BT/out",ind_sim,"_s",ind_str,
                     "_e",nSites_e,"_k",nSites_k,"_m",multiplicity,".rda")
      if (!(file.exists(fnam))){
        out <- NULL
        save(out, file=fnam)
        message(sprintf('i: %d  -  sim: %d  -  eSites: %d  -  kSites: %d  -  multiplicity: %d \n',
                        ind_sim,ind_str,nSites_e,nSites_k,multiplicity),appendLF=F)
        
        sampling.sites_e <- samplingStrategy_eDNA(OCN, nSites_e)
        sampling.sites_k <- samplingStrategy_kicknet(OCN, nSites_k)
        
        if (multiplicity==3){
          ind.sites_e <- c(sampling.sites_e, nNodes + sampling.sites_e, 2*nNodes + sampling.sites_e)
          ss_e <- rep(sampling.sites_e, 3)
        } else {
          ind.sites_e <- ss_e <- sampling.sites_e}
        
        data <- data.frame(ID=c(ss_e, sampling.sites_k),
                           values=c(C_obs[ind.sites_e], dens_obs[sampling.sites_k]),
                           type=c(rep("e",length(ss_e)), rep("k",nSites_k)))
        
        set.seed(ind_sim) 
        
        out <- run_eDITH_mixed_BT_nb(data, OCN, ll.type="lnorm", no.det=T,
                                     par.AEM = list(weight="exponential"),
                                     tau.prior = list(spec="unif",min=1, max=24), # same boundaries as tau_sol
                                     log_p0.prior = list(spec="unif",min=-5, max=5),
                                     n.AEM = 20, 
                                     verbose=T) # Bayesian version
        
        out$nSites_e <- nSites_e
        out$nSites_k <- nSites_k
        out$p_sol <- p_sol
        out$C_obs <- C_obs
        out$C_sol <- C_sol
        out$k_mean <- k_mean
        out$k_cov <- k_cov
        out$scaleDist <- scaleDistance
        out$tau_sol <- tau_sol
        out$ind_OCN <- ind_OCN
        out$nNodes <- nNodes
        out$ind_str <- ind_str
        out$omega_k <- omega_k
        out$Cstar <- Cstar
        out$sdlog_e <- sdlog_e
        
        out$covariates <- NULL # save space
        out$source.area <- NULL
        out$probDet <- NULL 
        out$outMCMC <- NULL
        save(out,file=fnam)
      } else {
        message("   Already done! \n", appendLF=F)
      }
      message("\n",appendLF=F)
    }
  }
}

# read simulation output ####
if (!file.exists("results/results_multiplicity.rda")){ 
  k <- 0
  ID <- nNodes <- scaleDist <- ind_OCN <- numeric(120)
  nSites_e <- nSites_k <- onTarget_e <- onTarget_k <- sim <- numeric(120)
  log_p0 <- tau <- tau_sol <- k_mean <- k_cov <- scaleDist <-  numeric(120)
  multiplicity <- omega_k <- Cstar <- sdlog_e <- MAE_p <- numeric(120)
  
  
  for (i in 1:60){
    nTotSites <- 40
    nS_e <- 40
    nS_k <- 0
    ind_sim <- 2
    
    for (ind_m in c(1,3)){
      fnam <- paste0("results/multiplicity_BT/out",i,"_s",ind_sim,"_e",nS_e,"_k",nS_k,"_m",ind_m,".rda")
      if (file.exists(fnam)){
        load(fnam)
        if (!is.null(out)){
          if (length(out$outMCMC)>0){
            out$outMCMC <- NULL
            save(out,file=fnam)}
          k <- k + 1
          ID[k] <- i
          sim[k] <- out$ind_str
          ind_OCN[k] <- out$ind_OCN
          nSites_e[k] <- out$nSites_e
          nSites_k[k] <- out$nSites_k
          tau_sol[k] <- out$tau_sol
          tau[k] <- out$param[1]
          log_p0[k] <- out$param[2]
          k_mean[k] <- out$k_mean
          k_cov[k] <- out$k_cov
          tau_sol[k] <- out$tau_sol
          MAE_p[k] <- mean(abs(out$p_map-out$p_sol))
          nNodes[k] <- out$nNodes
          scaleDist[k] <- out$scaleDist
          multiplicity[k] <- ind_m
          Cstar[k] <- out$Cstar
          omega_k[k] <- out$omega_k
          sdlog_e[k] <- out$sdlog_e
          
        }
      }
      fnam <- paste0("../results/multiplicity_optim/out",i,"_s",ind_sim,"_e",nS_e,"_k",nS_k,"_m",ind_m,".rda")
      if (file.exists(fnam)){
        load(fnam)
        if (!is.null(out)){
          MAE_p_optim[k] <- mean(abs(out$p-out$p_sol))
        }
      }
    }
    
  }
  save(MAE_p,MAE_p_optim,multiplicity, file="results/results_multiplicity.rda")
} else {load("results/results_multiplicity.rda")}

# Fig. S1 ####
pdf(file="Fig_Sx_draft.pdf",width=12/2.54, height=9/2.54)
boxplot(MAE_p[multiplicity==1],MAE_p[multiplicity==3], 
        MAE_p_optim[multiplicity==1],MAE_p_optim[multiplicity==3],outline=F,
        names= c("m1_BT","m3_BT","m1_optim","m3_optim"),
        col=rep(hcl.colors(2,"Blue-Red 2"),each=2),
        ylab="MAE")
dev.off()