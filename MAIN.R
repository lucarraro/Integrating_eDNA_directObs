rm(list=ls())

library(OCNet)
library(rivnet)
library(LaplacesDemon)
library(eDITH)

if (packageVersion("eDITH") <= "0.3.0"){
  stop("please install eDITH v1.0.0 or higher")
}

job <- 4 # set from 1 to 8 to run all simulations
nuFluvial <- 0.5 

vec_sim <- ((job-1)*100+1):(job*100)

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

  # load OCN and create river object
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
  
  sitesUpstream <- which(OCN$AG$A <= median(OCN$AG$A))
  sitesDownstream <- which(OCN$AG$A > median(OCN$AG$A))

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

  conc <- run_eDITH_single(param, OCN, covariates, Z.normalize=FALSE) 
  C_sol <- conc$C

  # we generate C_obs 3 times, and only pick the first set if multiplicity=1
  C_obs <- rep(C_sol, 3)
  C_obs[runif(length(C_obs)) < exp(-C_obs/Cstar)] <- 0 # non detection probability
  C_obs <- exp(rnorm(length(C_obs),log(C_obs), sdlog_e))

  dens_obs <- rnbinom(nNodes, size=p_sol/(omega_k-1), prob=1/omega_k)

  for (nTotSites in c(20,40)){
    for (ind_e in 1:5){
      nSites_e <- (ind_e-1)/4*nTotSites
      nSites_k <- nTotSites - nSites_e
      for (ind_str in 1:5){
        for (multiplicity in c(1,3)){
          fnam <- paste0("results/out",ind_sim,"_s",ind_str,
                         "_e",nSites_e,"_k",nSites_k,"_m",multiplicity,".rda")
          if (!(file.exists(fnam))){
            out <- NULL
            save(out, file=fnam)
            message(sprintf('i: %d  -  sim: %d  -  eSites: %d  -  kSites: %d  -  multiplicity: %d \n',
                            ind_sim,ind_str,nSites_e,nSites_k,multiplicity),appendLF=F)
            if (ind_str==1){ #FULLY RANDOM
              sampling.sites_e <- sample(nNodes,nSites_e)
              sampling.sites_k <- sample(nNodes,nSites_k)
            } else if (ind_str==2){ # EXPECTED BEST
              sampling.sites_e <- sampling_strategy_eDNA(OCN, nSites_e)
              sampling.sites_k <- sampling_strategy_direct(OCN, nSites_k)
            } else if (ind_str==3){ # BEST EDNA, RANDOM KICKNET
              sampling.sites_e <- sampling_strategy_eDNA(OCN, nSites_e)
              sampling.sites_k <- sample(nNodes,nSites_k)
            } else if (ind_str==4){ # RANDOM EDNA, BEST KICKNET
              sampling.sites_e <- sample(nNodes,nSites_e)
              sampling.sites_k <- sampling_strategy_direct(OCN, nSites_k)
            } else if (ind_str==5){ # BALANCED EDNA & KICKNET
              sseu <- sample(sitesUpstream, floor(nSites_e/2))
              ssed <- sample(sitesDownstream, nSites_e-length(sseu)) # one more site downstream if uneven
              sampling.sites_e <- c(sseu, ssed)
              ssku <- sample(sitesUpstream, floor(nSites_k/2))
              sskd <- sample(sitesDownstream, nSites_k-length(ssku)) # one more site downstream if uneven
              sampling.sites_k <- c(ssku, sskd)
            }
            if (multiplicity==3){
              ind.sites_e <- c(sampling.sites_e, nNodes + sampling.sites_e, 2*nNodes + sampling.sites_e)
              ss_e <- rep(sampling.sites_e, 3)
            } else {
              ind.sites_e <- ss_e <- sampling.sites_e}

            data <- data.frame(ID=c(ss_e, sampling.sites_k),
                               values=c(C_obs[ind.sites_e], dens_obs[sampling.sites_k]),
                               type=c(rep("e",length(ss_e)), rep("d",nSites_k)))

            set.seed(ind_sim) # for each map, we fit in the same way (same strategy -> same result)
            out <- run_eDITH_optim_joint(data, OCN, ll.type="lnorm", no.det=T,
                                         par.AEM = list(weight="exponential"),
                                         tau.prior = list(spec="unif",min=1, max=24), # same boundaries as tau_sol
                                         log_p0.prior = list(spec="unif",min=-5, max=5),
                                         n.AEM = 20,
                                         alpha.prior = list(spec="unif",min=1-1e-6,max=1+1e-6),
                                         verbose=T, n.attempts = 20) # attempts are 4x5

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

            #out$sampling.sites <- sampling.sites # it's already within data!
            out$covariates <- NULL # save space
            out$source.area <- NULL
            out$probDet <- NULL # can be recalculated if needed
            save(out,file=fnam)
          } else {
            message("   Already done! \n", appendLF=F)
          }
          message("\n",appendLF=F)
        }
      }
    }
  }
}





