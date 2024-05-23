rm(list=ls())

library(OCNet)
library(car)

source("eval_wilcox_paired.R")

lastSim <- 60

if (!file.exists("results/results_df.rda")){
sV <- 1:lastSim
simVec <- c(sV,100+sV,200+sV,300+sV,400+sV,500+sV,600+sV,700+sV)
k <- 0
ID <- nNodes <- scaleDist <- ind_OCN <- numeric(48000)
nSites_e <- nSites_k <- sim <- numeric(48000)
tau_sol <- k_mean <- k_cov <- scaleDist <-  numeric(48000)
multiplicity <- omega_k <- Cstar <- sdlog_e <- numeric(48000)
MAE_p <- numeric(48000)

for (i in simVec){
  for (nTotSites in c(20,40)){
    for (ind_e in 1:5){
      nS_e <- (ind_e-1)/4*nTotSites
      nS_k <- nTotSites - nS_e
      for (ind_sim in 1:5){
        for (ind_m in c(1,3)){
          fnam <- paste0("results/out",i,"_s",ind_sim,"_e",nS_e,"_k",nS_k,"_m",ind_m,".rda")
          if (file.exists(fnam)){
            load(fnam)
            if (!is.null(out)){
              k <- k + 1
              ID[k] <- i
              sim[k] <- out$ind_str
              ind_OCN[k] <- out$ind_OCN
              nSites_e[k] <- out$nSites_e
              nSites_k[k] <- out$nSites_k
              tau_sol[k] <- out$tau_sol

              k_mean[k] <- out$k_mean
              k_cov[k] <- out$k_cov
              tau_sol[k] <- out$tau_sol

              MAE_p[k] <- mean(abs(out$p-out$p_sol))

              nNodes[k] <- out$nNodes
              scaleDist[k] <- out$scaleDist
              multiplicity[k] <- ind_m
              Cstar[k] <- out$Cstar
              omega_k[k] <- out$omega_k
              sdlog_e[k] <- out$sdlog_e
            }
          }
        }
      }
    }
  }
}

ID <- ID[1:k]; nSites_e <- nSites_e[1:k]; 
nSites_k <- nSites_k[1:k]; sim <- sim[1:k]; nSites_tot <- nSites_e + nSites_k
k_mean <- k_mean[1:k]; k_cov <- k_cov[1:k]; tau_sol <- tau_sol[1:k]
MAE_p <- MAE_p[1:k];
scaleDist <- scaleDist[1:k]; nNodes <- nNodes[1:k]; ind_OCN <- ind_OCN[1:k]
Cstar <- Cstar[1:k]; omega_k <- omega_k[1:k]; sdlog_e <- sdlog_e[1:k]
multiplicity <- multiplicity[1:k]; 

frac_e <- nSites_e/nSites_tot

df <- data.frame(ID=ID, sim=sim, multiplicity=multiplicity,
                 nSites_e=nSites_e, nSites_k=nSites_k, nSites_tot=nSites_tot,
                 frac_e=frac_e, k_mean=k_mean, k_cov=k_cov,
                 scaleDist=scaleDist, omega_k=omega_k, Cstar=Cstar, sdlog_e=sdlog_e,
                 tau_sol=tau_sol, MAE_p=MAE_p)

} else {load("results/results_df.rda")}

for(i in 1:length(df)) {assign(names(df)[i], df[[i]])} # unwrap results

kol <- hcl.colors(5,"Temps")
kol1 <- colorRampPalette(c("black",kol[1],"white")); kol1 <- kol1(8); kol1 <- kol1[-c(1,2,8)]
kol2 <- colorRampPalette(c("black",kol[2],"white")); kol2 <- kol2(8); kol2 <- kol2[-c(1,2,8)]
kol3 <- colorRampPalette(c("black",kol[3],"white")); kol3 <- kol3(8); kol3 <- kol3[-c(1,2,8)]
kol4 <- colorRampPalette(c("black",kol[4],"white")); kol4 <- kol4(8); kol4 <- kol4[-c(1,2,8)]
kol5 <- colorRampPalette(c("black",kol[5],"white")); kol5 <- kol5(8); kol5 <- kol5[-c(1,2,8)]

# Best strategy ####
pdf(file="Fig3_draft.pdf",width=20/2.54,height=16/2.54)
par(mfrow=c(2,2))
boxplot(MAE_p ~ sim, outline=F, log="y",col=kol,
        names=c("eR-dR","eT-dT","eS-dR","eR-dT","eB-dB"),# R: random, S: strategic, B: balanced
        ylim=c(0.01,1000))
best_median(df, 1:nrow(df), 0, "sim", 1:5)
wt <- eval_wilcox_paired(df, 1:nrow(df), "sim", 1:5)
draw_significance(wt,0,coeff=1,height=1000,offset=50)

boxplot(MAE_p ~  frac_e + sim,  outline=F, log="y",
        col=c(kol1,kol2,kol3,kol4,kol5),
        names=rep(c("e0","e25","e50","e75","e100"),5),ylim=c(0.01,1000),
        main = "eR-dR    /    eT-dT     /    eT-dR    /    eR-dT    /    eB-dB")
abline(v=c(5.5,10.5,15.5,20.5))
best_median(df, sim==1, 0, "frac_e", (0:4)/4); best_median(df, sim==2, 5, "frac_e", (0:4)/4)
best_median(df, sim==3, 10, "frac_e", (0:4)/4); best_median(df, sim==4, 15, "frac_e", (0:4)/4)
best_median(df, sim==5, 20, "frac_e", (0:4)/4)
wt1 <- eval_wilcox_paired(df, sim==1, "frac_e", (0:4)/4)
draw_significance(wt1,1,coeff=4,height=1000,offset=100)
wt2 <- eval_wilcox_paired(df, sim==2, "frac_e", (0:4)/4)
draw_significance(wt2,6,coeff=4,height=1000,offset=100)
wt3 <- eval_wilcox_paired(df, sim==3, "frac_e", (0:4)/4)
draw_significance(wt3,11,coeff=4,height=1000,offset=100)
wt4 <- eval_wilcox_paired(df, sim==4, "frac_e", (0:4)/4)
draw_significance(wt4,16,coeff=4,height=1000,offset=100)
wt5 <- eval_wilcox_paired(df, sim==5, "frac_e", (0:4)/4)
draw_significance(wt5,21,coeff=4,height=1000,offset=100)

boxplot(MAE_p ~  nSites_tot + sim,log="y",outline=F,
        col=rep(kol,each=2),#c(kol1[c(2,4)],kol2[c(2,4)],kol3[c(2,4)],kol4[c(2,4)],kol5[c(2,4)]),
        ylim=c(0.01,100),names=rep(c("20","40"),5),
        main = c("n Tot Sites = 20     /    n Tot sites = 40"))
abline(v=c(2.5,4.5,6.5,8.5))
best_median(df, sim==1, 0,  "nSites_tot", c(20,40)); best_median(df, sim==2, 2,  "nSites_tot", c(20,40))
best_median(df, sim==3, 4,  "nSites_tot", c(20,40)); best_median(df, sim==4, 6,  "nSites_tot", c(20,40))
best_median(df, sim==5, 8,  "nSites_tot", c(20,40))
wt1 <- eval_wilcox_paired(df, sim==1, "nSites_tot", c(20,40))
wt2 <- eval_wilcox_paired(df, sim==2, "nSites_tot", c(20,40))
wt3 <- eval_wilcox_paired(df, sim==3, "nSites_tot", c(20,40))
wt4 <- eval_wilcox_paired(df, sim==4, "nSites_tot", c(20,40))
wt5 <- eval_wilcox_paired(df, sim==5, "nSites_tot", c(20,40))
draw_significance(wt1,0,coeff=0.05,height=100,offset=20)
draw_significance(wt2,2,coeff=0.05,height=100,offset=20)
draw_significance(wt3,4,coeff=0.05,height=100,offset=20)
draw_significance(wt4,6,coeff=0.05,height=100,offset=20)
draw_significance(wt5,8,coeff=0.05,height=100,offset=20)

boxplot(MAE_p ~  multiplicity + sim,log="y",outline=F,
        col=rep(kol,each=2),#c(kol1[c(2,4)],kol2[c(2,4)],kol3[c(2,4)],kol4[c(2,4)],kol5[c(2,4)]),
        ylim=c(0.01,100),names=rep(c("1","3"),5),
        main = c("multiplicity = 1     /    multiplicity = 3"))
abline(v=c(2.5,4.5,6.5,8.5))
best_median(df, sim==1, 0,"multiplicity", c(1,3)); best_median(df, sim==2, 2,"multiplicity", c(1,3))
best_median(df, sim==3, 4,"multiplicity", c(1,3)); best_median(df, sim==4, 6,"multiplicity", c(1,3))
best_median(df, sim==5, 8,"multiplicity", c(1,3))
wt1 <- eval_wilcox_paired(df, sim==1, "multiplicity", c(1,3))
wt2 <- eval_wilcox_paired(df, sim==2, "multiplicity", c(1,3))
wt3 <- eval_wilcox_paired(df, sim==3, "multiplicity", c(1,3))
wt4 <- eval_wilcox_paired(df, sim==4, "multiplicity", c(1,3))
wt5 <- eval_wilcox_paired(df, sim==5, "multiplicity", c(1,3))
draw_significance(wt1,0.5,coeff=0.5,height=100,offset=20)
draw_significance(wt2,2.5,coeff=0.5,height=100,offset=20)
draw_significance(wt3,4.5,coeff=0.5,height=100,offset=20)
draw_significance(wt4,6.5,coeff=0.5,height=100,offset=20)
draw_significance(wt5,8.5,coeff=0.5,height=100,offset=20)
dev.off()

# Prediction power for different contexts ####
pdf(file="Fig4_draft.pdf",width=20/2.54,height=16/2.54)
par(mfrow=c(3,3))
boxplot(MAE_p ~ (k_mean>median(k_mean)) + sim,
        log="y",outline=F,
        col=rep(kol,each=2),ylim=c(0.01,100),
        names=rep(c("low","high"),5),
        main = c("low k_mean     /    high k_mean"))
abline(v=c(2.5,4.5,6.5,8.5))
best_median(df,sim==1,0,"k_mean",split.median = TRUE); best_median(df,sim==2,2,"k_mean",split.median = TRUE)
best_median(df,sim==3,4,"k_mean",split.median = TRUE); best_median(df,sim==4,6,"k_mean",split.median = TRUE)
best_median(df,sim==5,8,"k_mean",split.median = TRUE)

boxplot(MAE_p ~ (k_cov>median(k_cov)) + sim,
        log="y",outline=F,
         col=rep(kol,each=2),ylim=c(0.01,100),
        names=rep(c("low","high"),5),
        main = c("low k_cov     /    high k_cov"))
abline(v=c(2.5,4.5,6.5,8.5))
best_median(df,sim==1,0,"k_cov",split.median = TRUE); best_median(df,sim==2,2,"k_cov",split.median = TRUE)
best_median(df,sim==3,4,"k_cov",split.median = TRUE); best_median(df,sim==4,6,"k_cov",split.median = TRUE)
best_median(df,sim==5,8,"k_cov",split.median = TRUE)

boxplot(MAE_p ~ (scaleDist>median(scaleDist)) + sim,
        log="y",outline=F,
         col=rep(kol,each=2),ylim=c(0.01,100),
        names=rep(c("low","high"),5),
        main = c("low scaleDist     /    high scaleDist"))
abline(v=c(2.5,4.5,6.5,8.5))
best_median(df,sim==1,0,"scaleDist",split.median = TRUE); best_median(df,sim==2,2,"scaleDist",split.median = TRUE)
best_median(df,sim==3,4,"scaleDist",split.median = TRUE); best_median(df,sim==4,6,"scaleDist",split.median = TRUE)
best_median(df,sim==5,8,"scaleDist",split.median = TRUE)

boxplot(MAE_p ~  (omega_k>median(omega_k)) + sim,
        log="y",outline=F,
         col=rep(kol,each=2),ylim=c(0.01,100),
        names=rep(c("low","high"),5),
        main = c("low omega_k     /    high omega_k"))
abline(v=c(2.5,4.5,6.5,8.5))
best_median(df,sim==1,0,"omega_k",split.median = TRUE); best_median(df,sim==2,2,"omega_k",split.median = TRUE)
best_median(df,sim==3,4,"omega_k",split.median = TRUE); best_median(df,sim==4,6,"omega_k",split.median = TRUE)
best_median(df,sim==5,8,"omega_k",split.median = TRUE)

boxplot(MAE_p ~ (sdlog_e>median(sdlog_e)) + sim,
        log="y",outline=F,
         col=rep(kol,each=2),ylim=c(0.01,100),
        names=rep(c("low","high"),5),
        main = c("low sdlog_e     /    high sdlog_e"))
abline(v=c(2.5,4.5,6.5,8.5))
best_median(df,sim==1,0,"sdlog_e",split.median = TRUE); best_median(df,sim==2,2,"sdlog_e",split.median = TRUE)
best_median(df,sim==3,4,"sdlog_e",split.median = TRUE); best_median(df,sim==4,6,"sdlog_e",split.median = TRUE)
best_median(df,sim==5,8,"sdlog_e",split.median = TRUE)

boxplot(MAE_p ~ (Cstar>median(Cstar)) + sim,
        log="y",outline=F,
        col=rep(kol,each=2),ylim=c(0.01,100),
        names=rep(c("low","high"),5),
        main = c("low Cstar     /    high Cstar"))
abline(v=c(2.5,4.5,6.5,8.5))
best_median(df,sim==1,0,"Cstar",split.median = TRUE); best_median(df,sim==2,2,"Cstar",split.median = TRUE)
best_median(df,sim==3,4,"Cstar",split.median = TRUE); best_median(df,sim==4,6,"Cstar",split.median = TRUE)
best_median(df,sim==5,8,"Cstar",split.median = TRUE)

boxplot(MAE_p ~ (tau_sol>median(tau_sol)) + sim,
        log="y",outline=F,
        col=rep(kol,each=2),ylim=c(0.01,100),
        names=rep(c("low","high"),5),
        main = c("low tau_sol     /    high tau_sol"))
abline(v=c(2.5,4.5,6.5,8.5))
best_median(df,sim==1,0,"tau_sol",split.median = TRUE); best_median(df,sim==2,2,"tau_sol",split.median = TRUE)
best_median(df,sim==3,4,"tau_sol",split.median = TRUE); best_median(df,sim==4,6,"tau_sol",split.median = TRUE)
best_median(df,sim==5,8,"tau_sol",split.median = TRUE)
dev.off()

# table with ratios of medians
tbl1 <- tbl2 <- matrix(0,5,7)
pars <- c("k_mean","k_cov","scaleDist","omega_k","sdlog_e","Cstar","tau_sol")
for (i in 1:5){
  for (j in 1:7){
    gg <- pars[j]
    tbl1[i,j] <- median(MAE_p[sim==i & df[[gg]]<=median(df[[gg]])]/MAE_p[sim==i & df[[gg]]>median(df[[gg]])])
    tbl2[i,j] <- median(MAE_p[sim==i & df[[gg]]<=median(df[[gg]])])/median(MAE_p[sim==i & df[[gg]]>median(df[[gg]])])
  }
}
round(1000*tbl2)/1000

# Focus on Strategy 2 (remove frac_e = 0) ####
pdf(file="Fig5_draft.pdf",width=20/2.54,height=16/2.54)
par(mfrow=c(3,3))
wt1 <- eval_wilcox_paired(df, (k_mean<median(k_mean)) & sim==2)
wt2 <- eval_wilcox_paired(df, (k_mean>median(k_mean)) & sim==2)
boxplot(MAE_p ~ df$frac_e  +(k_mean>median(k_mean)), outline=F,
        subset=(sim==2 & frac_e!=0), col=rep(kol2[2:5],2),ylim=c(0,10),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low k_mean     /    high k_mean"))
abline(v=4.5)
best_median(df, (k_mean<median(k_mean)) & sim==2, 0); best_median(df, (k_mean>median(k_mean)) & sim==2, 4)
draw_significance(wt1,0); draw_significance(wt2,4)

wt1 <- eval_wilcox_paired(df, (k_cov<median(k_cov)) & sim==2)
wt2 <- eval_wilcox_paired(df, (k_cov>median(k_cov)) & sim==2)
boxplot(MAE_p ~ df$frac_e + (k_cov>median(k_cov)), outline=F,
        subset=(sim==2 & frac_e!=0), col=rep(kol2[2:5],2),ylim=c(0,10),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low k_cov     /    high k_cov"))
abline(v=4.5)
best_median(df, (k_cov<median(k_cov)) & sim==2, 0); best_median(df, (k_cov>median(k_cov)) & sim==2, 4)
draw_significance(wt1,0); draw_significance(wt2,4)

wt1 <- eval_wilcox_paired(df, (scaleDist<median(scaleDist)) & sim==2)
wt2 <- eval_wilcox_paired(df, (scaleDist>median(scaleDist)) & sim==2)
boxplot(MAE_p ~ df$frac_e + (scaleDist>median(scaleDist)),outline=F,
        subset=(sim==2 & frac_e!=0), col=rep(kol2[2:5],2),ylim=c(0,10),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low scaleDist     /    high scaleDist"))
abline(v=4.5)
best_median(df, (scaleDist<median(scaleDist)) & sim==2, 0); best_median(df, (scaleDist>median(scaleDist)) & sim==2, 4)
draw_significance(wt1,0); draw_significance(wt2,4)

wt1 <- eval_wilcox_paired(df, (omega_k<median(omega_k)) & sim==2)
wt2 <- eval_wilcox_paired(df, (omega_k>median(omega_k)) & sim==2)
boxplot(MAE_p ~ df$frac_e + (omega_k>median(omega_k)),outline=F,
        subset=(sim==2 & frac_e!=0), col=rep(kol2[2:5],2),ylim=c(0,10),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low omega_k     /    high omega_k"))
abline(v=4.5)
best_median(df, (omega_k<median(omega_k)) & sim==2, 0); best_median(df, (omega_k>median(omega_k)) & sim==2, 4)
draw_significance(wt1,0); draw_significance(wt2,4)

wt1 <- eval_wilcox_paired(df, (sdlog_e<median(sdlog_e)) & sim==2)
wt2 <- eval_wilcox_paired(df, (sdlog_e>median(sdlog_e)) & sim==2)
boxplot(MAE_p ~ frac_e + (sdlog_e>median(sdlog_e)),outline=F,
        subset=(sim==2 & frac_e!=0), col=rep(kol2[2:5],2),ylim=c(0,10),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low sdlog_e     /    high sdlog_e"))
abline(v=4.5)
best_median(df, (sdlog_e<median(sdlog_e)) & sim==2, 0); best_median(df, (sdlog_e>median(sdlog_e)) & sim==2, 4)
draw_significance(wt1,0); draw_significance(wt2,4)

wt1 <- eval_wilcox_paired(df, (Cstar<median(Cstar)) & sim==2)
wt2 <- eval_wilcox_paired(df, (Cstar>median(Cstar)) & sim==2)
boxplot(MAE_p ~ frac_e + (Cstar>median(Cstar)),
        outline=F,ylim=c(0,10),
        subset=(sim==2 & frac_e!=0), col=rep(kol2[2:5],2),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low Cstar     /    high Cstar"))
abline(v=4.5)
best_median(df, (Cstar<median(Cstar)) & sim==2, 0); best_median(df, (Cstar>median(Cstar)) & sim==2, 4)
draw_significance(wt1,0); draw_significance(wt2,4)

wt1 <- eval_wilcox_paired(df, (tau_sol<median(tau_sol)) & sim==2)
wt2 <- eval_wilcox_paired(df, (tau_sol>median(tau_sol)) & sim==2)
boxplot(MAE_p ~ frac_e + (tau_sol>median(tau_sol)),
        outline=F,ylim=c(0,10),
        subset=(sim==2 & frac_e!=0), col=rep(kol2[2:5],2),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low tau_sol     /    high tau_sol"))
abline(v=4.5)
best_median(df, (tau_sol<median(tau_sol)) & sim==2, 0); best_median(df, (tau_sol>median(tau_sol)) & sim==2, 4)
draw_significance(wt1,0); draw_significance(wt2,4)

wt1 <- eval_wilcox_paired(df, nSites_tot==20 & sim==2)
wt2 <- eval_wilcox_paired(df, nSites_tot==40 & sim==2)
boxplot(MAE_p ~  frac_e + nSites_tot,outline=F,
        subset=(sim==2 & frac_e!=0), col=rep(kol2[2:5],2),ylim=c(0,10),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("n Tot Sites = 20     /    n Tot sites = 40"))
abline(v=4.5)
best_median(df, nSites_tot==20 & sim==2, 0); best_median(df, nSites_tot==40 & sim==2, 4)
draw_significance(wt1,0); draw_significance(wt2,4)

wt1 <- eval_wilcox_paired(df, multiplicity==1 & sim==2)
wt2 <- eval_wilcox_paired(df, multiplicity==3 & sim==2)
boxplot(MAE_p ~ frac_e + multiplicity,outline=F,
        subset=(sim==2 & frac_e!=0), col=rep(kol2[2:5],2),ylim=c(0,10),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("multiplicity = 1     /    multiplicity = 3"))
abline(v=4.5)
best_median(df, multiplicity==1 & sim==2, 0); best_median(df, multiplicity==3 & sim==2, 4)
draw_significance(wt1,0); draw_significance(wt2,4)
dev.off()

## Focus on Strategy 1 ####
pdf(file="FigS_s1_draft.pdf",width=20/2.54,height=16/2.54)
par(mfrow=c(3,3))
wt1 <- eval_wilcox_paired(df, (k_mean<median(k_mean)) & sim==1)
wt2 <- eval_wilcox_paired(df, (k_mean>median(k_mean)) & sim==1)
boxplot(MAE_p ~ df$frac_e  +(k_mean>median(k_mean)), outline=F,
        subset=(sim==1 & frac_e!=0), col=rep(kol1[2:5],2),ylim=c(0,30),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low k_mean     /    high k_mean"))
abline(v=4.5)
best_median(df, (k_mean<median(k_mean)) & sim==1, 0); best_median(df, (k_mean>median(k_mean)) & sim==1, 4)
draw_significance(wt1,0, height=30); draw_significance(wt2,4, height=30)

wt1 <- eval_wilcox_paired(df, (k_cov<median(k_cov)) & sim==1)
wt2 <- eval_wilcox_paired(df, (k_cov>median(k_cov)) & sim==1)
boxplot(MAE_p ~ df$frac_e + (k_cov>median(k_cov)), outline=F,
        subset=(sim==1 & frac_e!=0), col=rep(kol1[2:5],2),ylim=c(0,30),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low k_cov     /    high k_cov"))
abline(v=4.5)
best_median(df, (k_cov<median(k_cov)) & sim==1, 0); best_median(df, (k_cov>median(k_cov)) & sim==1, 4)
draw_significance(wt1,0, height=30); draw_significance(wt2,4, height=30)

wt1 <- eval_wilcox_paired(df, (scaleDist<median(scaleDist)) & sim==1)
wt2 <- eval_wilcox_paired(df, (scaleDist>median(scaleDist)) & sim==1)
boxplot(MAE_p ~ df$frac_e + (scaleDist>median(scaleDist)),outline=F,
        subset=(sim==1 & frac_e!=0), col=rep(kol1[2:5],2),ylim=c(0,30),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low scaleDist     /    high scaleDist"))
abline(v=4.5)
best_median(df, (scaleDist<median(scaleDist)) & sim==1, 0); best_median(df, (scaleDist>median(scaleDist)) & sim==1, 4)
draw_significance(wt1,0, height=30); draw_significance(wt2,4, height=30)

wt1 <- eval_wilcox_paired(df, (omega_k<median(omega_k)) & sim==1)
wt2 <- eval_wilcox_paired(df, (omega_k>median(omega_k)) & sim==1)
boxplot(MAE_p ~ df$frac_e + (omega_k>median(omega_k)),outline=F,
        subset=(sim==1 & frac_e!=0), col=rep(kol1[2:5],2),ylim=c(0,30),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low omega_k     /    high omega_k"))
abline(v=4.5)
best_median(df, (omega_k<median(omega_k)) & sim==1, 0); best_median(df, (omega_k>median(omega_k)) & sim==1, 4)
draw_significance(wt1,0, height=23); draw_significance(wt2,4, height=30)

wt1 <- eval_wilcox_paired(df, (sdlog_e<median(sdlog_e)) & sim==1)
wt2 <- eval_wilcox_paired(df, (sdlog_e>median(sdlog_e)) & sim==1)
boxplot(MAE_p ~ frac_e + (sdlog_e>median(sdlog_e)),outline=F,
        subset=(sim==1 & frac_e!=0), col=rep(kol1[2:5],2),ylim=c(0,30),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low sdlog_e     /    high sdlog_e"))
abline(v=4.5)
best_median(df, (sdlog_e<median(sdlog_e)) & sim==1, 0); best_median(df, (sdlog_e>median(sdlog_e)) & sim==1, 4)
draw_significance(wt1,0, height=30); draw_significance(wt2,4, height=30)

wt1 <- eval_wilcox_paired(df, (Cstar<median(Cstar)) & sim==1)
wt2 <- eval_wilcox_paired(df, (Cstar>median(Cstar)) & sim==1)
boxplot(MAE_p ~ frac_e + (Cstar>median(Cstar)),
        outline=F,ylim=c(0,30),
        subset=(sim==1 & frac_e!=0), col=rep(kol1[2:5],2),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low Cstar     /    high Cstar"))
abline(v=4.5)
best_median(df, (Cstar<median(Cstar)) & sim==1, 0); best_median(df, (Cstar>median(Cstar)) & sim==1, 4)
draw_significance(wt1,0, height=30); draw_significance(wt2,4, height=30)

wt1 <- eval_wilcox_paired(df, (tau_sol<median(tau_sol)) & sim==1)
wt2 <- eval_wilcox_paired(df, (tau_sol>median(tau_sol)) & sim==1)
boxplot(MAE_p ~ frac_e + (tau_sol>median(tau_sol)),
        outline=F,ylim=c(0,30),
        subset=(sim==1 & frac_e!=0), col=rep(kol1[2:5],2),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("low tau_sol     /    high tau_sol"))
abline(v=4.5)
best_median(df, (tau_sol<median(tau_sol)) & sim==1, 0); best_median(df, (tau_sol>median(tau_sol)) & sim==1, 4)
draw_significance(wt1,0, height=30); draw_significance(wt2,4, height=30)

wt1 <- eval_wilcox_paired(df, nSites_tot==20 & sim==1)
wt2 <- eval_wilcox_paired(df, nSites_tot==40 & sim==1)
boxplot(MAE_p ~  frac_e + nSites_tot,outline=F,
        subset=(sim==1 & frac_e!=0), col=rep(kol1[2:5],2),ylim=c(0,30),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("n Tot Sites = 20     /    n Tot sites = 40"))
abline(v=4.5)
best_median(df, nSites_tot==20 & sim==1, 0); best_median(df, nSites_tot==40 & sim==1, 4)
draw_significance(wt1,0, height=30); draw_significance(wt2,4, height=30)

wt1 <- eval_wilcox_paired(df, multiplicity==1 & sim==1)
wt2 <- eval_wilcox_paired(df, multiplicity==3 & sim==1)
boxplot(MAE_p ~ frac_e + multiplicity,outline=F,
        subset=(sim==1 & frac_e!=0), col=rep(kol1[2:5],2),ylim=c(0,30),
        names=rep(c("e25","e50","e75","e100"),2),
        main = c("multiplicity = 1     /    multiplicity = 3"))
abline(v=4.5)
best_median(df, multiplicity==1 & sim==1, 0); best_median(df, multiplicity==3 & sim==1, 4)
draw_significance(wt1,0, height=30); draw_significance(wt2,4, height=30)
dev.off()

