########################### BayesKin ###########################
## Bayesian estimates of rigid body kinematics
## Andy Pohl
## UofC - Faculty of Kinesiology

## Written by: Andy Pohl
## UofC - Faculty of Kinesiology
## June-Dec 2020
## Revision 1: April 2020
################################################################

##################### process_results.r  ###################
## processes_results table for a batch of simulations
#############################################################

# preliminaries 
rm(list = ls()) # clear workspace

library('rjags') # load rstan package to analysie rstan objects
library('coda') # load coda package for mcmc analysis

RESULTS_DIR = "" # location of results directory
POSE_DIRS = c("SingleLink_Mdl", "DoubleLink_Mdl", "TripleLink_Mdl")  # Folder names for the various model configurations
MODELS = c("LS", "Bayes_p1", "Bayes_p2", "Bayes_p3", "Bayes_p4", "Bayes_p5") # Names of Models used
file_ends = c('LSModel.rda', 'BayesModel_p1.rda', 'BayesModel_p2.rda', 'BayesModel_p3.rda', 'BayesModel_p4.rda', 'BayesModel_p5.rda')

PARMS = c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma') # Parameters of interest

nfiles = 1000 # the number of results files to process
nmodels = length(MODELS)

################################################################################
# Single Link
################################################################################
results_single = data.frame(pose = rep('SingleLink', nmodels*4*nfiles),
                            seed = rep(1:nfiles, each = nmodels*4),
                            model = rep(rep(MODELS, each = 4),nfiles),
                            parameter = rep(c('r1', 'r2', 'theta1', 'sigma'), nfiles*nmodels),
                            est = NA,
                            quantile = matrix(NA, nfiles*4*nmodels, 5),
                            neff = NA,
                            rhat = NA,
                            ts_se = NA, 
                            run_time = NA)

files = list.files(path = paste0(RESULTS_DIR, "SingleLink_Mdl/Bayes_p3"))
seeds = as.numeric(gsub("(.+?)(\\_.*)", "\\1", files))
seeds = seeds[order(seeds)]
rowcnt = 1
for(i in 1:nfiles){
    nlinks = 1
    seed = seeds[i]
    print(sprintf("Processing Single Link file %.0f of %.0f", i, nfiles))
    
    # process LS results
    lsmdl = readRDS(paste0(RESULTS_DIR, "SingleLink_Mdl/LS/", seed, "_LSModel.rda"))
    est = c(lsmdl$r.hat, lsmdl$theta.hat, lsmdl$sigma.hat)
    int = lsmdl$intervals[1:(nlinks+2),]
    results_single[rowcnt:(rowcnt+3), c('est')] = est
    results_single[rowcnt:(rowcnt+2), c('quantile.1', 'quantile.5')] = int
    results_single[rowcnt:(rowcnt+3), c('run_time')] = lsmdl$time
    rowcnt = rowcnt + 4
    
    # process Bayes P1
    p1mdl =  readRDS(paste0(RESULTS_DIR, "SingleLink_Mdl/Bayes_p1/", seed, "_BayesModel_p1.rda"))
    p1sum = summary(p1mdl$fit)
    results_single[rowcnt:(rowcnt+3),c('est', 'ts_se')] = p1sum$statistics[c('r[1]', 'r[2]', 'theta', 'sigma'),c('Mean', 'Time-series SE')]
    results_single[rowcnt:(rowcnt+3),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p1sum$quantiles[c('r[1]', 'r[2]', 'theta', 'sigma'),]
    results_single[rowcnt:(rowcnt+3), c('neff')] = effectiveSize(p1mdl$fit)[c('r[1]', 'r[2]', 'theta', 'sigma')]
    results_single[rowcnt:(rowcnt+3), c('rhat')] = gelman.diag(p1mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta', 'sigma'),1]
    results_single[rowcnt:(rowcnt+3), c('run_time')] = p1mdl$time 
    rowcnt = rowcnt + 4
    
    # process Bayes P2
    p2mdl =  readRDS(paste0(RESULTS_DIR, "SingleLink_Mdl/Bayes_p2/", seed, "_BayesModel_p2.rda"))
    p2sum = summary(p2mdl$fit)
    results_single[rowcnt:(rowcnt+3),c('est', 'ts_se')] = p2sum$statistics[c('r[1]', 'r[2]', 'theta', 'sigma'),c('Mean', 'Time-series SE')]
    results_single[rowcnt:(rowcnt+3),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p2sum$quantiles[c('r[1]', 'r[2]', 'theta', 'sigma'),]
    results_single[rowcnt:(rowcnt+3), c('neff')] = effectiveSize(p2mdl$fit)[c('r[1]', 'r[2]', 'theta', 'sigma')]
    results_single[rowcnt:(rowcnt+3), c('rhat')] = gelman.diag(p2mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta', 'sigma'),1]
    results_single[rowcnt:(rowcnt+3), c('run_time')] = p2mdl$time 
    rowcnt = rowcnt + 4
    
    # process Bayes P3
    p3mdl =  readRDS(paste0(RESULTS_DIR, "SingleLink_Mdl/Bayes_p3/", seed, "_BayesModel_p3.rda"))
    p3sum = summary(p3mdl$fit)
    results_single[rowcnt:(rowcnt+3),c('est', 'ts_se')] = p3sum$statistics[c('r[1]', 'r[2]', 'theta[1]', 'sigma'),c('Mean', 'Time-series SE')]
    results_single[rowcnt:(rowcnt+3),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p3sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'sigma'),]
    results_single[rowcnt:(rowcnt+3), c('neff')] = effectiveSize(p3mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'sigma')]
    results_single[rowcnt:(rowcnt+3), c('rhat')] = gelman.diag(p3mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'sigma'),1]
    results_single[rowcnt:(rowcnt+3), c('run_time')] = p3mdl$time 
    rowcnt = rowcnt + 4
    
    # process Bayes P4
    p4mdl =  readRDS(paste0(RESULTS_DIR, "SingleLink_Mdl/Bayes_p4/", seed, "_BayesModel_p4.rda"))
    p4sum = summary(p4mdl$fit)
    results_single[rowcnt:(rowcnt+3),c('est', 'ts_se')] = p4sum$statistics[c('r[1]', 'r[2]', 'theta', 'sigma'),c('Mean', 'Time-series SE')]
    results_single[rowcnt:(rowcnt+3),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p4sum$quantiles[c('r[1]', 'r[2]', 'theta', 'sigma'),]
    results_single[rowcnt:(rowcnt+3), c('neff')] = effectiveSize(p4mdl$fit)[c('r[1]', 'r[2]', 'theta', 'sigma')]
    results_single[rowcnt:(rowcnt+3), c('rhat')] = gelman.diag(p4mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta', 'sigma'),1]
    results_single[rowcnt:(rowcnt+3), c('run_time')] = p4mdl$time 
    rowcnt = rowcnt + 4
    
    # process Bayes P5
    p5mdl =  readRDS(paste0(RESULTS_DIR, "SingleLink_Mdl/Bayes_p5/", seed, "_BayesModel_p5.rda"))
    p5sum = summary(p5mdl$fit)
    results_single[rowcnt:(rowcnt+3),c('est', 'ts_se')] = p5sum$statistics[c('r[1]', 'r[2]', 'theta', 'sigma'),c('Mean', 'Time-series SE')]
    results_single[rowcnt:(rowcnt+3),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p5sum$quantiles[c('r[1]', 'r[2]', 'theta', 'sigma'),]
    results_single[rowcnt:(rowcnt+3), c('neff')] = effectiveSize(p5mdl$fit)[c('r[1]', 'r[2]', 'theta', 'sigma')]
    results_single[rowcnt:(rowcnt+3), c('rhat')] = gelman.diag(p5mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta', 'sigma'),1]
    results_single[rowcnt:(rowcnt+3), c('run_time')] = p5mdl$time 
    rowcnt = rowcnt + 4
    
}

################################################################################
# Double Link
################################################################################
results_double = data.frame(pose = rep('DoubleLink', nmodels*5*nfiles),
                            seed = rep(1:nfiles, each = nmodels*5),
                            model = rep(rep(MODELS, each = 5),nfiles),
                            parameter = rep(c('r1', 'r2', 'theta1', 'theta2', 'sigma'), nfiles*nmodels),
                            est = NA,
                            quantile = matrix(NA, nfiles*5*nmodels, 5),
                            neff = NA,
                            rhat = NA,
                            ts_se = NA,
                            run_time = NA)

files = list.files(path = paste0(RESULTS_DIR, "DoubleLink_Mdl/Bayes_p3"))
seeds = as.numeric(gsub("(.+?)(\\_.*)", "\\1", files))
seeds = seeds[order(seeds)]
rowcnt = 1
for(i in 1:nfiles){
    nlinks = 2
    seed = seeds[i]
    print(sprintf("Processing Double Link file %.0f of %.0f", i, nfiles))
    
    # process LS results
    lsmdl = readRDS(paste0(RESULTS_DIR, "DoubleLink_Mdl/LS/", seed, "_LSModel.rda"))
    est = c(lsmdl$r.hat, lsmdl$theta.hat, lsmdl$sigma.hat)
    int = lsmdl$intervals[1:(nlinks+2),]
    results_double[rowcnt:(rowcnt+4), c('est')] = est
    results_double[rowcnt:(rowcnt+3), c('quantile.1', 'quantile.5')] = int
    results_double[rowcnt:(rowcnt+4), c('run_time')] = lsmdl$time
    rowcnt = rowcnt + 5
    
    # process Bayes P1
    p1mdl =  readRDS(paste0(RESULTS_DIR, "DoubleLink_Mdl/Bayes_p1/", seed, "_BayesModel_p1.rda"))
    p1sum = summary(p1mdl$fit)
    results_double[rowcnt:(rowcnt+4),c('est', 'ts_se')] = p1sum$statistics[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),c('Mean', 'Time-series SE')]
    results_double[rowcnt:(rowcnt+4),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p1sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),]
    results_double[rowcnt:(rowcnt+4), c('neff')] = effectiveSize(p1mdl$fit)[c('r[1]', 'r[2]', 'theta[1]','theta[2]', 'sigma')]
    results_double[rowcnt:(rowcnt+4), c('rhat')] = gelman.diag(p1mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'), 1]
    results_double[rowcnt:(rowcnt+4), c('run_time')] = p1mdl$time
    rowcnt = rowcnt + 5
    
    # process Bayes P2
    p2mdl =  readRDS(paste0(RESULTS_DIR, "DoubleLink_Mdl/Bayes_p2/", seed, "_BayesModel_p2.rda"))
    p2sum = summary(p2mdl$fit)
    results_double[rowcnt:(rowcnt+4),c('est', 'ts_se')] = p2sum$statistics[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),c('Mean', 'Time-series SE')]
    results_double[rowcnt:(rowcnt+4),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p2sum$quantiles[c('r[1]', 'r[2]', 'theta[1]','theta[2]', 'sigma'),]
    results_double[rowcnt:(rowcnt+4), c('neff')] = effectiveSize(p2mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma')]
    results_double[rowcnt:(rowcnt+4), c('rhat')] = gelman.diag(p2mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),1]
    results_double[rowcnt:(rowcnt+4), c('run_time')] = p2mdl$time
    rowcnt = rowcnt + 5
    
    # process Bayes P3
    p3mdl =  readRDS(paste0(RESULTS_DIR, "DoubleLink_Mdl/Bayes_p3/", seed, "_BayesModel_p3.rda"))
    p3sum = summary(p3mdl$fit)
    results_double[rowcnt:(rowcnt+4),c('est', 'ts_se')] = p3sum$statistics[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),c('Mean', 'Time-series SE')]
    results_double[rowcnt:(rowcnt+4),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p3sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),]
    results_double[rowcnt:(rowcnt+4), c('neff')] = effectiveSize(p3mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma')]
    results_double[rowcnt:(rowcnt+4), c('rhat')] = gelman.diag(p3mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),1]
    results_double[rowcnt:(rowcnt+4), c('run_time')] = p3mdl$time
    rowcnt = rowcnt + 5
    
    # process Bayes P4
    p4mdl =  readRDS(paste0(RESULTS_DIR, "DoubleLink_Mdl/Bayes_p4/", seed, "_BayesModel_p4.rda"))
    p4sum = summary(p4mdl$fit)
    results_double[rowcnt:(rowcnt+4),c('est', 'ts_se')] = p4sum$statistics[c('r[1]', 'r[2]','theta[1]', 'theta[2]', 'sigma'),c('Mean', 'Time-series SE')]
    results_double[rowcnt:(rowcnt+4),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p4sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),]
    results_double[rowcnt:(rowcnt+4), c('neff')] = effectiveSize(p4mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]','sigma')]
    results_double[rowcnt:(rowcnt+4), c('rhat')] = gelman.diag(p4mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]','sigma'),1]
    results_double[rowcnt:(rowcnt+4), c('run_time')] = p4mdl$time 
    rowcnt = rowcnt + 5
    
    # process Bayes P5
    p5mdl =  readRDS(paste0(RESULTS_DIR, "DoubleLink_Mdl/Bayes_p5/", seed, "_BayesModel_p5.rda"))
    p5sum = summary(p5mdl$fit)
    results_double[rowcnt:(rowcnt+4),c('est', 'ts_se')] = p5sum$statistics[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),c('Mean', 'Time-series SE')]
    results_double[rowcnt:(rowcnt+4),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p5sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),]
    results_double[rowcnt:(rowcnt+4), c('neff')] = effectiveSize(p5mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma')]
    results_double[rowcnt:(rowcnt+4), c('rhat')] = gelman.diag(p5mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'sigma'),1]
    results_double[rowcnt:(rowcnt+4), c('run_time')] = p5mdl$time 
    rowcnt = rowcnt + 5
    
}

################################################################################
# Triple Link
################################################################################
results_triple = data.frame(pose = rep('TripleLink', nmodels*6*nfiles),
                            seed = rep(1:nfiles, each = nmodels*6),
                            model = rep(rep(MODELS, each = 6),nfiles),
                            parameter = rep(c('r1', 'r2', 'theta1', 'theta2', 'theta3', 'sigma'), nfiles*nmodels),
                            est = NA,
                            quantile = matrix(NA, nfiles*6*nmodels, 5),
                            neff = NA,
                            rhat = NA,
                            ts_se = NA,
                            run_time = NA)

files = list.files(path = paste0(RESULTS_DIR, "TripleLink_Mdl/Bayes_p3"))
seeds = as.numeric(gsub("(.+?)(\\_.*)", "\\1", files))
seeds = seeds[order(seeds)]
rowcnt = 1
for(i in 1:nfiles){
    nlinks = 3
    seed = seeds[i]
    print(sprintf("Processing Triple Link file %.0f of %.0f", i, nfiles))
    
    # process LS results
    lsmdl = readRDS(paste0(RESULTS_DIR, "TripleLink_Mdl/LS/", seed, "_LSModel.rda"))
    est = c(lsmdl$r.hat, lsmdl$theta.hat, lsmdl$sigma.hat)
    int = lsmdl$intervals[1:(nlinks+2),]
    results_triple[rowcnt:(rowcnt+5), c('est')] = est
    results_triple[rowcnt:(rowcnt+4), c('quantile.1', 'quantile.5')] = int
    results_triple[rowcnt:(rowcnt+5), c('run_time')] = lsmdl$time
    rowcnt = rowcnt + 6
    
    # process Bayes P1
    p1mdl =  readRDS(paste0(RESULTS_DIR, "TripleLink_Mdl/Bayes_p1/", seed, "_BayesModel_p1.rda"))
    p1sum = summary(p1mdl$fit)
    results_triple[rowcnt:(rowcnt+5),c('est', 'ts_se')] = p1sum$statistics[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),c('Mean', 'Time-series SE')]
    results_triple[rowcnt:(rowcnt+5),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p1sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),]
    results_triple[rowcnt:(rowcnt+5), c('neff')] = effectiveSize(p1mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma')]
    results_triple[rowcnt:(rowcnt+5), c('rhat')] = gelman.diag(p1mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'), 1]
    results_triple[rowcnt:(rowcnt+5), c('run_time')] = p1mdl$time
    rowcnt = rowcnt + 6
    
    # process Bayes P2
    p2mdl =  readRDS(paste0(RESULTS_DIR, "TripleLink_Mdl/Bayes_p2/", seed, "_BayesModel_p2.rda"))
    p2sum = summary(p2mdl$fit)
    results_triple[rowcnt:(rowcnt+5),c('est', 'ts_se')] = p2sum$statistics[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),c('Mean', 'Time-series SE')]
    results_triple[rowcnt:(rowcnt+5),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p2sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),]
    results_triple[rowcnt:(rowcnt+5), c('neff')] = effectiveSize(p2mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma')]
    results_triple[rowcnt:(rowcnt+5), c('rhat')] = gelman.diag(p2mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),1]
    results_triple[rowcnt:(rowcnt+5), c('run_time')] = p2mdl$time
    rowcnt = rowcnt + 6
    
    # process Bayes P2
    p3mdl =  readRDS(paste0(RESULTS_DIR, "TripleLink_Mdl/Bayes_p3/", seed, "_BayesModel_p3.rda"))
    p3sum = summary(p3mdl$fit)
    results_triple[rowcnt:(rowcnt+5),c('est', 'ts_se')] = p3sum$statistics[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),c('Mean', 'Time-series SE')]
    results_triple[rowcnt:(rowcnt+5),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p3sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),]
    results_triple[rowcnt:(rowcnt+5), c('neff')] = effectiveSize(p3mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma')]
    results_triple[rowcnt:(rowcnt+5), c('rhat')] = gelman.diag(p3mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),1]
    results_triple[rowcnt:(rowcnt+5), c('run_time')] = p3mdl$time
    rowcnt = rowcnt + 6
    
    # process Bayes P4
    p4mdl =  readRDS(paste0(RESULTS_DIR, "TripleLink_Mdl/Bayes_p4/", seed, "_BayesModel_p4.rda"))
    p4sum = summary(p4mdl$fit)
    results_triple[rowcnt:(rowcnt+5),c('est', 'ts_se')] = p4sum$statistics[c('r[1]', 'r[2]','theta[1]', 'theta[2]', 'theta[3]', 'sigma'),c('Mean', 'Time-series SE')]
    results_triple[rowcnt:(rowcnt+5),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p4sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),]
    results_triple[rowcnt:(rowcnt+5), c('neff')] = effectiveSize(p4mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]','sigma')]
    results_triple[rowcnt:(rowcnt+5), c('rhat')] = gelman.diag(p4mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]','sigma'),1]
    results_triple[rowcnt:(rowcnt+5), c('run_time')] = p4mdl$time 
    rowcnt = rowcnt + 6
    
    # process Bayes P5
    p5mdl =  readRDS(paste0(RESULTS_DIR, "TripleLink_Mdl/Bayes_p5/", seed, "_BayesModel_p5.rda"))
    p5sum = summary(p5mdl$fit)
    results_triple[rowcnt:(rowcnt+5),c('est', 'ts_se')] = p5sum$statistics[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),c('Mean', 'Time-series SE')]
    results_triple[rowcnt:(rowcnt+5),c('quantile.1', 'quantile.2','quantile.3','quantile.4','quantile.5')] = p5sum$quantiles[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),]
    results_triple[rowcnt:(rowcnt+5), c('neff')] = effectiveSize(p5mdl$fit)[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma')]
    results_triple[rowcnt:(rowcnt+5), c('rhat')] = gelman.diag(p5mdl$fit)$psrf[c('r[1]', 'r[2]', 'theta[1]', 'theta[2]', 'theta[3]', 'sigma'),1]
    results_triple[rowcnt:(rowcnt+5), c('run_time')] = p5mdl$time 
    rowcnt = rowcnt + 6
    
}

################################################################################
# Combine and save
################################################################################
print("Combing files")
results = rbind(results_single, results_double, results_triple)
print("Saving Results...")

saveRDS(results, file = paste0(RESULTS_DIR, 'results/results_', format(Sys.time(), "%d%m%Y"),".rda"))
write.table(results, file = paste0(RESULTS_DIR, 'results/results_', format(Sys.time(), "%d%m%Y"),".csv"), sep = ',', col.names = TRUE, row.names = FALSE)
print("Complete!")
