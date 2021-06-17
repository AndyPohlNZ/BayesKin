########################### BayesKin ###########################
## Bayesian estimates of rigid body kinematics
## Andy Pohl
## UofC - Faculty of Kinesiology

## Written by: Andy Pohl
## UofC - Faculty of Kinesiology
## June-Dec 2020
## Revision 1: April 2021
################################################################

##################### generate_results.r  ###################
## processes_results table for a batch of simulations
#############################################################


# preliminaries 
rm(list = ls()) # clear workspace

library('dplyr')
library('latex2exp')
library('xtable')
library('coda')
WORKING_DIR = "" # Replace with the location of the BayesKin/src directory on the local PC.
setwd(WORKING_DIR)
RESULTS_DIR = paste0(WORKING_DIR, "/src/") # location of results directory

setwd(RESULTS_DIR)
# Load results file

FILENAME = "" # Filename of .rda result file produced by process_results.r
results = readRDS(FILENAME)
results = data.frame(results, stringsAsFactors = F)
results$model = factor(results$model, levels = c('LS', "Bayes_p1", "Bayes_p2", "Bayes_p3", "Bayes_p4", "Bayes_p5"))
results$parameter = factor(results$parameter, levels = c('r1', "r2", "theta1", "theta2", "theta3", "sigma"))
# Convert r/sigma to mm
results[results$parameter %in% c('r1', 'r2', 'sigma'),c('est', 'quantile.1', 'quantile.2', 'quantile.3', 'quantile.4', 'quantile.5')] = results[results$parameter %in% c('r1', 'r2', 'sigma'),c('est', 'quantile.1', 'quantile.2', 'quantile.3', 'quantile.4', 'quantile.5')]*1000

# Set constants
POSES = c("SingleLink", "DoubleLink", "TripleLink")
MODELS = c("LS", "Bayes_p1", "Bayes_p2", "Bayes_p3", "Bayes_p4", "Bayes_p5")
PARMS = c("r1", "r2", "theta1", "theta2", "theta3", "sigma")
PARMLABS = list(r1=TeX("$r_1 \\,  (mm)$"), r2=TeX("$r_2  \\, (mm)$"), theta1=TeX("$\\theta_1 \\, (^o)$"), 
                theta2=TeX("$\\theta_2 \\, (^o)$"), theta3=TeX("$\\theta_3\\, (^o) "), sigma=TeX("$\\sigma \\, (mm)$"))
MODELLABS = c("LS", "Bayes P1", "Bayes P2", "Bayes P3", "Bayes P4", "Bayes P5")

TRUEVALS = list(r1 = 0.07*1000, r2 = 0.03*1000, 
                theta1 = -55, theta2 = -110, theta3 = -10, sigma = 1.5/1000 *1000)

################################################################################
## FILTERING Poor MCMC convergence
################################################################################
# identify poor MCMC convergence
bad_seeds = results %>% group_by(pose) %>% filter(rhat>1.1) %>% select(pose, model, seed) %>% distinct()
bad_seeds %>% summarise(n_distinct(seed))

results = results[!((results$seed %in% bad_seeds$seed[bad_seeds$pose == 'SingleLink']) & (results$pose == 'SingleLink')),]
results = results[!((results$seed %in% bad_seeds$seed[bad_seeds$pose == 'DoubleLink']) & (results$pose == 'DoubleLink')),]
results = results[!((results$seed %in% bad_seeds$seed[bad_seeds$pose == 'TripleLink']) & (results$pose == 'TripleLink')),]

# filter to only 1000 iterations per pose
results_temps = results[results$pose =='SingleLink',]
temp_seeds = unique(results_temps$seed)
results_temps = results_temps[results_temps$seed %in% temp_seeds[1:1000],]

results_tempd = results[results$pose =='DoubleLink',]
temp_seeds = unique(results_tempd$seed)
results_tempd = results_tempd[results_tempd$seed %in% temp_seeds[1:1000],]

results_tempt = results[results$pose =='TripleLink',]
temp_seeds = unique(results_tempt$seed)
results_tempt = results_tempt[results_tempt$seed %in% temp_seeds[1:1000],]

results = rbind(results_temps, results_tempd, results_tempt)

results %>% group_by(pose) %>% summarise(n_distinct(seed))

################################################################################
## MCMC convergence summary
################################################################################

MCMC_summary = results %>% filter(model != "LS") %>% 
    group_by(pose, model, parameter) %>% 
    summarise( mean_neff = mean(rhat), sd_neff = sd(rhat))

xtable(MCMC_summary)


################################################################################
## Table 1: Bias, variance, mse and RMSE
################################################################################

tab1 = results %>% group_by(pose, model, parameter) %>% summarize(bias = mean(est), var = sd(est)^2)

parms = c("r1", "r2", "theta1", "theta2", "theta3", "sigma")
for(i in 1:length(parms)){
    tab1$bias[tab1$parameter == parms[i]] = tab1$bias[tab1$parameter == parms[i]] - TRUEVALS[[parms[i]]]
}
tab1$bias[tab1$parameter %in% c('r1', 'r2', 'sigma')] = tab1$bias[tab1$parameter %in% c('r1', 'r2', 'sigma')] 
tab1$mse = tab1$bias^2 + tab1$var
tab1$rmse = sqrt(tab1$mse)


print(xtable(tab1[,c('pose', 'model', 'parameter', 'rmse')], digits = 3), include.rowname=F)

################################################################################
## FIGURE 1 strip charts
################################################################################
deltadeg = 2
deltar = .02 *1000
axislimits = list(r1 = TRUEVALS[[1]] +  c(-deltar, deltar),
                  r2 = TRUEVALS[[2]] +  c(-deltar, deltar),
                  theta1 = TRUEVALS[[3]] + c(-deltadeg, deltadeg),
                  theta2 = TRUEVALS[[4]] + c(-deltadeg, deltadeg),
                  theta3 = TRUEVALS[[5]] + c(-deltadeg, deltadeg),
                  sigma = c(0, 0.01))
axislimits$sigma = axislimits$sigma*1000

colors = list(grey = rgb(97/255, 97/255,97/255, 0.2),
              orange = rgb(241/255, 136/255,5/255, 1),
              blue = rgb(18/255, 78/255,120/255,  0.2),
              red = rgb(162/255, 0/255,33/255, 0.2),
              green = rgb(76/255, 159/255, 112/255, 0.2),
              purple = rgb(217/255, 187/255, 249/255, 0.2),
              pink = rgb(217/255, 3/255, 104/255, 0.2))

point_cols = list(colors$pink, colors$blue, colors$green, colors$purple)

# Adjust the following to choose what is plotted in the nparmsx4 grid of plots.
MODELS = c("LS", "Bayes_p1", "Bayes_p2", "Bayes_p3")
MODELLABS = c("LS", "Bayes P1", "Bayes P2", "Bayes P3")
################################################################################
## Single Link
################################################################################
sl_results = results %>% filter(pose == 'SingleLink')
par(mfcol = c(4,4), mar = c(3, 1, 1,1), mgp = c(2,0.5,0), cex.lab = 1.2, cex.main = 1.2)
for(i in 1:length(MODELS)){
    m_results = sl_results %>% filter(model == MODELS[i])
    ps = unique(m_results$parameter)
    for(j in 1:length(ps)){
        est = m_results %>% filter(parameter == ps[j]) %>% select(est)
        if(j ==1){
            plot(x = est[,1], y = runif(length(est[,1]), 0, 0.1), 
                 xlim = get(as.character(ps[j]),axislimits), ylim = c(-0.03,0.25), 
                 pch = 16, col = point_cols[[i]],
                 xlab = get(as.character(ps[j]),PARMLABS),
                 ylab =NA, yaxt='n', bty = 'n',
                 main = MODELLABS[i])
        }else{
            plot(x = est[,1], y = runif(length(est[,1]), 0, 0.1), 
                 xlim = get(as.character(ps[j]),axislimits), ylim = c(-0.03,0.25), 
                 pch = 16, col =  point_cols[[i]],
                 xlab = get(as.character(ps[j]),PARMLABS),  bty = 'n',
                 ylab =NA, yaxt='n')
        }
        
        segments(x0=TRUEVALS[[ps[j]]], y0=-0.02, x1=TRUEVALS[[ps[j]]], y1=0.12, col = colors$orange, lwd=1.5)
    }
}

################################################################################
## Double Link
################################################################################
sl_results = results %>% filter(pose == 'DoubleLink')
par(mfcol = c(5,4), mar = c(3, 1, 1,1), mgp = c(2,0.5,0), cex.lab = 1.2, cex.main = 1.2)
for(i in 1:length(MODELS)){
    m_results = sl_results %>% filter(model == MODELS[i])
    ps = unique(m_results$parameter)
    for(j in 1:length(ps)){
        est = m_results %>% filter(parameter == ps[j]) %>% select(est)
        if(j ==1){
            plot(x = est[,1], y = runif(length(est[,1]), 0, 0.1), 
                 xlim = get(as.character(ps[j]),axislimits), ylim = c(-0.03,0.25), 
                 pch = 16, col = point_cols[[i]],
                 xlab = get(as.character(ps[j]),PARMLABS), ylab =NA, yaxt='n', bty = 'n',
                 main = MODELLABS[i])
        }else{
            plot(x = est[,1], y = runif(length(est[,1]), 0, 0.1), 
                 xlim = get(as.character(ps[j]),axislimits),ylim = c(-0.03,0.25), 
                 pch = 16, col = point_cols[[i]],
                 xlab = get(as.character(ps[j]),PARMLABS), ylab =NA, yaxt='n', bty = 'n')
        }
        
        segments(x0=TRUEVALS[[ps[j]]], y0=-0.02, x1=TRUEVALS[[ps[j]]], y1=0.12, col = colors$orange, lwd=1.5)
    }
}

################################################################################
## Triple Link
################################################################################
sl_results = results %>% filter(pose == 'TripleLink')
par(mfcol = c(6,4), mar = c(3, 1, 1,1), mgp = c(2,0.5,0), cex.lab = 1.2, cex.main = 1.2)
for(i in 1:length(MODELS)){
    m_results = sl_results %>% filter(model == MODELS[i])
    ps = unique(m_results$parameter)
    for(j in 1:length(ps)){
        est = m_results %>% filter(parameter == ps[j]) %>% select(est)
        if(j ==1){
            plot(x = est[,1], y = runif(length(est[,1]), 0, 0.1), 
                 xlim = get(as.character(ps[j]),axislimits), ylim = c(-0.03,0.25), 
                 pch = 16, col = point_cols[[i]],
                 xlab = get(as.character(ps[j]),PARMLABS), ylab =NA, yaxt='n', bty = 'n',
                 main = MODELLABS[i])
        }else{
            plot(x = est[,1], y = runif(length(est[,1]), 0, 0.1),
                 xlim = get(as.character(ps[j]),axislimits), ylim = c(-0.03,0.25), 
                 pch = 16, col = point_cols[[i]],
                 xlab = get(as.character(ps[j]),PARMLABS), ylab =NA, yaxt='n', bty = 'n')
        }
        segments(x0=TRUEVALS[[ps[j]]], y0=-0.02, x1=TRUEVALS[[ps[j]]], y1=0.12, col = colors$orange, lwd=1.5)
    }
}


################################################################################
## Table 2: EQULIIVANCE
################################################################################

## Prior 1
p1 = results %>% filter(model =='Bayes_p1' )
Bayes_perf = data.frame(seed = numeric(), pose = character(), model = character(), parameter = character(), 
                        equilivance_ls = numeric(), equilivance_true = numeric(),outperformed = numeric())
for(i in 1:length(POSES)){
    temp = p1 %>% filter(pose == POSES[i])
    ps = unique(temp$parameter)
    for(j in 1:length(ps)){
        print(sprintf("Computing Pose %s, parameter: %s", POSES[i], ps[j]))
        temp1 = temp %>% filter(parameter == ps[j])
        seeds = unique(temp1$seed)
        for(k in 1:length(seeds)){
            tmpseed = seeds[k]
            print(sprintf("    seed: %s", as.character(tmpseed)))
            
            ls_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'LS') %>% select(est)
            bayes_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'Bayes_p1') %>% select(est)
            
            # Equilivance of Bayes and LS
            if((temp$est[k] - (2* temp1$ts_se[k]) < ls_est)  &(temp1$est[k] + (2* temp1$ts_se[k]) > ls_est)){
                tmpcvgls = 1
            } else{
                tmpcvgls =0
            }
            
            # equilivance of Bayes and True Val
            if((temp$est[k] - (2* temp1$ts_se[k]) < TRUEVALS[[ps[j]]])  &(temp1$est[k] + (2* temp1$ts_se[k]) > TRUEVALS[[ps[j]]])){
                tmpcvgtrue = 1
            } else{
                tmpcvgtrue =0
            }
            
            
            tmpimprv = as.numeric((abs(bayes_est - TRUEVALS[[ps[j]]])) < (abs(ls_est  - TRUEVALS[[ps[j]]])))
            
            toappend = data.frame(seed = tmpseed, pose = POSES[i], model = "Bayes_p1", parameter = ps[j], 
                                  equilivance_ls = tmpcvgls, equilivance_true = tmpcvgtrue, outperformed = tmpimprv)
            Bayes_perf = rbind(Bayes_perf, toappend)
            
        }
        
    }
}

Bayes_perf_p1 = Bayes_perf %>% group_by(pose, model, parameter) %>% summarize(equilivance_ls = mean(equilivance_ls), equilivance_true = mean(equilivance_true),
                                                                              perf = mean(outperformed))

## Prior 2
p2 = results %>% filter(model =='Bayes_p2' )
Bayes_perf = data.frame(seed = numeric(), pose = character(), model = character(), parameter = character(), 
                        equilivance_ls = numeric(), equilivance_true = numeric(),outperformed = numeric())
for(i in 1:length(POSES)){
    temp = p2 %>% filter(pose == POSES[i])
    ps = unique(temp$parameter)
    for(j in 1:length(ps)){
        print(sprintf("Computing Pose %s, parameter: %s", POSES[i], ps[j]))
        temp1 = temp %>% filter(parameter == ps[j])
        seeds = unique(temp1$seed)
        for(k in 1:length(seeds)){
            tmpseed = seeds[k]
            print(sprintf("    seed: %s", as.character(tmpseed)))
            
            ls_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'LS') %>% select(est)
            bayes_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'Bayes_p2') %>% select(est)
            
            # Equilivance of Bayes and LS
            if((temp$est[k] - (2* temp1$ts_se[k]) < ls_est)  &(temp1$est[k] + (2* temp1$ts_se[k]) > ls_est)){
                tmpcvgls = 1
            } else{
                tmpcvgls =0
            }
            
            # equilivance of Bayes and True Val
            if((temp$est[k] - (2* temp1$ts_se[k]) < TRUEVALS[[ps[j]]])  &(temp1$est[k] + (2* temp1$ts_se[k]) > TRUEVALS[[ps[j]]])){
                tmpcvgtrue = 1
            } else{
                tmpcvgtrue =0
            }
            
            
            tmpimprv = as.numeric((abs(bayes_est - TRUEVALS[[ps[j]]])) < (abs(ls_est  - TRUEVALS[[ps[j]]])))
            
            toappend = data.frame(seed = tmpseed, pose = POSES[i], model = "Bayes_p2", parameter = ps[j], 
                                  equilivance_ls = tmpcvgls, equilivance_true = tmpcvgtrue, outperformed = tmpimprv)
            Bayes_perf = rbind(Bayes_perf, toappend)
            
        }
        
    }
}

Bayes_perf_p2 = Bayes_perf %>% group_by(pose, model, parameter) %>% summarize(equilivance_ls = mean(equilivance_ls), equilivance_true = mean(equilivance_true),
                                                                              perf = mean(outperformed))

## Prior 3
p3 = results %>% filter(model =='Bayes_p3' )
Bayes_perf = data.frame(seed = numeric(), pose = character(), model = character(), parameter = character(), 
                        equilivance_ls = numeric(), equilivance_true = numeric(),outperformed = numeric())
for(i in 1:length(POSES)){
    temp = p3 %>% filter(pose == POSES[i])
    ps = unique(temp$parameter)
    for(j in 1:length(ps)){
        print(sprintf("Computing Pose %s, parameter: %s", POSES[i], ps[j]))
        temp1 = temp %>% filter(parameter == ps[j])
        seeds = unique(temp1$seed)
        for(k in 1:length(seeds)){
            tmpseed = seeds[k]
            print(sprintf("    seed: %s", as.character(tmpseed)))
            
            ls_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'LS') %>% select(est)
            bayes_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'Bayes_p3') %>% select(est)
            
            # Equilivance of Bayes and LS
            if((temp$est[k] - (2* temp1$ts_se[k]) < ls_est)  &(temp1$est[k] + (2* temp1$ts_se[k]) > ls_est)){
                tmpcvgls = 1
            } else{
                tmpcvgls =0
            }
            
            # equilivance of Bayes and True Val
            if((temp$est[k] - (2* temp1$ts_se[k]) < TRUEVALS[[ps[j]]])  &(temp1$est[k] + (2* temp1$ts_se[k]) > TRUEVALS[[ps[j]]])){
                tmpcvgtrue = 1
            } else{
                tmpcvgtrue =0
            }
            
            
            tmpimprv = as.numeric((abs(bayes_est - TRUEVALS[[ps[j]]])) < (abs(ls_est  - TRUEVALS[[ps[j]]])))
            
            toappend = data.frame(seed = tmpseed, pose = POSES[i], model = "Bayes_p3", parameter = ps[j], 
                                  equilivance_ls = tmpcvgls, equilivance_true = tmpcvgtrue, outperformed = tmpimprv)
            Bayes_perf = rbind(Bayes_perf, toappend)
            
        }
        
    }
}

Bayes_perf_p3 = Bayes_perf %>% group_by(pose, model, parameter) %>% summarize(equilivance_ls = mean(equilivance_ls), equilivance_true = mean(equilivance_true),
                                                                              perf = mean(outperformed))

## Prior 4
p4 = results %>% filter(model =='Bayes_p4' )
Bayes_perf = data.frame(seed = numeric(), pose = character(), model = character(), parameter = character(), 
                        equilivance_ls = numeric(), equilivance_true = numeric(),outperformed = numeric())
for(i in 1:length(POSES)){
    temp = p4 %>% filter(pose == POSES[i])
    ps = unique(temp$parameter)
    for(j in 1:length(ps)){
        print(sprintf("Computing Pose %s, parameter: %s", POSES[i], ps[j]))
        temp1 = temp %>% filter(parameter == ps[j])
        seeds = unique(temp1$seed)
        for(k in 1:length(seeds)){
            tmpseed = seeds[k]
            print(sprintf("    seed: %s", as.character(tmpseed)))
            
            ls_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'LS') %>% select(est)
            bayes_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'Bayes_p4') %>% select(est)
            
            # Equilivance of Bayes and LS
            if((temp$est[k] - (2* temp1$ts_se[k]) < ls_est)  &(temp1$est[k] + (2* temp1$ts_se[k]) > ls_est)){
                tmpcvgls = 1
            } else{
                tmpcvgls =0
            }
            
            # equilivance of Bayes and True Val
            if((temp$est[k] - (2* temp1$ts_se[k]) < TRUEVALS[[ps[j]]])  &(temp1$est[k] + (2* temp1$ts_se[k]) > TRUEVALS[[ps[j]]])){
                tmpcvgtrue = 1
            } else{
                tmpcvgtrue =0
            }
            
            
            tmpimprv = as.numeric((abs(bayes_est - TRUEVALS[[ps[j]]])) < (abs(ls_est  - TRUEVALS[[ps[j]]])))
            
            toappend = data.frame(seed = tmpseed, pose = POSES[i], model = "Bayes_p4", parameter = ps[j], 
                                  equilivance_ls = tmpcvgls, equilivance_true = tmpcvgtrue, outperformed = tmpimprv)
            Bayes_perf = rbind(Bayes_perf, toappend)
            
        }
        
    }
}

Bayes_perf_p4 = Bayes_perf %>% group_by(pose, model, parameter) %>% summarize(equilivance_ls = mean(equilivance_ls), equilivance_true = mean(equilivance_true),
                                                                              perf = mean(outperformed))

## Prior 5
p5 = results %>% filter(model =='Bayes_p5' )
Bayes_perf = data.frame(seed = numeric(), pose = character(), model = character(), parameter = character(), 
                        equilivance_ls = numeric(), equilivance_true = numeric(),outperformed = numeric())
for(i in 1:length(POSES)){
    temp = p5 %>% filter(pose == POSES[i])
    ps = unique(temp$parameter)
    for(j in 1:length(ps)){
        print(sprintf("Computing Pose %s, parameter: %s", POSES[i], ps[j]))
        temp1 = temp %>% filter(parameter == ps[j])
        seeds = unique(temp1$seed)
        for(k in 1:length(seeds)){
            tmpseed = seeds[k]
            print(sprintf("    seed: %s", as.character(tmpseed)))
            
            ls_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'LS') %>% select(est)
            bayes_est = results %>% filter(seed == tmpseed, parameter == ps[j], pose == POSES[i], model == 'Bayes_p5') %>% select(est)
            
            # Equilivance of Bayes and LS
            if((temp$est[k] - (2* temp1$ts_se[k]) < ls_est)  &(temp1$est[k] + (2* temp1$ts_se[k]) > ls_est)){
                tmpcvgls = 1
            } else{
                tmpcvgls =0
            }
            
            # equilivance of Bayes and True Val
            if((temp$est[k] - (2* temp1$ts_se[k]) < TRUEVALS[[ps[j]]])  &(temp1$est[k] + (2* temp1$ts_se[k]) > TRUEVALS[[ps[j]]])){
                tmpcvgtrue = 1
            } else{
                tmpcvgtrue =0
            }
            
            
            tmpimprv = as.numeric((abs(bayes_est - TRUEVALS[[ps[j]]])) < (abs(ls_est  - TRUEVALS[[ps[j]]])))
            
            toappend = data.frame(seed = tmpseed, pose = POSES[i], model = "Bayes_p5", parameter = ps[j], 
                                  equilivance_ls = tmpcvgls, equilivance_true = tmpcvgtrue, outperformed = tmpimprv)
            Bayes_perf = rbind(Bayes_perf, toappend)
            
        }
        
    }
}

Bayes_perf_p5 = Bayes_perf %>% group_by(pose, model, parameter) %>% summarize(equilivance_ls = mean(equilivance_ls), equilivance_true = mean(equilivance_true),
                                                                              perf = mean(outperformed))
Bayes_perf = rbind(Bayes_perf_p1, Bayes_perf_p2, Bayes_perf_p3, Bayes_perf_p4, Bayes_perf_p5)

Bayes_perf


################################################################################
## Computational time.
################################################################################
comp_time = results %>% select(pose, seed, model, run_time) %>% group_by(pose, model) 
comp_time = unique(comp_time)
comp_time = comp_time %>% summarise(run_time_mean = mean(run_time), run_time_sd = sd(run_time))

#   temp = Bayes_perf %>% group_by(pose, model, parameter) %>% summarise(equilivance = mean(coverage), perf = mean(outperformed))
xtable(comp_time, digits = 2)


#####################################################################################
## Plot likelihoods
#####################################################################################
source('./library.r')

set.seed(4)

## 1) Specify parameters
n.links = 3               # Number of links.
seg.length = c(0.45, 0.35, 0.25)         # Length of segment specified in m. NB: Pataky pg 2.
r_true = c(0.07, 0.03)      # Origin location in m.
theta_true = c(-55, -110, -10)          # Rotation angle in deg. NB: specified to avoid 0/360 issue
# and mirroring solutions.
sigma_true = 1.5/1000       # Measurement noise in m. NB Noise specified as 1.5mm
# midpoint of range explored by Pataky et al.
true_vals = c(r_true = r_true, theta_true = theta_true, sigma_true = sigma_true)


# Generate posture
links = list(gen_link(seg.length = seg.length[1],
                      plate.center = 0.7*seg.length[1]),
             gen_link(seg.length = seg.length[2],
                      plate.center = 0.5*seg.length[2]),
             gen_link(seg.length = seg.length[3],
                      plate.center = 0.5*seg.length[3]))
posture = gen_posture(links, r = r_true, theta = theta_true)
y = gen_obs(posture, sigma_true)

LS_result = LS_soln(y, links, 
                    inits = true_vals,
                    init_type = 'random' )
inits = c(0.0446749866008758, 0.012961977859959, -40.1004006620497, 88.534386176616, 142.351554371417)
# lik for r
rs = seq(-1, 1, by = 0.01)
costs = matrix(NA, nrow = length(rs), ncol = length(rs))
for(i in 1:length(rs)){
    for(j in 1:length(rs)){
        costs[i,j] = cost(params = c(rs[i], rs[j], true_vals[3:5]), y =y, links = links)
    }
}

filled.contour(x=rs, y = rs, z=costs, 
               color = function(n) rev(hcl.colors(n, "Spectral")),
               main = 'Cost for R', xlab = 'r1', ylab = 'r2', nlevels = 30,
               plot.axes={points(true_vals[1], true_vals[2], pch = 3, col = 'green')
                   points(LS_result$r.hat[1], LS_result$r.hat[2], pch = 3, col = 'blue')
                   points(inits[1], inits[2], pch = 3, col = 'black')
                   axis(1, seq(-1, 1, by = 0.1))
                   axis(2, seq(-1, 1, by = 0.1))})

# lik for theta 1 vs 2
thetas = seq(-360, 360, by = 1)
costs = matrix(NA, nrow = length(thetas), ncol = length(thetas))
for(i in 1:length(thetas)){
    for(j in 1:length(thetas)){
        costs[i,j] = cost(params = c(true_vals[1:2], thetas[i], thetas[j], true_vals[5]), y =y, links = links)
    }
}

filled.contour(x=thetas, y = thetas, z=costs, 
               color = function(n) rev(hcl.colors(n, "Spectral")),
               main = 'Cost for theta', xlab = 'theta1', ylab = 'theta2', nlevels = 30,
               plot.axes={points(true_vals[3], true_vals[4], pch = 3, col = 'green')
                   points(LS_result$theta.hat[1], LS_result$theta.hat[2], pch = 3, col = 'blue')
                   points(inits[3], inits[4], pch = 3, col = 'black')
                   axis(1, seq(-360, 360, by = 45))
                   axis(2, seq(-360, 360, by = 45))})


# lik for theta 2 vs 3
thetasi = seq(-135, 135, length.out = 500)
thetasj = seq(-45, 270, length.out = 500)

costs = matrix(NA, nrow = length(thetasi), ncol = length(thetasj))
for(i in 1:length(thetasi)){
    for(j in 1:length(thetasj)){
        costs[i,j] = cost(params = c(true_vals[1:3], thetasi[i], thetasj[j]), y=y, links = links)
    }
}

filled.contour(x=thetasi, y = thetasj, z=costs, 
               color = function(n) rev(hcl.colors(n, "Spectral")),
               main = "LS Cost", xlab = TeX('$\\theta_2$'), ylab = TeX('$\\theta_3$'), nlevels = 30,
               plot.axes={points(true_vals[4], true_vals[5], pch = 19, col = rgb(241/255, 136/255,5/255, 1))
                   points(LS_result$theta.hat[2], LS_result$theta.hat[3], pch = 19, col = rgb(217/255, 3/255, 104/255, 1))
                   points(inits[4], inits[5], pch = 19, col = 'black')
                   axis(1, seq(-135, 135, by = 45))
                   axis(2, seq(-45, 270, by = 45))})



## lik for sigma
sigmas = seq(1e-3, 1, length.out = 10000)
lls = rep(NA, length(sigmas))
for(i in 1:length(sigmas)){
    lls[i] = loglikelihood(hat = c(true_vals[1:5], log(sigmas[i])), y=y, links = links)
}

plot(x = log(sigmas), y=lls, type = 'l', main = 'LL for sigma', xlab = 'log(sigma)')
abline(v=log(true_vals[6]), col = 'red')



plot_system(posture, y)
posture_hat = gen_posture(links = links, r = LS_result$r.hat, theta = LS_result$theta.hat)
plot_system(posture_hat, y, add=T)

