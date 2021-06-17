########################### BayesKin ###########################
## Bayesian estimates of rigid body kinematics
## Andy Pohl
## UofC - Faculty of Kinesiology

## Written by: Andy Pohl
## UofC - Faculty of Kinesiology
## June-Dec 2020
## Revision 1: May 2020
################################################################

##################### SupplementE.r  ###################
## Performs additional analysis examining the robustness of Bayesian models to
## increases in measurement noise
################################################################

## Preliminaries
rm(list = ls())

WORKING_DIR = "/media/data/ActiveProjects/P01_BayesianKinematics/PublishedCode/BayesKin/src/SupplementE" # Replace with appropriate location of working directory on local PC
SAVE_MDL = TRUE # assume save
INIT_TYPE = 'true_vals' # either 'true_vals', 'ls_est', 'random
NITS = 1000 # number of iterations.

setwd(WORKING_DIR)

# External libraries
library('dplyr')
library('latex2exp')

source('./library.r') # Retrieve function library

set.seed(1)  # Set seed for reproducibility

########################### Section A  #########################
## Specify simulation parameters
################################################################

## 1) Specify parameters
n.links = 1                 # Number of links.
seg.length = c(0.45)        # Length of segment specified in m. NB: Pataky pg 2.
r_true = c(0.07, 0.03)      # Origin location in m.
theta_true = c(-55)         # Rotation angle in deg. NB: specified to avoid 0/360 issue
# and mirroring solutions.
sigma_true = c(0.1, 1.5, 5, 10, 20,40)*(1/1000) # Measurement Noise in m       
links = list(gen_link(seg.length = seg.length[1],
                      plate.center = 0.7*seg.length[1])) # Link geometry


########################### Section B  #########################
## Run simulations 
    # run NITS iterations per sigma level
################################################################
results = matrix(NA, nrow = length(sigma_true)*NITS*5, ncol = 7)
row_idx =1
for(i in 1:length(sigma_true)){
    true_vals = list(r_true = r_true, theta_true = theta_true, sigma_true = sigma_true[i])
    
    # Generate posture
    posture = gen_posture(links, r = true_vals$r_true, theta = true_vals$theta_true)
    seeds = 1:NITS # Generate seed values for each iteration
    
    for(j in 1:NITS){
        seed = seeds[j]
        set.seed(seed)
        print(sprintf("--------------------Iteration %.0f of %.0f  -  seed = %.0f---------------------\n",j, NITS, seed))
        # gen obs
        y = gen_obs(posture, true_vals$sigma_true)
        #plot_system(posture = posture, y=y)
        print("Compute LS solution.")
        LS_result = LS_soln(y, links,
                            inits = c(true_vals$r_true, true_vals$theta_true, true_vals$sigma_true),
                            init_type = 'true_vals')
        results[row_idx,] = c(true_vals$sigma_true, 1, seed, LS_result$r.hat[1], LS_result$r.hat[2], LS_result$theta.hat[1],LS_result$sigma.hat)
        row_idx = row_idx+1
        
        Bayes1_result = Bayes_soln(y=y,
                                   links = links,
                                   mdlfile = './BayesModel_vague1.jags', # Vague 1
                                   init_type = INIT_TYPE,
                                   true_vals = c(true_vals$r_true,true_vals$theta_true, true_vals$sigma_true),
                                   ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))  
        Bayes1_result = summary(Bayes1_result$fit)
        results[row_idx,] = c(true_vals$sigma_true, 2, seed, Bayes1_result$statistics['r[1]',1], Bayes1_result$statistics['r[2]',1], Bayes1_result$statistics['theta',1],  Bayes1_result$statistics['sigma',1])
        row_idx = row_idx +1
        
        Bayes2_result = Bayes_soln(y=y,
                                   links = links,
                                   mdlfile = './BayesModel_vague2.jags', # Vague 2
                                   init_type = INIT_TYPE,
                                   true_vals = c(true_vals$r_true,true_vals$theta_true, true_vals$sigma_true),
                                   ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))  
        Bayes2_result = summary(Bayes2_result$fit)
        results[row_idx,] = c(true_vals$sigma_true, 3, seed, Bayes2_result$statistics['r[1]',1], Bayes2_result$statistics['r[2]',1], Bayes2_result$statistics['theta',1],  Bayes2_result$statistics['sigma',1])
        row_idx = row_idx +1
 
        Bayes4_result = Bayes_soln(y=y,
                                   links = links,
                                   mdlfile = './BayesModel_p3.jags', # Weakly Informative
                                   init_type = INIT_TYPE,
                                   true_vals = c(true_vals$r_true,true_vals$theta_true, true_vals$sigma_true),
                                   ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))  
        Bayes4_result = summary(Bayes4_result$fit)
        results[row_idx,] = c(true_vals$sigma_true, 4, seed, Bayes4_result$statistics['r[1]',1], Bayes4_result$statistics['r[2]',1], Bayes4_result$statistics['theta[1]',1],  Bayes4_result$statistics['sigma',1])
        row_idx = row_idx +1
        
        Bayes5_result = Bayes_soln(y=y,
                                   links = links,
                                   mdlfile = './BayesModel_p2.jags', # Informative
                                   init_type = INIT_TYPE,
                                   true_vals = c(true_vals$r_true,true_vals$theta_true, true_vals$sigma_true),
                                   ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))
        
        Bayes5_result = summary(Bayes5_result$fit)
        results[row_idx,] = c(true_vals$sigma_true, 5, seed, Bayes5_result$statistics['r[1]',1], Bayes5_result$statistics['r[2]',1], Bayes5_result$statistics['theta',1],Bayes5_result$statistics['sigma',1])
        row_idx = row_idx +1
        
    }
}

if(SAVE_MDL){
    saveRDS(results, "./DeltaVarResults.rds")
}
results = readRDS('./DeltaVarResults.rds')

########################### Section C  #########################
## Process results and generate figures
################################################################
models = c('LS', 'Bayes Vague 1', "Bayes Vague 2", 
           'Bayes Weakly Informative', 'Bayes Informative')

abs_error = matrix(NA, nrow = length(sigma_true)*length(models)*4, ncol = 6)
abs_error = data.frame(abs_error)
names(abs_error) = c('sigma', 'Model', 'Parameter', 'lwr', 'med', 'upr')

perf = matrix(NA, nrow = 3*length(sigma_true)*length(models), ncol = 4)
perf = data.frame(perf)
names(perf) = c('sigma', 'Model', 'Parameter', 'Perf')
row_idx = 1
row_idx2 = 1
for(i in 1:length(sigma_true)){
    sigma_result = results[results[,1]==sigma_true[i],]
    for(j in 1:length(models)){ # Summarize abs error of each model
        mdl_result = sigma_result[sigma_result[,2] == j,]
        abs_error[row_idx,1] = unique(mdl_result[,1])
        abs_error[row_idx,2] = unique(mdl_result[,2])
        abs_error[row_idx,3] = 'rx'
        abs_error[row_idx, 4:6] = quantile(abs(mdl_result[,4] - r_true[1]), probs = c(0.025, 0.5, 0.975))*1000
        row_idx = row_idx+1
        
        abs_error[row_idx,1] = unique(mdl_result[,1])
        abs_error[row_idx,2] = unique(mdl_result[,2])
        abs_error[row_idx,3] = 'ry'
        abs_error[row_idx, 4:6] = quantile(abs(mdl_result[,5] - r_true[2]), probs = c(0.025, 0.5, 0.975))*1000
        row_idx = row_idx+1
        
        abs_error[row_idx,1] = unique(mdl_result[,1])
        abs_error[row_idx,2] = unique(mdl_result[,2])
        abs_error[row_idx,3] = 'theta1'
        abs_error[row_idx, 4:6] = quantile(abs(mdl_result[,6] - theta_true[1]), probs = c(0.025, 0.5, 0.975))
        row_idx = row_idx+1

        abs_error[row_idx,1] = unique(mdl_result[,1])
        abs_error[row_idx,2] = unique(mdl_result[,2])
        abs_error[row_idx,3] = 'sigma'
        abs_error[row_idx, 4:6] = quantile(abs(mdl_result[,7] - sigma_true[i]), probs = c(0.025, 0.5, 0.975))*1000
        row_idx = row_idx +1
        
        if(j >=2){ # Compute performance of Bayes vs LS.
            perf[row_idx2,1] = unique(mdl_result[,1])
            perf[row_idx2,2] = unique(mdl_result[,2])
            perf[row_idx2,3] = 'rx'
            perf[row_idx2,4] = mean(abs(mdl_result[,4] - r_true[1]) < abs(sigma_result[sigma_result[,2] == 1,4] - r_true[1]))
            row_idx2 = row_idx2 +1
            
            perf[row_idx2,1] = unique(mdl_result[,1])
            perf[row_idx2,2] = unique(mdl_result[,2])
            perf[row_idx2,3] = 'ry'
            perf[row_idx2,4] = mean(abs(mdl_result[,5] - r_true[2]) < abs(sigma_result[sigma_result[,2] == 1,5] - r_true[2]))
            row_idx2 = row_idx2 +1
            
            perf[row_idx2,1] = unique(mdl_result[,1])
            perf[row_idx2,2] = unique(mdl_result[,2])
            perf[row_idx2,3] = 'theta1'
            perf[row_idx2,4] = mean(abs(mdl_result[,6] - theta_true[1]) < abs(sigma_result[sigma_result[,2] == 1,6] - theta_true[1]))
            row_idx2 = row_idx2 +1

            perf[row_idx2,1] = unique(mdl_result[,1])
            perf[row_idx2,2] = unique(mdl_result[,2])
            perf[row_idx2,3] = 'sigma'
            perf[row_idx2,4] = mean(abs(mdl_result[,7] - sigma_true[i]) < abs(sigma_result[sigma_result[,2] == 1,7] - sigma_true[i]))
            row_idx2 = row_idx2 +1
        }
    }
}
abs_error$sigma = abs_error$sigma*1000
perf$sigma = perf$sigma*1000

########################### Section D  #########################
## Generate figures and tables
################################################################

# Color scheme
colors = list(grey = rgb(97/255, 97/255,97/255, 0.2),
              orange = rgb(241/255, 136/255,5/255, 1),
              blue = rgb(18/255, 78/255,120/255,  0.2),
              light_blue = rgb(93/255, 173/255,226/255,  0.2),
              red = rgb(162/255, 0/255,33/255, 0.2),
              green = rgb(76/255, 159/255, 112/255, 0.2),
              purple = rgb(217/255, 187/255, 249/255, 0.2),
              pink = rgb(217/255, 3/255, 104/255, 0.2))


colors2 = list(grey = rgb(97/255, 97/255,97/255, 0.2),
               orange = rgb(241/255, 136/255,5/255, 1),
               blue = rgb(18/255, 78/255,120/255, 1),
               light_blue = rgb(93/255, 173/255,226/255,  1),
               red = rgb(162/255, 0/255,33/255, 1),
               green = rgb(76/255, 159/255, 112/255, 1),
               purple = rgb(217/255, 187/255, 249/255, 1),
               pink = rgb(217/255, 3/255, 104/255, 1))
point_cols = list(colors$pink, colors$blue,colors$light_blue, colors$purple,colors$green)
point_cols2 = list(colors2$pink, colors2$blue,colors2$light_blue, colors2$purple,colors2$green)

## -----------------------------------------------------------------------------
# Figure 1 'Line' plot of absolute error accross strain.
## -----------------------------------------------------------------------------

ylabs = c(TeX("Absolute Error $\\hat{r}_x \\, (mm)$"), TeX("Absolute Error $\\hat{r}_y \\, (mm)$"), TeX("Absolute Error $\\hat{\\theta}_1 \\, (deg)$"), TeX("Absolute Error $\\hat{\\sigma} \\, (mm)$"))
titles = c(TeX("$r_x$"), TeX("$r_y$"), TeX('$\\theta_1'), TeX("$\\sigma$"))
titles = c("(a)","(b)",'(c)',"(d)")

parameters = c('rx', 'ry', 'theta1', 'sigma')
ylims = list(c(0, 250), c(0,250), c(0, 60), c(0, 40))
# j = 1:4 parameters
# i = 1:5 models

par(mfrow=c(1,4))
for(j in 1:4){
    for(i in 1:5){
        parm_abs_error = abs_error %>% filter(Parameter ==parameters[j])

        mdl_abs_error = parm_abs_error %>% filter(Model ==i)
        jitter = (2-(5-i)) * 0.05*mdl_abs_error[,1]
        
        plot_col2 = point_cols2[[i]]
        plot_col1 = point_cols[[i]]
        if(i ==1){
            plot(mdl_abs_error$sigma + jitter, mdl_abs_error$med,type = 'b', pch = 16, col =plot_col2,
                 xlab = TeX("Noise Level: $\\sigma \\, (mm)$"), ylab = ylabs[j], bty='n',
                 xlim = c(-0.02, 50), ylim = ylims[[j]],
                 xaxt ='n')
            segments(x0=mdl_abs_error$sigma + jitter, y0 =  mdl_abs_error$lwr, x1 = mdl_abs_error$sigma + jitter, y1 = mdl_abs_error$upr, lty = 3, col = plot_col2)
            axis(side = 1, at = (sigma_true*1000), as.character(sigma_true*1000), las=2)
        }else{
            lines(mdl_abs_error$sigma + jitter, mdl_abs_error$med,type = 'b', pch = 16, col =plot_col2)
            segments(x0=mdl_abs_error$sigma + jitter, y0 =  mdl_abs_error$lwr, x1 = mdl_abs_error$sigma + jitter, y1 = mdl_abs_error$upr, lty = 3, col = plot_col2)
            
        }
        title(titles[j], font.main=1)
        
    }
}
legend(0, 40, legend = models[1:5], 
       col = c(point_cols2[[1]],point_cols2[[2]], point_cols2[[3]], point_cols2[[4]], point_cols2[[5]]),
       lty = 1, bty='n')


## -----------------------------------------------------------------------------
# Figure 2 'Line' plot of proportion of sims where Bayes had < error than LS
## -----------------------------------------------------------------------------

titles = c(TeX("(a) \t $r_x$"),TeX("(b) \t $r_y$"),TeX('(c) \t $\\theta_1$'),TeX("(d) \t $\\sigma$"))
par(mfrow=c(1,4))
for(j in 1:4){
    for(i in 2:5){
        parm_perf = perf %>% filter(Parameter ==parameters[j])
        mdl_parm_perf = parm_perf %>% filter(Model ==i)
        jitter = 0
        
        plot_col2 = point_cols2[[i]]
        if(i ==2){
            plot(mdl_parm_perf$sigma + jitter, mdl_parm_perf$Perf,type = 'b', pch = 16, col =plot_col2,
                 xlab = TeX("Noise Level: $\\sigma \\, (mm)$"), ylab = "Proportion of Simulations", bty='n',
                 xlim = c(-0.02, 50), ylim = c(0.4,1),
                 xaxt ='n')
            axis(side = 1, at = (sigma_true*1000), as.character(sigma_true*1000), las=2)
        }else{
            lines(mdl_parm_perf$sigma + jitter, mdl_parm_perf$Perf,type = 'b', pch = 16, col =plot_col2)
            
        }
        title(titles[j], main.font=1)
        
    }
}
legend(0, 0.5, legend = models[2:5], 
       col = c(point_cols2[[2]], point_cols2[[3]], point_cols2[[4]], point_cols2[[5]]),
       lty = 1, bty='n')

## -----------------------------------------------------------------------------
# Figure 3 Scatter Plot of LS estimator vs Bayes estimator
## -----------------------------------------------------------------------------
axis_lims = list(c(-0.2*1000, 0.5*1000), c(-0.5*1000, 0.2*1000), c(-150,100), c(0*1000, 0.12*1000))
true_vals_vec = c(r_true*1000, theta_true, sigma_true[6]*1000)
results_40 = results[results[,1]==40/1000,] # This can be changed to explore other noise levels
results_40[,c(4,5,7)] = results_40[,c(4,5,7)]*1000
ylabs = c("Vague 1", "Vague 2", "Weakly Informative", "Informative")
titles = c(TeX("$r_x \\, (mm)$"), TeX("$r_y \\, (mm)$"), TeX('$\\theta_1 \\, (deg)'), TeX("$\\sigma \\, (mm)$"))

par(mfrow=c(4,4), mai = c(0.6, 0.6, 0.0, 0))
for(i in 2:5){
    for(j in 1:4){
        ls_result = results_40[results_40[,2]==1, 3+j]
        mdl_result = results_40[results_40[,2]==i, 3+j]
        b_ls = sqrt((mdl_result-true_vals_vec[j])^2) - sqrt((ls_result - true_vals_vec[j])^2)

        b_ls = ifelse(b_ls<0, 1,2)
        cols = list(rgb(52/255, 152/255, 219/255, 0.2),
                    colors$pink)
        mycol = rep("", length(b_ls))
        for(ii in 1:length(b_ls)){
            mycol[ii] = cols[[b_ls[ii]]]    
        }
        if(j ==1 & i ==5){
            plot(ls_result, mdl_result, pch = 16, col = mycol,
                 #asp = 1,
                 xlim = axis_lims[[j]], ylim = axis_lims[[j]],
                 xlab = "LS est",
                 ylab = ylabs[i-1],
                 bty='n')
        }else if (j >1 & i ==5 ){
            plot(ls_result, mdl_result, pch = 16, col = mycol,
                 #asp = 1,
                 xlim = axis_lims[[j]], ylim = axis_lims[[j]],
                 xlab = "LS est",
                 ylab = NA,
                 bty='n')
        }else if(j==1 & i < 5){
            plot(ls_result, mdl_result, pch = 16, col = mycol,
                 #asp = 1,
                 xlim = axis_lims[[j]], ylim = axis_lims[[j]],
                 xlab = NA,
                 ylab = ylabs[i-1],
                 bty='n') 
        }else{
            plot(ls_result, mdl_result, pch = 16, col = mycol,
                 #asp = 1,
                 xlim = axis_lims[[j]], ylim = axis_lims[[j]],
                 xlab = NA, ylab = NA,
                 bty='n')
        }
        if(i ==2){
            title(titles[j], line =-1)
        }
        
        abline(0,1, col = colors$orange)
        points(true_vals_vec[j],true_vals_vec[j], pch=18, cex = 2,col = colors2$orange)
        
    }
}

########################### Section E  #########################
## Additional Figures not referenced in supplement
################################################################

# Scatter plot of error distribution for each model/noise level
ylabs = list(r1=TeX("$r_1 \\,  (mm)$"), r2=TeX("$r_2  \\, (mm)$"), theta1=TeX("$\\theta_1 \\, (^o)$"), 
             sigma=TeX("$\\sigma \\, (mm)$"))
ylims = list(c(-0.1,0.6), c(-0.6, 0.4), c(-180, 60), c(-0.03, 0.12))

par(mfrow = c(4, 6),mai = c(0, 0.6, 0.0, 0), oma=c(0,0,0,0))
for(j in 1:4){
    for(i in 1:6){
        # i= 1 # noise level (1:6)
        # j = 1 # parameter (1:4)
        true_val_tmp = c(r_true, theta_true, sigma_true[i])
        
        sigma_results = results[(results[,1]==sigma_true[i]), ]
        plot_col = rep("", nrow(sigma_results))
        for(ii in 1:nrow(sigma_results)){
            plot_col[ii] = point_cols[[sigma_results[ii,2]]]
        }
        
        if(i != 1){
            plot(sigma_results[,2]+rnorm(nrow(sigma_results), 0, 0.1), sigma_results[,3+j],
                 pch = 16, col = plot_col,
                 ylab = NA, yaxt = 'n',
                 ylim =ylims[[j]],
                 xlab = NA, bty='n', xaxt = 'n')
        }else{
            plot(sigma_results[,2]+rnorm(nrow(sigma_results), 0, 0.1), sigma_results[,3+j],
                 pch = 16, col = plot_col,
                 ylab = ylabs[[j]],
                 ylim =ylims[[j]],
                 xlab = NA, bty='n', xaxt = 'n')
        }
        abline(h = true_val_tmp[j], col = colors$orange)
        
        if(j ==1){
            title(TeX(sprintf("$\\sigma$ = %.1fmm$", c(sigma_true[i]*1000))), line = -2)
        }
    }
}
