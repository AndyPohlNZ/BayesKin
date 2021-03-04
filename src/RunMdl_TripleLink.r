########################### BayesKin ###########################
## Bayesian estimates of rigid body kinematics
## Andy Pohl
## UofC - Faculty of Kinesiology
################################################################

##################### RunMdl_TripleLink.r  ###################
## Performs analysis comparing Bayesian to Least squares approaches for
## single link kinematic chains.
################################################################
rm(list = ls())

WORKING_DIR = "" # Replace with appropiate location of working directory on local PC
SAVE_MDL = TRUE # assume save

INIT_TYPE = 'random'

setwd(WORKING_DIR)

source('library.r') # Retrieve function library
set.seed(1)  # Set seed for reproducibility

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


nits = 1000 # number of iterations.
seeds = 1:nits
i = 1
while (length(list.files('./TripleLink_Mdl/Bayes_p3/')) < nits){
    files = list.files(path = paste0("./TripleLink_Mdl/Bayes_p3/"))
    processed_seeds = as.numeric(gsub("(.+?)(\\_.*)", "\\1", files))
    #set seed
    seed = seeds[i]
    
    if(!(seed %in% processed_seeds)){
        set.seed(seed)
        print(sprintf("Iteration %.0f of %.0f  -  seed = %.0f",length(list.files('./TripleLink_Mdl/Bayes_p3/')), nits, seed))
        # gen obs
        y = gen_obs(posture, sigma_true)
        
        #Compute LS solution.
        LS_result = LS_soln(y, links, 
                            init_type = INIT_TYPE,
                            inits = true_vals)
        
        print("LS Solution Complete")
        r_mdl_file_name = paste0('./TripleLink_Mdl/LS/', as.character(as.numeric(seed)))
        r_mdl_file_name = paste0(r_mdl_file_name, '_')
        r_mdl_file_name = paste0(r_mdl_file_name, 'LSModel')
        r_mdl_file_name = paste0(r_mdl_file_name, '.rda')
        print(sprintf("Saving file %s", r_mdl_file_name))
        saveRDS(LS_result, file = r_mdl_file_name)
        
        # Compute Bayes_P1
            Bayes_p1_result = Bayes_soln(y=y,
                                         links = links,
                                         mdlfile = './BayesModel_p1.jags',
                                         init_type = INIT_TYPE,
                                         true_vals = true_vals,
                                         ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))
        print("Bayes P1 Solution Complete")
        mdl_type = 'BayesModel_p1'
        r_mdl_file_name = paste0('./TripleLink_Mdl/Bayes_p1/', as.character(as.numeric(seed)))
        r_mdl_file_name = paste0(r_mdl_file_name, '_')
        r_mdl_file_name = paste0(r_mdl_file_name, mdl_type)
        r_mdl_file_name = paste0(r_mdl_file_name, '.rda')
        print(sprintf("Saving file %s", r_mdl_file_name))
        saveRDS(Bayes_p1_result, file = r_mdl_file_name)
        
        # Compute Bayes_P2
            Bayes_p2_result = Bayes_soln(y=y,
                                         links = links,
                                         mdlfile = './BayesModel_p2.jags',
                                         init_type = INIT_TYPE,
                                         true_vals = true_vals,
                                         ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))
        
        print("Bayes P2 Solution Complete")
        mdl_type = 'BayesModel_p2'
        r_mdl_file_name = paste0('./TripleLink_Mdl/Bayes_p2/', as.character(as.numeric(seed)))
        r_mdl_file_name = paste0(r_mdl_file_name, '_')
        r_mdl_file_name = paste0(r_mdl_file_name, mdl_type)
        r_mdl_file_name = paste0(r_mdl_file_name, '.rda')
        print(sprintf("Saving file %s", r_mdl_file_name))
        saveRDS(Bayes_p2_result, file = r_mdl_file_name)
        
        # Compute Bayes_P3
            Bayes_p3_result = Bayes_soln(y=y,
                                         links = links,
                                         mdlfile = './BayesModel_p3.jags',
                                         init_type = INIT_TYPE,
                                         true_vals = true_vals,
                                         ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))
        
        print("Bayes P3 Solution Complete")
        mdl_type = 'BayesModel_p3'
        r_mdl_file_name = paste0('./TripleLink_Mdl/Bayes_p3/', as.character(as.numeric(seed)))
        r_mdl_file_name = paste0(r_mdl_file_name, '_')
        r_mdl_file_name = paste0(r_mdl_file_name, mdl_type)
        r_mdl_file_name = paste0(r_mdl_file_name, '.rda')
        print(sprintf("Saving file %s", r_mdl_file_name))
        saveRDS(Bayes_p3_result, file = r_mdl_file_name)
    }
    i = i+1
    print('#######################################################################')
}


