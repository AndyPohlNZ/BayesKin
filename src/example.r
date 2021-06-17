########################### BayesKin ###########################
## Bayesian estimates of rigid body kinematics
## Andy Pohl
## UofC - Faculty of Kinesiology

## Written by: Andy Pohl
## UofC - Faculty of Kinesiology
## June-Dec 2020
## Revision 1: June 2021
################################################################

##################### example.r  ###################
## Performs analysis comparing Bayesian to Least squares approaches for
## 3 link kinematic chains.
################################################################
rm(list = ls())
WORKING_DIR = "" # Replace with the location of the BayesKin/src directory on the local PC.

setwd(WORKING_DIR)
INIT_TYPE = 'ls_est' # specify the type of initial values to be used for MCMC sampling - either 'true_vals', 'ls_est', 'random
source('library.r') # Retrieve function library
set.seed(1)  # Set seed for reproducibility


## 1) Specify parameters
n.links = 3                             # Number of links.
seg.length = c(0.45, 0.35, 0.25)        # Length of segment 
r_true = c(0.07, 0.03)                  # Origin location in m.
theta_true = c(-55, -110, -10)          # Rotation angle in deg. 
sigma_true = 1.5/1000                   # Measurement noise in m. 
true_vals = c(r_true = r_true, theta_true = theta_true, sigma_true = sigma_true)

# Generate posture
links = list(gen_link(seg.length = seg.length[1],
                      plate.center = 0.7*seg.length[1]),
             gen_link(seg.length = seg.length[2],
                      plate.center = 0.5*seg.length[2]),
             gen_link(seg.length = seg.length[3],
                      plate.center = 0.5*seg.length[3]))
posture = gen_posture(links, r = r_true, theta = theta_true)

# posture is a list with entries $r containing the location of the origin, $J containing the location of each joint
# where J[[1]] is the origin and J[[4]] is the termination of the kinematic chain. $alpha gives the actual position markers

# Generate observation
y = gen_obs(posture, sigma_true)
# y is a list of length nlinks with each list item containing a 2xnmarkers array of observed marker locations in the 
# global frame.

plot_system(posture, y=y) # Provides a visualisation of the pose with observed markers in orange and true positions in blue

## 2) Apply LS or Bayesian models
#Compute LS solution.
LS_result = LS_soln(y, links,
                    inits = c(r_true, theta_true, sigma_true),
                    init_type = 'true_vals')

LS_result
# The output of LS_result contains estimates of r and theta_i given by r.hat, theta.hat and $sigma.hat
# along with 95% confidence intervals for each parameter in $intervals.  $time denotes the computation time of the LS solution
# $value provides the value of the cost function at the LS estimate.

# Compute a Bayesian solution
Bayes_result = Bayes_soln(y=y,
                          links = links,
                          mdlfile = './BayesModel_p3.jags', # choose what model to use by specifying the appropiate .jags file here.
                          init_type = INIT_TYPE,
                          true_vals = true_vals,
                          ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))
# Bayes_result provides a list with $fit being a 'coda' object containing the MCMC samples for the fitted model
# $time provides the computational 
summary(Bayes_result$fit) # provides posterior means for each pose parameter along with their 95% credible interval (via 2.5 and 97.5 quantiles)
plot(Bayes_result$fit) # provides traceplots and density estimates for each parameter.
