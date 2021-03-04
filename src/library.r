########################### BayesKin ###########################
## Code for a comparision of Bayesian models for inverse kinematic problems
## with Maximum likelihood estimation

## Written by: Andy Pohl
## UofC - Faculty of Kinesiology
## June-Dec 2020
################################################################

########################### library.r  #########################
## Contains source functions for BayesKin package
################################################################
##### Preliminaries
library(coda)
library(numDeriv)
library(rjags)

# set number of parallel cores
options(mc.cores=4)

########################### Section A  #########################
## Miscellaneous helper functions
################################################################
rt_ls <- function(n, df, mu, a){
  # Generates a random sample from the location scale student-t distribution with mean mu, degrees of freedom
  # df and scale a
  
  return(rt(n,df)*a + mu)
}

mod_180 = function(angles){
  # Rescales an angle on the (-inf, inf) scale to the interval [-180, 180)
  temp = sapply(angles, function(ang){
    if(abs(ang)<=180){
      return(ang)
    }else if(ang >180){
      remainder = ang-180
      return(-180 + remainder)
    }else if (ang< -180){
      remainder = ang + 180
      return(180 + remainder)
    }
  })
  return(temp)
}
########################### Section B  #########################
## Links, posture and observations
################################################################

gen_link = function(seg.length,
                    plate.center,
                    marker.plate.dim = c(0.1, 0.05)){
  # Creates a rigid link of length seg.length with a marker plate at plate.center
  # with dimensions plate.placement.length*seg.length and markers at each corner
  # of the marker plate.
  plate.loc = c(plate.center,0)
  
  
  markers = cbind(c(plate.loc[1]-marker.plate.dim[1]/2, marker.plate.dim[2]/2),
                  c(plate.loc[1]+marker.plate.dim[1]/2, marker.plate.dim[2]/2),
                  c(plate.loc[1]+marker.plate.dim[1]/2, -marker.plate.dim[2]/2),
                  c(plate.loc[1]-marker.plate.dim[1]/2, -marker.plate.dim[2]/2))
  
  link = list(length = seg.length,
              x = markers)
  return(link)
  
}

gen_posture = function(links,
                       r,
                       theta,
                       degrees = TRUE){
  # Applies rigid body transformations translation r and rotation by theta to a
  # collection of links.  Note that links must be specified in order (link 1 - 
  # link 3) and the length of theta should be the same as the number of links
  
  if(degrees){
    theta = theta * pi/180
  }
  
  nlinks = length(links)
  
  #  Generate rotation matrix for eac link
  Gam = vector('list', nlinks)
  for(i in 1:nlinks){
    Gam[[i]] = matrix(c(cos(theta[i]), sin(theta[i]), 
                        -sin(theta[i]), cos(theta[i])), 2, 2)
  }
  
  # Compute the locations of joint centers for each link in the kinematic chain
  J = vector('list', nlinks+1)
  J[[1]] = matrix(r, 2, 1)
  for(i in 2:(nlinks+1)){
    J[[i]] =  J[[i-1]] + Gam[[i-1]] %*% c(links[[i-1]]$length, 0)
  }
  
  # Compute actual position of markers alpha in the global coordinate system
  alpha = vector(mode = 'list', length = nlinks)
  for(i in 1:nlinks){
    nmarkers = ncol(links[[i]]$x)
    alpha[[i]] = matrix(NA, 2, nmarkers)
    for(j in 1:nmarkers){
      x_ij = links[[i]]$x[,j]
      alpha[[i]][,j] = J[[i]] + Gam[[i]]%*%x_ij
    }
  }
  
  # Return posture object (list with origin r, joint position J and expected marker position alpha)
  return(list(r = r, J = J, alpha = alpha))
}

gen_obs = function(posture, sigma){
  # Generates an observation of a given posture subject to random 0 mean gaussain
  # noise with standard deviation sigma
  nlinks = length(posture$alpha)
  
  y = vector(mode = 'list', length = nlinks)
  for(i in 1:nlinks){
    nmarkers = ncol(links[[i]]$x)
    y[[i]] = matrix(NA, ncol = nmarkers, nrow = 2)
    for(j in 1:nmarkers){
      y[[i]][1, j] = rnorm(1, mean = posture$alpha[[i]][1, j], sd= sigma)
      y[[i]][2, j] = rnorm(1, mean = posture$alpha[[i]][2, j], sd= sigma)
    }
  }
  return(y)
}

plot_system= function(posture, y=NA, ...){
  # Plots the Kinematic chain and observed markers y
  plot(NULL, NULL,
       xlim = c(0, 1), ylim = c(-1, 0.1),
       xlab = 'x', ylab = 'y')
  points(t(do.call(cbind,posture$J)), type = "o", pch = 20,...) # mu_tmp
  orange = rgb(255/255, 140/255, 0, 0.5)
  blue = rgb(0, 0, 1, 0.5)
  points(t(do.call(cbind,posture$alpha)), pch = 20, col = orange)
  if(!is.na(y)){
    points(t(do.call(cbind,y)), pch = 20, col = blue)
  }
}

########################### Section C  #########################
## Maximum likelihood/Least Squares solutions.
################################################################

cost = function(params, y, links){
  # The cost function (equation ... in paper) given a vector
  # params = vector of parameters (rx, ry, theta1, theta2, theta3 ), observations 
  # y and known positions of markers within the LCS of each link links.
  nlinks = length(links)
  r.hat = params[1:2]
  theta.hat = params[3:(nlinks+2)]
  
  obs.hat = gen_posture(links, r = r.hat, theta = theta.hat)
  y.hat = obs.hat$alpha
  
  cost = 0
  for(i in 1:nlinks){
    nmarkers = ncol(y.hat[[i]])
    for(j in 1:nmarkers){
      diff = y.hat[[i]][,j] - y[[i]][,j]
      cost = cost+ (sqrt(sum(diff*diff)))
    }
  }
  return(cost)
}


costgrad_analytical = function(x, y, links){
  # Generates the analytical gradient of the cost function.  Originally generated 
  # via MATLAB and the symbolic toolbox.
  # x is parameter vector (r1, r2, theta1, theta2, theta3)
  # y[[nlinks]] of observations
  # links contains the known LCS information of each link
  
  nlinks = length(links)
  r1 = x[1]
  r2 = x[2]
  x111 = links[[1]]$x[1,1]
  x112 = links[[1]]$x[2,1]
  x121 = links[[1]]$x[1,2]
  x122 = links[[1]]$x[2,2]
  x131 = links[[1]]$x[1,3]
  x132 = links[[1]]$x[2,3]
  x141 = links[[1]]$x[1,4]
  x142 = links[[1]]$x[2,4]
  
  y111 = y[[1]][1,1]
  y112 = y[[1]][2,1]
  y121 = y[[1]][1,2]
  y122 = y[[1]][2,2]
  y131 = y[[1]][1,3]
  y132 = y[[1]][2,3]
  y141 = y[[1]][1,4]
  y142 = y[[1]][2,4]
  
  
  if(nlinks ==1){
    theta1 = x[3]
    gradient = rep(NA, 3)
    for(i in 1:length(x)){
        gradtxtfile = sprintf("./AnalyticalCostGrad/SingleLink/CostGradientChar%.0f.txt", i)
        gradtxt = read.table(file = gradtxtfile, sep=',' , stringsAsFactors = FALSE)
        if((dim(gradtxt)[1] != 1) & (dim(gradtxt)[2] !=1)){stop("Error in reading gradient text file!")}
        gradient[i] = eval(parse(text = gradtxt[1,1]))
      }
    }else if(nlinks ==2){
    theta1 = x[3]
    theta2 = x[4]
    
    L1 = links[[1]]$length
    
    x211 = links[[2]]$x[1,1]
    x212 = links[[2]]$x[2,1]
    x221 = links[[2]]$x[1,2]
    x222 = links[[2]]$x[2,2]
    x231 = links[[2]]$x[1,3]
    x232 = links[[2]]$x[2,3]
    x241 = links[[2]]$x[1,4]
    x242 = links[[2]]$x[2,4]
    
    y211 = y[[2]][1,1]
    y212 = y[[2]][2,1]
    y221 = y[[2]][1,2]
    y222 = y[[2]][2,2]
    y231 = y[[2]][1,3]
    y232 = y[[2]][2,3]
    y241 = y[[2]][1,4]
    y242 = y[[2]][2,4]
    
    gradient = rep(NA, 4)
    for(i in 1:length(x)){
      gradtxtfile = sprintf("./AnalyticalCostGrad/DoubleLink/CostGradientChar%.0f.txt", i)
      gradtxt = read.table(file = gradtxtfile, sep=',' , stringsAsFactors = FALSE)
      if((dim(gradtxt)[1] != 1) & (dim(gradtxt)[2] !=1)){stop("Error in reading gradient text file!")}
      gradient[i] = eval(parse(text = gradtxt[1,1]))
    }
    
  }else if(nlinks ==3){
    theta1 = x[3]
    theta2 = x[4]
    theta3 = x[5]
    
    L1 = links[[1]]$length
    L2 = links[[2]]$length
    
    x211 = links[[2]]$x[1,1]
    x212 = links[[2]]$x[2,1]
    x221 = links[[2]]$x[1,2]
    x222 = links[[2]]$x[2,2]
    x231 = links[[2]]$x[1,3]
    x232 = links[[2]]$x[2,3]
    x241 = links[[2]]$x[1,4]
    x242 = links[[2]]$x[2,4]
    
    y211 = y[[2]][1,1]
    y212 = y[[2]][2,1]
    y221 = y[[2]][1,2]
    y222 = y[[2]][2,2]
    y231 = y[[2]][1,3]
    y232 = y[[2]][2,3]
    y241 = y[[2]][1,4]
    y242 = y[[2]][2,4]
    
    x311 = links[[3]]$x[1,1]
    x312 = links[[3]]$x[2,1]
    x321 = links[[3]]$x[1,2]
    x322 = links[[3]]$x[2,2]
    x331 = links[[3]]$x[1,3]
    x332 = links[[3]]$x[2,3]
    x341 = links[[3]]$x[1,4]
    x342 = links[[3]]$x[2,4]
    
    y311 = y[[3]][1,1]
    y312 = y[[3]][2,1]
    y321 = y[[3]][1,2]
    y322 = y[[3]][2,2]
    y331 = y[[3]][1,3]
    y332 = y[[3]][2,3]
    y341 = y[[3]][1,4]
    y342 = y[[3]][2,4]
    gradient = rep(NA, 5)
    for(i in 1:length(x)){
      gradtxtfile = sprintf("./AnalyticalCostGrad/TripleLink/CostGradientChar%.0f.txt", i)
      gradtxt = read.table(file = gradtxtfile, sep=',' , stringsAsFactors = FALSE)
      if((dim(gradtxt)[1] != 1) & (dim(gradtxt)[2] !=1)){stop("Error in reading gradient text file!")}
      gradient[i] = eval(parse(text = gradtxt[1,1]))
    }
  }
  return(gradient)
}

hess_analytical = function(x, y, links){
  # Generates the analytical hessian of the cost function.  Originally generated 
  # via MATLAB and the symbolic toolbox.
  # x is parameter vector (r1, r2, theta1, theta2, theta3)
  # y[[nlinks]] of observations
  # links contains the known LCS information of each link
  
  nlinks = length(links)
  r1 = x[1]
  r2 = x[2]
  lsigma = x[length(x)]
  
  x111 = links[[1]]$x[1,1]
  x112 = links[[1]]$x[2,1]
  x121 = links[[1]]$x[1,2]
  x122 = links[[1]]$x[2,2]
  x131 = links[[1]]$x[1,3]
  x132 = links[[1]]$x[2,3]
  x141 = links[[1]]$x[1,4]
  x142 = links[[1]]$x[2,4]
  
  y111 = y[[1]][1,1]
  y112 = y[[1]][2,1]
  y121 = y[[1]][1,2]
  y122 = y[[1]][2,2]
  y131 = y[[1]][1,3]
  y132 = y[[1]][2,3]
  y141 = y[[1]][1,4]
  y142 = y[[1]][2,4]
  
  
  if(nlinks ==1){
    theta1 = x[3]
    hess = matrix(NA, nrow = length(x), ncol = length(x))
    for(i in 1:length(x)){
      for(j in 1:length(x)){
        hesstxtfile = sprintf("./AnalyticalHessian/SingleLink/hesschar%.0f_%.0f.txt", i,j)
        hesstxt = read.table(file = hesstxtfile, sep=',' , stringsAsFactors = FALSE)
        if((dim(hesstxt)[1] != 1) & (dim(hesstxt)[2] !=1)){stop("Error in reading hessian text file!")}
        hess[i,j] = eval(parse(text = hesstxt[1,1]))
      }
    }
    
  }else if(nlinks ==2){
    theta1 = x[3]
    theta2 = x[4]
    L1 = links[[1]]$length
    
    x211 = links[[2]]$x[1,1]
    x212 = links[[2]]$x[2,1]
    x221 = links[[2]]$x[1,2]
    x222 = links[[2]]$x[2,2]
    x231 = links[[2]]$x[1,3]
    x232 = links[[2]]$x[2,3]
    x241 = links[[2]]$x[1,4]
    x242 = links[[2]]$x[2,4]
    
    y211 = y[[2]][1,1]
    y212 = y[[2]][2,1]
    y221 = y[[2]][1,2]
    y222 = y[[2]][2,2]
    y231 = y[[2]][1,3]
    y232 = y[[2]][2,3]
    y241 = y[[2]][1,4]
    y242 = y[[2]][2,4]
    
    hess = matrix(NA, nrow = length(x), ncol = length(x))
    for(i in 1:length(x)){
      for(j in 1:length(x)){
        hesstxtfile = sprintf("./AnalyticalHessian/DoubleLink/hesschar%.0f_%.0f.txt", i,j)
        hesstxt = read.table(file = hesstxtfile, sep=',' , stringsAsFactors = FALSE)
        if((dim(hesstxt)[1] != 1) & (dim(hesstxt)[2] !=1)){stop("Error in reading hessian text file!")}
        hess[i,j] = eval(parse(text = hesstxt[1,1]))
      }
    }
  }else if(nlinks ==3){
    theta1 = x[3]
    theta2 = x[4]
    theta3 = x[5]
    
    L1 = links[[1]]$length
    L2 = links[[2]]$length
    
    x211 = links[[2]]$x[1,1]
    x212 = links[[2]]$x[2,1]
    x221 = links[[2]]$x[1,2]
    x222 = links[[2]]$x[2,2]
    x231 = links[[2]]$x[1,3]
    x232 = links[[2]]$x[2,3]
    x241 = links[[2]]$x[1,4]
    x242 = links[[2]]$x[2,4]
    
    y211 = y[[2]][1,1]
    y212 = y[[2]][2,1]
    y221 = y[[2]][1,2]
    y222 = y[[2]][2,2]
    y231 = y[[2]][1,3]
    y232 = y[[2]][2,3]
    y241 = y[[2]][1,4]
    y242 = y[[2]][2,4]
    
    x311 = links[[3]]$x[1,1]
    x312 = links[[3]]$x[2,1]
    x321 = links[[3]]$x[1,2]
    x322 = links[[3]]$x[2,2]
    x331 = links[[3]]$x[1,3]
    x332 = links[[3]]$x[2,3]
    x341 = links[[3]]$x[1,4]
    x342 = links[[3]]$x[2,4]
    
    y311 = y[[3]][1,1]
    y312 = y[[3]][2,1]
    y321 = y[[3]][1,2]
    y322 = y[[3]][2,2]
    y331 = y[[3]][1,3]
    y332 = y[[3]][2,3]
    y341 = y[[3]][1,4]
    y342 = y[[3]][2,4]
    
    hess = matrix(NA, nrow = length(x), ncol = length(x))
    for(i in 1:length(x)){
      for(j in 1:length(x)){
        hesstxtfile = sprintf("./AnalyticalHessian/TripleLink/hesschar%.0f_%.0f.txt", i,j)
        hesstxt = read.table(file = hesstxtfile, sep=',' , stringsAsFactors = FALSE)
        if((dim(hesstxt)[1] != 1) & (dim(hesstxt)[2] !=1)){stop("Error in reading hessian text file!")}
        hess[i,j] = eval(parse(text = hesstxt[1,1]))
      }
    }
  }
  return(hess)
}


LS_soln = function(y, links,
                   init_type = 'random',
                   inits = NA){
  # Computes the least squares solution given observations y.
  # y = observaitons
  # links = known LCS geometry
  # init_type = type of initial values to specify (random = randomly generated 
  #             true_vals = true values of the simulation or speicified = custom
  #             specification of initial values.)
  # inits = specified initis if true_vals is specified.
  nlinks = length(links)
  
  if(init_type =='random'){
    inits.r = c(runif(1, -0.1, 0.1),
                runif(1, -0.1, 0.1))
    inits.theta = runif(nlinks, -180, 180)
    print(paste('Computing LS solution with initial Values set to the RANDOM VALUES:', paste(c(inits.r, inits.theta), collapse = ', ' )))
    
    
  }else if(init_type == 'true_vals'){
    if (is.na(inits)){stop("If init_type is not random then inits must be specified")}
    print('Computing LS solution with initial Values set to the TRUE VALUES')
    inits.r = inits[1:2]
    inits.theta = inits[3:(nlinks +2)]
  }else if (init_type == "specified"){
    if (is.na(inits)){stop("If init_type is not random then inits must be specified")}
    print('Computing LS solution with initial Values set to the SPECIFIED VALUES')
    inits.r = inits[1:2]
    inits.theta = inits[3:(nlinks +2)]
  }
  inits = c(inits.r, inits.theta)
  
  # run optimisation using BFGS method.
  t1 = Sys.time()
  temp = optim(par=inits,fn = cost, gr = costgrad_analytical,
               y = y, links = links,
               method = "BFGS", control = list(trace=FALSE, maxit = 10000, reltol = 1e-20))
  t2 = Sys.time()
  return_list = list(r.hat = temp$par[1:2],
                     theta.hat = temp$par[3:(nlinks+2)],
                     value = temp$val, time = as.numeric(t2-t1))
  
  # generate residuals and use for corresponding estimate of sigma_hat
  nmarkers = 0
  obs.hat = gen_posture(links, return_list$r.hat, return_list$theta.hat)
  residual_sum = 0
  for(i in 1:nlinks){
    nmarkers = ncol(y[[i]])
    for(j in 1:nmarkers){
      residual_sum = residual_sum +
        (y[[i]][1, j]- obs.hat$alpha[[i]][1,j])^2 + # xresidual
        (y[[i]][2, j]- obs.hat$alpha[[i]][2, j])^2 # yresidual
      nmarkers = nmarkers +1
    }
    
  }
  p = length(return_list$r.hat) + length(return_list$theta.hat)
  return_list$sigma.hat = sqrt(residual_sum / (2*nmarkers)) 
  
  # Generate CI estimate via hessian matrix
  hess = hess_analytical(x = c(return_list$r.hat, return_list$theta.hat, log(return_list$sigma.hat)),
                         links = links,
                         y = y)
  fisher_info = solve(-hess)
  prop_sigma = diag(diag(sqrt(diag(fisher_info))))
  intervals = matrix(NA, nrow=nlinks+3, ncol =2)
  intervals[,1] = c(return_list$r.hat, return_list$theta.hat, return_list$sigma.hat) - 1.96*prop_sigma
  intervals[,2] = c(return_list$r.hat, return_list$theta.hat, return_list$sigma.hat) + 1.96*prop_sigma
  intervals[nrow(intervals), ] = exp(intervals[nrow(intervals),])
  return_list$intervals = intervals
  
  return(return_list)
}

########################### Section D  #########################
## Compute Bayesian solution.
################################################################

Bayes_soln = function(y, links, mdlfile,
                      true_vals,
                      init_type = 'random',
                      ls_est = NA,
                      ...){
  # Applies the Bayesian model to observations y and known data links.
  # y = observaitons
  # links = known LCS geometry
  # mdlfile = specified jags model file
  # true_vals = true values of the specified simulation.
  # init_type = type of initial values to specify (random = randomly generated 
  #             true_vals = true values of the simulation or speicified = custom
  #             specification of initial values.)
  # inits = specified initis if true_vals is specified.
  
  
  nlinks = length(links)
  nmarkers = sapply(links, function(x){ncol(x$x)})
  M = 2 # dimension
  seg_lengths = sapply(links, function(x){x$length})
  x = lapply(links, function(x){x$x})
  
  # initilise x and y data for jags models
  xjags = yjags = array(NA, dim = c(nlinks, max(nmarkers), 2))
  for(i in 1:nlinks){
    xjags[i,,] = t(x[[i]])
    yjags[i,,] = t(y[[i]])
  }
  
  # input data for jags model
  jags_input = list(nlinks = nlinks,
                    nmarkers = nmarkers,
                    ndim = 2,
                    x = xjags,
                    y = yjags,
                    leng = cbind(seg_lengths, rep(0,nlinks)),
                    true_vals = true_vals)
  
  # Generate initial values
  if(init_type =='true_vals'){
    if(grepl('p3', mdlfile)){
      print('Computing Bayes P3 solution with initial Values set to the TRUE VALUES')
      inits = list(hip = theta_true[1],
                   knee = theta_true[2] - theta_true[1],
                   ankle = -90 -theta_true[2] + theta_true[3],
                   r = r_true,
                   sigma_mm = sigma_true*1000)
    }else if(grepl('p2', mdlfile)){
      print('Computing Bayes P2 solution with initial Values set to the TRUE VALUES')
      inits = list(thetar = theta_true * pi/180,
                   r = r_true,
                   tau = 1/sigma_true/sigma_true)
    }else{
      print('Computing Bayes P1 solution with initial Values set to the TRUE VALUES')
      inits = list(thetar = theta_true * pi/180,
                   r = r_true,
                   sigma = sigma_true)
    }
  }else if(init_type =='specified'){
    
    init_theta = true_vals[3:nlinks]
    init_r = true_vals[1:2]
    init_sigma = true_vals[nlinks +3]
    
    if(grepl('p3', mdlfile)){
      
      print('Computing Bayes P3 solution with initial Values set to the SPECIFIED VALUES')
      inits = list(hip = init_theta[1],
                   knee = init_theta[2] - init_theta[1],
                   ankle = -90 -init_theta[2] + init_theta[3],
                   r = init_r,
                   sigma_mm = init_sigma*1000)
    }else if(grepl('p2', mdlfile)){
      print('Computing Bayes P2 solution with initial Values set to the SPECIFIED VALUES')
      inits = list(thetar = init_theta * pi/180,
                   r = init_r,
                   tau = 1/init_sigma/init_sigma)
    }else{
      print('Computing Bayes P1 solution with initial Values set to the SPECIFIED VALUES')
      inits = list(thetar = init_theta * pi/180,
                   r = init_r,
                   sigma = init_sigma)
    }
  }else if(init_type == 'random'){
    if(grepl('p3', mdlfile)){ # if regularized prior specify inits in terms of hip, knee and ankle angles.
      print('Computing Bayes P3 solution with initial Values set to the RANDOM VALUES')
      r_random = rt_ls(2, df = 5, mu = 0, a = 1)
      hip_random = rt_ls(1, df = 5, mu=-45, a = 15)
      knee_random = rt_ls(1, df=5, mu=30, a=6)
      ankle_random = rt_ls(1, df=5, mu=0, a=7)
      sigma_mm_random = rgamma(1, shape = 1.2, scale=0.1)
      inits = list(hip = hip_random,
                   knee = knee_random,
                   ankle = ankle_random,
                   r = r_random,
                   sigma_mm = sigma_mm_random)
    }else if(grepl('p2', mdlfile)){
      print('Computing Bayes P2 solution with initial Values set to the RANDOM VALUES')
      r_random = c(rnorm(1, true_vals[1], 0.001), dnorm(1, true_vals[2], 0.001))
      theta_random = runif(nlinks, -180,180)
      tau_random = rnorm(1, 1/(true_vals[nlinks+3]*true_vals[nlinks+3]), 1)
      inits = list(thetar = theta_random * pi/180,
                   r = r_random,
                   tau = tau_random)
    }else{
      print('Computing Bayes P1 solution with initial Values set to the RANDOM VALUES')
      theta_random = runif(nlinks, -180,180)
      r_random = runif(2, -0.1, 0.1)
      inits = list(thetar = theta_random * pi/180,
                   r = r_random,
                   sigma = runif(1, 0.1*true_vals[nlinks+3], 10*true_vals[nlinks+3]))
    }
    
  }else if(init_type == 'ls_est'){
    if(sum(is.na(ls_est))>0){stop("If using 'ls_est' as initials please supply LS estimates")}
    if(grepl('p3', mdlfile)){ # if regularized prior specify inits in terms of hip, knee and ankle angles.
      print('Computing Bayes P3 solution with initial Values set to the LS VALUES')
      hip_random = rt_ls(1, df = 5, mu=-45, a = 15)
      knee_random = rt_ls(1, df=5, mu=30, a=6)
      ankle_random = rt_ls(1, df=5, mu=0, a=7)
      
      inits = list(hip = hip_random,
                   knee = knee_random,
                   ankle = ankle_random,
                   r = c(ls_est[1], ls_est[2]),
                   sigma_mm = tail(ls_est,1)*1000)
    }else if(grepl('p2', mdlfile)){
      print('Computing Bayes P2 solution with initial Values set to the LS VALUES')
      inits = list(thetar = c(ls_est[3:(nlinks+2)]) * pi/180,
                   r = c(ls_est[1], ls_est[2]),
                   tau = 1/tail(ls_est,1)/tail(ls_est,1))
    }else{
      print('Computing Bayes P1 solution with initial Values set to the LS VALUES')
      inits = list(thetar = c(ls_est[3:(nlinks+2)]) * pi/180,
                   r = c(ls_est[1], ls_est[2]),
                   sigma = tail(ls_est,1))
    }
  }
  
  # Compile jags model
  t1 = Sys.time()
  m1 = jags.model(file = mdlfile,
                  data = jags_input,
                  inits = inits,
                  n.chains = 4,
                  n.adapt = 10000
  )
  
  # Sample from model
  out1 = coda.samples(model = m1,
                      variable.names = c("theta", "r", "sigma"),
                      n.iter = 100000,
                      thin = 5)
  t2 = Sys.time()
  
  
  return(list(fit = out1, time = as.numeric(t2-t1)))
}



