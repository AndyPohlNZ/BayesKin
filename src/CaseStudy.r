############################################################

# ---------------------- Preliminaries -----------------------------------------
# preliminaries  # package installs
rm(list = ls()) # clear workspace
library('dplyr')
library('tidyr')
library('latex2exp')
library('xtable')
library('coda')
#library('rjags')
#load.module("vonmises")

SRC_DIR = "" # location of results directory
setwd(SRC_DIR)
source("library.r")

# Color plaette for plotting
colors = list(grey = rgb(97/255, 97/255,97/255, 0.2),
              orange = rgb(241/255, 136/255,5/255, 1),
              blue = rgb(18/255, 78/255,120/255,  1),
              red = rgb(162/255, 0/255,33/255, 1),
              green = rgb(76/255, 159/255, 112/255, 1),
              purple = rgb(217/255, 187/255, 249/255, 1),
              pink = rgb(217/255, 3/255, 104/255, 1))

# Formatted plot labels
PARMLABS = list(r1=TeX("$r_1 \\,  (m)$"), r2=TeX("$r_2  \\, (m)$"), theta1=TeX("$\\theta_1 \\, (^o)$"), 
                theta2=TeX("$\\theta_2 \\, (^o)$"), theta3=TeX("$\\theta_3\\, (^o) "), sigma=TeX("$\\sigma \\, (m)$"))
# sampling frequency
dt = 0.005

# Basic functions
gen_rotation_matrix = function(theta){
    # Generates a rotation matrix for angle theta specified in radians
    return(rbind(c(cos(theta), -sin(theta)), c(sin(theta), cos(theta))))
}
gen_vec_angle = function(v){
    # Generates the angle a vector v makes with the x axis
    angle = abs(atan2(v[2], v[1]))
    return(angle)
}

calc_knee_joint_angle = function(theta1, theta2){
    # Calculates the knee angle as per Buzeck 1990
    if (abs(theta2) < abs(theta1)){
        return(abs(theta2) - abs(theta1))
    }else{
        return(-(abs(theta1) - abs(theta2)))
    }
}

calc_ankle_joint_angle = function(theta2, theta3){
    # Calculates the ankle angle as per Buzeck 1990
    
    if ((theta2 < 0) & (theta3<0)){
        return(abs(theta3) - abs(theta2) + 90)
    }else if ((theta2<0) & (theta3 > 0)){
        return(180 - abs(theta2) - abs(theta3) - 90)
    }
}


# ---------------------- Data Loading/Processing --------------------------------
# Z-Y = sagittal plane...
# load data for random subject.
demographics = read.table("../data/s2_demographics.csv", sep=",", header=TRUE)
static = read.table("../data/s2_staticTrial.csv", sep=",", header = TRUE)
events = read.table("../data/s2_events.csv", sep=",", header = TRUE)

# Re-align coordinate system
static$Z = -static$Z # reflect
static[,c("X", "Y", "Z")] = static[,c("X", "Y", "Z")]/1000 # change to m

# dynamic trial
dynamic = read.table("../data/s2_dynamicTrial.csv", sep = ",", header = TRUE)
dynamic$Z = -dynamic$Z # reflect
dynamic[,c("X", "Y", "Z")] = dynamic[,c("X", "Y", "Z")]/1000 #change to m

# Comptue Joint Centers
HJC = static[static$marker=='RHip', c("Z", "Y")]
KJC = static[static$marker=='LatKnee', c("Z", "Y")] + (static[static$marker=='MedKnee', c("Z", "Y")] - static[static$marker=='LatKnee', c("Z", "Y")])/2
AJC = static[static$marker=='LatAnkle', c("Z", "Y")] + (static[static$marker=='MedAnkle', c("Z", "Y")] - static[static$marker=='LatAnkle', c("Z", "Y")])/2
TOE = static[static$marker=='Toe', c("Z", "Y")]

# Compute LCSs
thigh = rbind(static[static$marker == "Thigh1", c('Z', "Y")],
              static[static$marker == "Thigh2", c('Z', "Y")],
              static[static$marker == "Thigh3", c('Z', "Y")],
              static[static$marker == "Thigh4", c('Z', "Y")])


shank = rbind(static[static$marker == "Shank1", c('Z', "Y")],
              static[static$marker == "Shank2", c('Z', "Y")],
              static[static$marker == "Shank3", c('Z', "Y")],
              static[static$marker == "Shank4", c('Z', "Y")])


foot = rbind(static[static$marker == "Foot1", c('Z', "Y")],
             static[static$marker == "Foot2", c('Z', "Y")],
             static[static$marker == "Foot3", c('Z', "Y")],
             static[static$marker == "Foot4", c('Z', "Y")])


# Plot of static trial.
plot(rbind(HJC, KJC), type = 'o', pch = 16, col = rgb(18/255, 78/255,120/255,1),
     xlim = c(-0.6, 0.0), ylim = c(0,1), asp = 1,
     xlab = 'x (m)', ylab = 'y (m)',
     bty='n')
points(thigh, pch = 4,col = rgb(18/255, 78/255,120/255,1))
lines(rbind(KJC, AJC), type = 'o', pch = 16, col = rgb(76/255, 159/255, 112/255, 1))
points(shank, pch = 4, col = rgb(76/255, 159/255, 112/255, 1))
lines(rbind(AJC, TOE), type = 'o', pch = 16, col = rgb(217/255, 3/255, 104/255, 1))
points(foot, pch = 4, col = rgb(217/255, 3/255, 104/255, 1))


# Realign to the horizontal static trial expected by Bayesian Model
thigh_lcs_origin = c(0,0)
thigh_lcs_termination = as.numeric(KJC - HJC)
angle = gen_vec_angle(thigh_lcs_termination)
rot.mat = gen_rotation_matrix(angle)
thigh_lcs_termination=  rot.mat %*% thigh_lcs_termination
thigh_lcs = do.call(rbind, apply(thigh, 1, function(x) return(x - HJC)))
thigh_lcs = apply(thigh_lcs, 1, function(x) return(rot.mat %*% x))

plot(rbind(thigh_lcs_origin, t(thigh_lcs_termination)), type = 'o', pch = 16, 
     xlim = c(-1, 1), ylim = c(-.8,1), asp = 1,
     main = "Thigh LCS", xlab = 'x', ylab = 'y')
points(t(thigh_lcs), pch = 4)

shank_lcs_origin = c(0,0)
shank_lcs_termination = as.numeric(AJC - KJC)
angle = gen_vec_angle(shank_lcs_termination)
rot.mat = gen_rotation_matrix(angle)
shank_lcs_termination=  rot.mat %*% shank_lcs_termination
shank_lcs = do.call(rbind, apply(shank, 1, function(x) return(x - KJC)))
shank_lcs = apply(shank_lcs, 1, function(x) return(rot.mat %*% x))

plot(rbind(shank_lcs_origin, t(shank_lcs_termination)), type = 'o', pch = 16, 
     xlim = c(-1, 1), ylim = c(-.8,1), asp = 1,
     main = "Shank LCS", xlab = 'x', ylab = 'y')
points(t(shank_lcs), pch = 4)

foot_lcs_origin = c(0,0)
foot_lcs_termination  = as.numeric(TOE - AJC)
angle = gen_vec_angle(foot_lcs_termination)
rot.mat = gen_rotation_matrix(angle)
foot_lcs_termination=  rot.mat %*% foot_lcs_termination
foot_lcs = do.call(rbind, apply(foot, 1, function(x) return(x - AJC)))
foot_lcs = apply(foot_lcs, 1, function(x) return(rot.mat %*% x))
plot(rbind(foot_lcs_origin, t(foot_lcs_termination)), type = 'o', pch = 16, 
     xlim = c(-1, 1), ylim = c(-.8,1), asp = 1,
     main = "Foot LCS", xlab = 'x', ylab = 'y')
points(t(foot_lcs), pch = 4)

# Add row/colnames
rownames(thigh_lcs) = c("X", "Y")
colnames(thigh_lcs) = c("thigh1", "thigh2", "thigh3", "thigh4")

rownames(shank_lcs) = c("X", "Y")
colnames(shank_lcs) = c("shank1", "shank2", "shank3", "shank4")

rownames(foot_lcs) = c("X", "Y")
colnames(foot_lcs) = c("foot1", "foot2", "foot3", "foot4")

# Compute Length of each segment
thigh_L = sqrt(sum(thigh_lcs_termination - thigh_lcs_origin)^2)
shank_L = sqrt(sum(shank_lcs_termination - shank_lcs_origin)^2)
foot_L = sqrt(sum(foot_lcs_termination - foot_lcs_origin)^2)


# Compute in format expected by bayesian model
link1 = list(length = thigh_L,
             x = thigh_lcs)
link2 = list(length = shank_L,
             x = shank_lcs)
link3 = list(length = foot_L,
             x = foot_lcs)

# links = list(thigh=link1, shank=link2, foot =link3)
links = list(thigh = link1, shank = link2, foot = link3)

# Get nstrides strides from data
nstrides = 20
events = events[complete.cases(events),]
events = events[10:(10+nstrides-1),]
events = events * dt # convert events to time

# Extract stance phases
stancephases = list()
for(i in 1:nrow(events)){
    stancei = dynamic[(dynamic$time > events[i,1]) & (dynamic$time < events[i, 2]),]
    stancephases[[i]] = stancei
}


# Normalise to 101 data points
stancephasenorm = list()
for(i in 1:length(stancephases)){
    stancephase = stancephases[[i]]
    markers = unique(stancephase$marker)
    stancephasenorm[[i]] = list()
    for(j in 1:length(markers)){
        markerj = markers[j]
        stancephasemarker = stancephase[stancephase$marker == markerj,]
        interpX = approx(stancephasemarker$time, stancephasemarker$Z, n = 101)
        interpY = approx(stancephasemarker$time, stancephasemarker$Y, n = 101)
        toadd = cbind(interpX$x, interpX$y, interpY$y)
        colnames(toadd) = c('time', "X","Y")
        stancephasenorm[[i]][[j]] =toadd
    }
    names(stancephasenorm[[i]]) = markers
    
}

# Compute x in the format expected for bayesian model for each frame
maxtime = nrow(stancephasenorm[[1]][[1]])
nstrides = length(stancephasenorm)
nlinks = length(links)
nmarkers = sapply(links, function(x){ncol(x$x)})
M = 2
seg_lengths = sapply(links, function(x){x$length})
x = lapply(links, function(x){x$x})
xjags  = array(NA, dim = c(nlinks, max(nmarkers), 2))
for(i in 1:nlinks){
    xjags[i,,] = t(x[[i]])
}

# links = list(thigh=link1, shank=link2, foot =link3)
links = list(thigh = link1, shank = link2, foot = link3)

################################################################################
# Run analysis across time 
################################################################################
maxtime = 5
dt = 0.005
bayes_results = matrix(nrow = maxtime/dt, ncol = 6*3)
bayes_results_knee = matrix(nrow = maxtime/dt, ncol =3)
bayes_results_ankle = matrix(nrow = maxtime/dt, ncol =3)

ls_results = matrix(nrow = maxtime/dt, ncol = 6*3)
ls_results_knee = matrix(nrow=maxtime/dt, 1)
ls_results_ankle = matrix(nrow=maxtime/dt, 1)

timei = 0
idx = 1
while(timei < maxtime){
    print(paste("Computing for solution for time:", timei, "s"))
    dynamic_t = dynamic %>% filter(time ==timei)
    
    y1 = t(rbind(dynamic_t[dynamic_t$marker == "Thigh1", c('Z', "Y")],
                 dynamic_t[dynamic_t$marker == "Thigh2", c('Z', "Y")],
                 dynamic_t[dynamic_t$marker == "Thigh3", c('Z', "Y")],
                 dynamic_t[dynamic_t$marker == "Thigh4", c('Z', "Y")]))
    colnames(y1) = c("thigh1","thigh2","thigh3","thigh4")
    rownames(y1) = c("X","Y")
    
    y2 = t(rbind(dynamic_t[dynamic_t$marker == "Shank1", c('Z', "Y")],
                 dynamic_t[dynamic_t$marker == "Shank2", c('Z', "Y")],
                 dynamic_t[dynamic_t$marker == "Shank3", c('Z', "Y")],
                 dynamic_t[dynamic_t$marker == "Shank4", c('Z', "Y")]))
    colnames(y2) = c("shank1","shank2","shank3","shank4")
    rownames(y2) = c("X","Y")
    
    y3 = t(rbind(dynamic_t[dynamic_t$marker == "Foot1", c('Z', "Y")],
                 dynamic_t[dynamic_t$marker == "Foot2", c('Z', "Y")],
                 dynamic_t[dynamic_t$marker == "Foot3", c('Z', "Y")],
                 dynamic_t[dynamic_t$marker == "Foot4", c('Z', "Y")]))
    colnames(y3) = c("foot1","foot2","foot3","foot4")
    rownames(y3) = c("X","Y")
    
    y = list(y1, y2, y3)
    
    LS_result = LS_soln(y, links, 
                        inits = c(0,0, -90, -90, -90), 
                        init_type = "specified")
    ls_results[idx, ] = c(LS_result$r.hat[1], LS_result$intervals[1,],
                          LS_result$r.hat[2], LS_result$intervals[2,],
                          LS_result$theta.hat[1],LS_result$intervals[3,],
                          LS_result$theta.hat[2], LS_result$intervals[4,],
                          LS_result$theta.hat[3], LS_result$intervals[5,],
                          LS_result$sigma.hat[1], LS_result$intervals[6,])
    ls_results_knee[idx] = calc_knee_joint_angle(LS_result$theta.hat[1], LS_result$theta.hat[2])
    ls_results_ankle[idx] =calc_ankle_joint_angle(LS_result$theta.hat[2], LS_result$theta.hat[3])
    
    
    Bayes_result = Bayes_soln(y=y,
                              links = links,
                              mdlfile = './BayesModel_p3.jags',
                              init_type = 'specified',
                              true_vals = c(0,0, -90, -90,-90, 2e-3),
                              ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))
    sampsmat = as.matrix(Bayes_result$fit)

    Bayes_hat = list(r.hat = summary(Bayes_result$fit)$statistics[c('r[1]', 'r[2]'), c("Mean")],
                     theta.hat = summary(Bayes_result$fit)$statistics[c('theta[1]', 'theta[2]', 'theta[3]'), c("Mean")],
                     sigma.hat =summary(Bayes_result$fit)$statistics[c('sigma'), c("Mean")])

    temp = summary(Bayes_result$fit)
    Bayes_ci = list(r.hat = temp$quantiles[c('r[1]', 'r[2]'), c('2.5%', '97.5%')],
                    theta.hat = temp$quantiles[c('theta[1]', 'theta[2]', 'theta[3]'), c('2.5%', '97.5%')],
                    sigma.hat = temp$quantiles[c('sigma'), c('2.5%', '97.5%')])

    
    bayes_results[idx,] = c(Bayes_hat$r.hat[1], Bayes_ci$r.hat[1,],
                            Bayes_hat$r.hat[2], Bayes_ci$r.hat[2,],
                            Bayes_hat$theta.hat[1], Bayes_ci$theta.hat[1,],
                            Bayes_hat$theta.hat[2], Bayes_ci$theta.hat[2,],
                            Bayes_hat$theta.hat[3], Bayes_ci$theta.hat[3,],
                            Bayes_hat$sigma.hat[1], Bayes_ci$sigma.hat)
    
    bayes_knee = apply(cbind(sampsmat[,'theta[1]'], sampsmat[,'theta[2]']), 1, function(x) calc_knee_joint_angle(x[1], x[2]))
    bayes_results_knee[idx, ] = c(mean(bayes_knee), quantile(bayes_knee, probs = c(0.025, 0.975)))
    
    bayes_ankle =apply(cbind(sampsmat[,'theta[2]'], sampsmat[,'theta[3]']), 1, function(x) calc_ankle_joint_angle(x[1], x[2]))
    bayes_results_ankle[idx, ] = c(mean(bayes_ankle), quantile(bayes_ankle, probs = c(0.025, 0.975)))
    
    posture_hat_ls = gen_posture(links = links, r = LS_result$r.hat, theta = LS_result$theta.hat)
    
    posture_hat_bayes = gen_posture(links = links, r = Bayes_hat$r.hat, theta = Bayes_hat$theta.hat)
    
    # Plots an 'animation of the case study data
    # png(paste0("./caseStudy/", timei*1000, ".png"), width=673, height = 528)
    # plot(NA, xlim = c(-0.5, 1), ylim = c(0, 1),
    #      main = paste("Time: ", timei, ('s')),
    #      xlab="X (m)", ylab = "Y (m)",
    #      bty = "n")
    # points(t(y1), pch = 3 )
    # points(t(y2), pch = 3 )
    # points(t(y3), pch = 3 )
    #
    # points(t(do.call(cbind,posture_hat_ls$J)), type = "o", pch = 20, col = colors$pink) # mu_tmp
    # points(t(do.call(cbind,posture_hat_bayes$J)), type = "o", pch = 20, col = colors$blue) # mu_tmp
    # points(t(do.call(cbind,posture_hat_ls$alpha)), pch = 20, col = colors$pink)
    # points(t(do.call(cbind,posture_hat_bayes$alpha)), pch = 20, col =colors$blue)
    # dev.off()
    timei = round(timei + dt, 3)
    idx = idx +1
}


################################################################################
# Analysis plots
################################################################################
# Parameters
timeseq = seq(0, 4.995, by = dt)
events = events[1:8,]
events[1,1] = 0

par(mfrow=c(5,1), mar=c(2, 4, 0, 2))
# r1
plot(NULL, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(0,5), ylim = c(0.1, 0.25),
    xlab = NA, ylab = PARMLABS[['r1']],
     bty = 'n', xaxt='n')

for(i in 1:nrow(events)){
    rect(xleft = events[i,1], ybottom = 0,  
         xright = events[i,2], ytop = 90, lty=0,
         col = rgb(237/255, 237/255, 237/255,1)) 
    if(i >1){
    abline(v = c(events[i,1], events[i,2]), col=rgb(241/255, 136/255,5/255, 1))
    }else{
        abline(v =  events[i,2], col=rgb(241/255, 136/255,5/255, 1))
        
    }
}
lines(timeseq, bayes_results[,1], type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,2], type = 'l', lty =2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,3], type = 'l', lty = 2, col = rgb(18/255, 78/255,120/255, 1))

lines(timeseq, ls_results[,1], type = 'l', lty =1, col = rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,2], type = 'l', lty =2, col =  rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,3], type = 'l', lty = 2, col =  rgb(217/255, 3/255, 104/255, 1))

#r2
plot(NULL, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(0,5), ylim = c(0.75, 0.90),
     xlab = NA, ylab = PARMLABS[['r2']],
     bty = 'n', xaxt='n')

for(i in 1:nrow(events)){
    rect(xleft = events[i,1], ybottom = 0,  
         xright = events[i,2], ytop = 90, lty=0,
         col = rgb(237/255, 237/255, 237/255,1)) 
    if(i >1){
        abline(v = c(events[i,1], events[i,2]), col=rgb(241/255, 136/255,5/255, 1))
    }else{
        abline(v =  events[i,2], col=rgb(241/255, 136/255,5/255, 1))
        
    }
}
     
lines(timeseq, bayes_results[,4], type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,5], type = 'l', lty =2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,6], type = 'l', lty = 2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, ls_results[,4], type = 'l', lty =1, col = rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,5], type = 'l', lty =2, col =  rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,6], type = 'l', lty = 2, col =  rgb(217/255, 3/255, 104/255, 1))

#theta1
plot(NULL, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(0,5), ylim = c(-120, -60),
     xlab = NA, ylab = PARMLABS[['theta1']],
     bty = 'n', xaxt='n')

for(i in 1:nrow(events)){
    rect(xleft = events[i,1], ybottom = -180,  
         xright = events[i,2], ytop = 180, lty=0,
         col = rgb(237/255, 237/255, 237/255,1)) 
    if(i >1){
        abline(v = c(events[i,1], events[i,2]), col=rgb(241/255, 136/255,5/255, 1))
    }else{
        abline(v =  events[i,2], col=rgb(241/255, 136/255,5/255, 1))
        
    }
}
lines(timeseq, bayes_results[,7], type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,8], type = 'l', lty =2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,9], type = 'l', lty = 2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, ls_results[,7], type = 'l', lty =1, col = rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,8], type = 'l', lty =2, col =  rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,9], type = 'l', lty = 2, col =  rgb(217/255, 3/255, 104/255, 1))

#theta2
plot(NULL, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(0,5), ylim = c(-175, -75),
     xlab = NA,, ylab = PARMLABS[['theta2']],
     bty = 'n', xaxt='n')
for(i in 1:nrow(events)){
    rect(xleft = events[i,1], ybottom = -180,  
         xright = events[i,2], ytop = 180, lty=0,
         col = rgb(237/255, 237/255, 237/255,1)) 
    if(i >1){
        abline(v = c(events[i,1], events[i,2]), col=rgb(241/255, 136/255,5/255, 1))
    }else{
        abline(v =  events[i,2], col=rgb(241/255, 136/255,5/255, 1))
        
    }
}
lines(timeseq, bayes_results[,10], type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,11], type = 'l', lty =2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,12], type = 'l', lty = 2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, ls_results[,10], type = 'l', lty =1, col = rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,11], type = 'l', lty =2, col =  rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,12], type = 'l', lty = 2, col =  rgb(217/255, 3/255, 104/255,))

#theta3
plot(NULL, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(0,5), ylim = c(-100, 0),
     xlab = NA, ylab = PARMLABS[['theta3']],
     bty = 'n')
mtext("Time (s)", side=1, line=1)

for(i in 1:nrow(events)){
    rect(xleft = events[i,1], ybottom = -180,  
         xright = events[i,2], ytop = 180, lty=0,
         col = rgb(237/255, 237/255, 237/255,1)) 
    if(i >1){
        abline(v = c(events[i,1], events[i,2]), col=rgb(241/255, 136/255,5/255, 1))
    }else{
        abline(v =  events[i,2], col=rgb(241/255, 136/255,5/255, 1))
        
    }
}
lines(timeseq, bayes_results[,13], type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,14], type = 'l', lty =2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,15], type = 'l', lty = 2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, ls_results[,13], type = 'l', lty =1, col = rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,14], type = 'l', lty =2, col =  rgb(217/255, 3/255, 104/255, 1))
lines(timeseq, ls_results[,15], type = 'l', lty = 2, col =  rgb(217/255, 3/255, 104/255,))


################################################################################
# Plot of Sigma
################################################################################
par(mfrow=c(1,1))
plot(NULL, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(0,5), ylim = c(1, 10),
     xlab = NA, ylab = TeX("$\\sigma \\, (mm)$"),
     bty = 'n')
mtext("Time (s)", side=1, line=1)
for(i in 1:nrow(events)){
    rect(xleft = events[i,1], ybottom = -180,  
         xright = events[i,2], ytop = 180, lty=0,
         col = rgb(237/255, 237/255, 237/255,1)) 
    if(i >1){
        abline(v = c(events[i,1], events[i,2]), col=rgb(241/255, 136/255,5/255, 1))
    }else{
        abline(v =  events[i,2], col=rgb(241/255, 136/255,5/255, 1))
        
    }
}
lines(timeseq, bayes_results[,16]*1000, type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,17]*1000, type = 'l', lty =2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results[,18]*1000, type = 'l', lty = 2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, ls_results[,16]*1000, type = 'l', lty =1, col = rgb(217/255, 3/255, 104/255, 1))

######################################################################################
## Plot knee/ankle joint angles
######################################################################################

kneeplot = bayes_results_knee[((timeseq>=events$RFS[4]) & (timeseq<=events$RFS[5])), ]
kneeplottime = timeseq[((timeseq>=events$RFS[4]) & (timeseq<=events$RFS[5]))]
kneeplotevents = events[4,]
plot(NA, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(min(kneeplottime),max(kneeplottime)), ylim = c(0, 90),
     xlab = "Time (s)", ylab = TeX("$\\theta_{knee} \\, (^o)$"),
     bty = 'n')
rect(xleft = kneeplotevents[,1], ybottom = 0,  
     xright = kneeplotevents[,2], ytop = 90, lty=0,
     col = rgb(237/255, 237/255, 237/255,1))
abline(v = c(kneeplotevents[,1], kneeplotevents[,2]), col=rgb(241/255, 136/255,5/255, 1))
lines(timeseq, bayes_results_knee[,1], type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results_knee[,2], type = 'l', lty =2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results_knee[,3], type = 'l', lty = 2, col = rgb(18/255, 78/255,120/255, 1))


## Ankle
ankleplot = bayes_results_ankle[((timeseq>=events$RFS[4]) & (timeseq<=events$RFS[5])), ]
ankleplottime = timeseq[((timeseq>=events$RFS[4]) & (timeseq<=events$RFS[5]))]
ankleplotevents = events[4,]
plot(NA, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(min(kneeplottime),max(kneeplottime)), ylim = c(-20, 50),
     xlab = "Time (s)", ylab = TeX("$\\theta_{ankle} \\, (^o)$"),
     bty = 'n')
rect(xleft = ankleplotevents[,1], ybottom = -20,  
     xright = ankleplotevents[,2], ytop = 90, lty=0,
     col = rgb(237/255, 237/255, 237/255,1))
abline(v = c(ankleplotevents[,1], ankleplotevents[,2]), col=rgb(241/255, 136/255,5/255, 1))
lines(timeseq, bayes_results_ankle[,1], type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results_ankle[,2], type = 'l', lty =2, col = rgb(18/255, 78/255,120/255, 1))
lines(timeseq, bayes_results_ankle[,3], type = 'l', lty = 2, col = rgb(18/255, 78/255,120/255, 1))
abline(v=2.215)

dynamic_t = dynamic %>% filter(time ==2.215) # 2.175 = max during stance, 2.545 = max during swing. # 2.61 = max ankle 2.215 = max ankle stance

y1 = t(rbind(dynamic_t[dynamic_t$marker == "Thigh1", c('Z', "Y")],
             dynamic_t[dynamic_t$marker == "Thigh2", c('Z', "Y")],
             dynamic_t[dynamic_t$marker == "Thigh3", c('Z', "Y")],
             dynamic_t[dynamic_t$marker == "Thigh4", c('Z', "Y")]))
colnames(y1) = c("thigh1","thigh2","thigh3","thigh4")
rownames(y1) = c("X","Y")

y2 = t(rbind(dynamic_t[dynamic_t$marker == "Shank1", c('Z', "Y")],
             dynamic_t[dynamic_t$marker == "Shank2", c('Z', "Y")],
             dynamic_t[dynamic_t$marker == "Shank3", c('Z', "Y")],
             dynamic_t[dynamic_t$marker == "Shank4", c('Z', "Y")]))
colnames(y2) = c("shank1","shank2","shank3","shank4")
rownames(y2) = c("X","Y")

y3 = t(rbind(dynamic_t[dynamic_t$marker == "Foot1", c('Z', "Y")],
             dynamic_t[dynamic_t$marker == "Foot2", c('Z', "Y")],
             dynamic_t[dynamic_t$marker == "Foot3", c('Z', "Y")],
             dynamic_t[dynamic_t$marker == "Foot4", c('Z', "Y")]))
colnames(y3) = c("foot1","foot2","foot3","foot4")
rownames(y3) = c("X","Y")

y = list(y1, y2, y3)
Bayes_result = Bayes_soln(y=y,
                          links = links,
                          mdlfile = './BayesModel_p3.jags',
                          init_type = 'specified',
                          true_vals = c(0,0, -90, -90,-90, 2e-3),
                          ls_est = c(LS_result$r.hat, LS_result$theta.hat, LS_result$sigma.hat))
sampsmat = as.matrix(Bayes_result$fit)

randidxs = sample(1:nrow(sampsmat))[1:200]

plot(NULL, NULL,
     xlim = c(-1, 1), ylim = c(0, 1),
     xlab = NA, ylab = NA, bty='n', xaxt='n', yaxt='n')
for(i in 1:200){
    randidxi = randidxs[i]
    posture_hat_bayes = gen_posture(links = links, r = sampsmat[randidxi, c(1,2)], theta = sampsmat[randidxi, c(4,5,6)])
    
    points(t(do.call(cbind,posture_hat_bayes$J)), type = "o", pch = 20, col = rgb(18/255, 78/255,120/255, 0.1)) # mu_tmp
}


######################################################################################
## Plot credible interval width
######################################################################################
plot(NA, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(min(kneeplottime),max(kneeplottime)), ylim = c(0, 20),
     xb = "Time (s)", ylab = TeX("95% credible interval width $\\theta_{knee}$ $(^o)$"),
     bty = 'n')
rect(xleft = ankleplotevents[,1], ybottom = -20,  
     xright = ankleplotevents[,2], ytop = 90, lty=0,
     col = rgb(237/255, 237/255, 237/255,1))
abline(v = c(ankleplotevents[,1], ankleplotevents[,2]), col=rgb(241/255, 136/255,5/255, 1))
lines(timeseq, bayes_results_knee[,3]-bayes_results_knee[,2], type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))


plot(NA, type = 'l', col = rgb(18/255, 78/255,120/255, 1),
     xlim = c(min(kneeplottime),max(kneeplottime)), ylim = c(0, 20),
     xlab = "Time (s)", ylab = TeX("95% credible interval width $\\theta_{ankle}$ $(^o)$"),
     bty = 'n')
rect(xleft = ankleplotevents[,1], ybottom = -20,  
     xright = ankleplotevents[,2], ytop = 90, lty=0,
     col = rgb(237/255, 237/255, 237/255,1))
abline(v = c(ankleplotevents[,1], ankleplotevents[,2]), col=rgb(241/255, 136/255,5/255, 1))
lines(timeseq, bayes_results_ankle[,3]-bayes_results_ankle[,2], type = 'l', lty =1, col = rgb(18/255, 78/255,120/255, 1))
