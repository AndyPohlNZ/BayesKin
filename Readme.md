# Bayesian Kinematics

<<<<<<< HEAD
This repository contains code and data to replicate analyses outlined in Pohl, Schofield and Ferber. (2019). Compoaring the performance of Bayesian and least-squares approaches to inverse kinematics.
=======
This repository contains code and data to replicate analyses outlined in Pohl, Schofield and Ferber. (2021). Examination of a Bayesian approach to inverse Kinematics.
>>>>>>> master

If you have any queries please contact the author at andrew.pohl@ucalgary.ca

## Requirements
The scripts contained in this repository require R version 3.6 or higher and the JAGS gibbs sampler.  These are freely avaiable and we refer users to the relevant documentation for instalation at:

R: https://www.r-project.org/
JAGS: https://sourceforge.net/projects/mcmc-jags/

In addition users may wish to install the R studio IDE which is freely avaliable from: https://www.rstudio.com/

The models and analysis code presented rely on the following user contributed R pacakges:
- coda: https://cran.r-project.org/web/packages/coda/index.html
- numDeriv: https://cran.r-project.org/web/packages/numDeriv/index.html
- rjags: https://cran.r-project.org/web/packages/rjags/index.html

In addition the following packages are used for visulisation:
- dplyr: https://cran.r-project.org/web/packages/dplyr/index.html
- latex2exp: https://cran.r-project.org/web/packages/latex2exp/index.html
- xtable: https://cran.r-project.org/web/packages/xtable/index.html


Packages can be installed from the Comprehensive R archive Network (CRAN) by running the following in an R terminal:

    # Required Packages
    install.packages('coda')
    install.packages('numDeriv')
    install.packages('rjags')
    # Optional packages for visulaisation
    install.packages('dplyr')
    install.packages('latex2exp')
    install.packages('xtable')


## Getting Started
A  `example.r` script has been included to guide readers through successfully simulating an observation of a 3 link kinematic chain and fitting the Bayesian or LS models described in the article.

In each script a user needs to specify the directory where this repository has been extracted so that the necessary `library.r' and JAGS models can be accessed.  This is performed by completing the line:
    WORKING_DIR = "<LOCATION OF THE BAYES KIN REPOSITORY>" 

`example.r` first simulates an observation of a 3 link kinematic pose with pose parameters given by `seg.length`, `r_true`, `theta_true` and `sigma_true`.  Then a

## Directory contents
In the `src` directory the following can be found:

- `library.r`: Contains helper functions for running simulations and fitting models.
- `RunMdl_<Single/Double/Triple>Link.r`: simulates and fit LS and Bayesian models to the 1000 single, double or triple link kinematic chains described in the manuscript. Results are saved in a custom model output format wtihin the relevent `./<SingleLink/DoubleLink/TripleLink>_Mdl/` directories.
- `process_results.r`: Processes model output for subsequent analysis
- `analyse_results.r`: Produces figures and tables found within the manuscript
- `BayesModel_<p1/p2/.../p5>.jags`: Contain JAGS code for fitting the Bayesian models outlined in the manuscript.
- The contents of `/AnalyticalCostGrad/` and `/AnalyticalHessian/` contain analytical expressions for the Gradient and Hessian matrix of the cost function used for LS optimisation in text format.  This improves accuracy and speed over numerical approximations to gradients and the Hessian requried for optimisation of the LS cost.
- `/SupplementD/`: Contains code and models for replicating the results provided in Supplement Material D.

The `SupplementaryMaterial` directory contains pdf's and raw latex files for the Supplementary Material referenced in the original manuscript.

## Replication of Results
To replicate the results provided in the main manuscript the scripts should be run in order either from within RStudio or via the command line:
1. `RunMdl_SingleLink.r`
2. `RunMdl_DoubleLink.r`
3. `RunMdl_TripleLink.r`
4. `process_results.r`
5. `analyse_results.r`

