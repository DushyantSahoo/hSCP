# hSCP

# Overview

MATLAB code for estimating Hierarchical Sparse Connectivity Patterns (hSCP) using fMRI data.

contact: sahoodushyant@gmail.com

# Main function

The current implementation requires a symmetric matrix, for example, a connectivity matrix as an input. Though the user can give any type of square matrix as an input, the analysis is done using a symmetric decomposition. Below is the main function, and description about input and output. 

```[W, lambda, error] = hscp_amsgrad(A,k,alpha,loop,eta,beta1,beta2,eps,tole,svd_check)```

Below are the inputs to the above function

1) A contains the input matrices in in a cell format, for example- A{1}, A{2}, ....

2) k contains the number of components in each hierarchy, for example k = [120,20,4] where 120 is the number of nodes in data or the size of input matrix
 20 is the number of components at level 1,
 4 is the number of components at level 2.
 k can also be [120,10] having only one level or [120,40,20,4] having
 three levels

3) alpha is the sparsity level of each components at each level, you will
 have to play with it a little bit to get a balance between sprasity in
 the components and good approximation error.
 If k = [120,20,4] then alpha can be [1,1] i.e. sparsity for each level

 4) loop is the number of iterations of gradient descent

 5) eta, beta1, beta2, eps are the hyperparameters for amsgrad

 6) tole is the % change in the error before gradient descent stops

 7) svd_check is used for initializing the algorithm with SVD, details are
 given in the main paper
 
 Below are the typical hyperparameter settings that would work-
 1) svd_check = 1
 2) loop = 6000
 3) eta = 0.1
 4) beta1 = 0.99
 5) beta2 = 0.999
 6) eps = 10^-8
 7) tol1 = 10^(-4)

 There are three outputs-
 1) W stores components at different level in cells, each cell of W will
 store components at each level
 2) lambda stores subject specific information in cell, each cell of
 lambda will store subject and level specific information
 3) error stores % information captured by the decomposition, ideally it
 should be decreasing with the iterations.

MATLAB for simulating data and running the code is also given which would give user an idea of the input parameters and how the output looks. Please refer to "Hierarchical extraction of functional connectivity components in human brain using resting-state fMRI" paper for more details. 

# Reference

D. Sahoo, T. D. Satterthwaite and C. Davatzikos, "Hierarchical extraction of functional connectivity components in human brain using resting-state fMRI," in IEEE Transactions on Medical Imaging, doi: 10.1109/TMI.2020.3042873.

Podosinnikova, Anastasia. Robust Principal Component Analysis as a Nonlinear Eigenproblem. Diss. Universität des Saarlandes Saarbrücken, 2013.

Reddi, Sashank J., Satyen Kale, and Sanjiv Kumar. "On the convergence of adam and beyond." arXiv preprint arXiv:1904.09237 (2019).
