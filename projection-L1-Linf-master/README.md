# Projection onto the intersection of the L1-ball and L-infinity-ball

This repository contains a Matlab/Mex/C++ implementation of a linear time algorithm for the projection onto the intersection of the L1-ball and L-infinity-ball (box constraint). 
The derivation and description of the algorithm can be found [here](https://github.com/anastasia-podosinnikova/projection-L1-Linf/files/13646/script.pdf).

If you use this code or the derivation of the algorithm for your research, please, cite:
*A. Podosinnikova. Robust Principal Component Analysis as a Nonlinear Eigenproblem. Master's Thesis, Saarland University, Department of Mathematics and Computer Science, 2013.*


## Using this code with Matlab

- make sure your Matlab recognizes your gcc compiler
```
mex -setup C++
```
- compile the ```proj_L1_box.cpp``` function
```
install.m
```
- see an example of the code's usage
```
example_projection.m
```


## Questions?
Do not hesitate to contact me: firstname.lastname@inria.fr (Anastasia Podosinnikova)



