Files in this folder is for further analyzing homogeneous states observed near lower left corner of bifurcation diagrams.

## TangentContinuationLDS_bjorn.m:
Tangent Continuation with 4 preset parameters:\
chkst: step number between plots for monitoring\
stopt: a value in mu set between 0 and 1\
counter: count how many times solution passes the vertical line of stopt (set at 0)\
presetnum: stop continuation when counter > presetnum

## hom_bif_mu.m:
Analytic calcuation for mu values at bifurcation point when with homogeneous states and the corresponding eigenvalue for the Laplace term.\
Parameter: N, m, d (number of nodes, number of connected neighbors each side, coupling strength)\
Return: mu1, mu2, eval (2 mu values at bifurcation and eigenvalue of Laplace)

## testbif.m:
Use mu1 returned from hom_bif_mu.m as initial mu value\
Check if initial condition is a solution for LDS_RHS\
3 figures:\
&nbsp;&nbsp;&nbsp;&nbsp; 1. bifurcation plot constantly updated for monitoring\
&nbsp;&nbsp;&nbsp;&nbsp; 2. subplot 1: bifurcation diagram with intial point circled in red and end point circled in magenta; subplot 2: first 3 spectra\
&nbsp;&nbsp;&nbsp;&nbsp; 3. 3 line plots for differences in node values at each step (1. maximum of absolute difference in values between adjacent nodes; 2. value of 10th node minus u-plus(mu); 3. value of 10th node minus u-minus(mu))\
Save data file\
Parameters:\
N: number of nodes\
p0(1): coupling strength d\
steps: continuation steps\
stepsize: continuation stepsize\
adapt_stepsize: flag for Tangent Continuation\
m: number of connected neighbors to right and left (error message returned if number specified exceeds possible values)\
k: number of initial nodes at value a

## display_pattern.m:
Plot configuration of the eigenvector (manually designated) of the Laplace term calculated numerically from testbif.m
