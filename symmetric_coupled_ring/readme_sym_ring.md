LDS_Driver_neighbor.m:\
Symmetrically coupled ring with nearest m-neighbors right and left;\
Use Secant continuation;\
2 Figures: 1. Plot bifurcation, 2. Display 9 configurations across time\
\
certainplot.m:\
View configurations at certain steps (upon specification) in sol obtained from LDS_Driver_neighbor.m\
\
circulant.m:\
Return circulant matrix (file downloaded from MathWorks https://www.mathworks.com/matlabcentral/fileexchange/22858-circulant-matrix on 6/16/2020)\
\
parameters:\
N: number of nodes\
p0(1): coupling strength d\
p0(2): value of mu\
steps: continuation steps\
stepsize: continuation stepsize (default: 0.0005)\
m: number of nearest neighbors to right and left (error message returned if number specified exceeds possible values)\
u0: initial conditions (N*1 vector)
