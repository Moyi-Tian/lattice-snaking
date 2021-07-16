Continuation methods:\
SecantContinuationLDS.m:\
Munually preset chkst (number of steps skipped between plots shown for monitoring purpose)\
Return eflag (steps where fsolve gives exit flag with value less than or equal to 0)\
TangentContinuationLDS.m:\
Munually preset chkst (number of steps skipped between plots shown for monitoring purpose)\
Adjust step size if pre-given flag == 1\
\
\
test.m:\
Use Secant Continuation\
Prepare data, generate movie of bifurcation plots and save movie\
\
testbif.m:\
Mually select between Secant Continuation and Tangent Continuation (if latter then also set adapt_stepsize flag)\
Prepare data\
Figure 1: for monitoring during continuation process\
Figure 2: subplot 1 - bifurcation plot (initial circled red, end circled magenta); subplot 2 - plot first 3 spectra (spec only returned if using Tangent Continuation)\
Print absolute step change in mu, norm and node value\
Save data file\
\
testdraw.m:\
Based on existed data, generate bifurcation plot movie and save\
\
LDS_network_data.m:\
Use Secant Continuation\
Generate data and save data file\
\
LDS_animation.m:\
Load data\
generate movies for bifurcation plots and configurations and save\
\
manualtrace_d.m:\
Based on existing data, choose new start point from mu vector and pattern vector\
Use Secant Continuation to continue in d\
\
\
Parameters:\
N: number of nodes\
p0(1): coupling strength d\
p0(2): value of mu\
steps: continuation steps\
stepsize: continuation stepsize\
m: number of nearest neighbors to right and left (error message returned if number specified exceeds possible values)\
u0: initial conditions (N*1 vector)\
k: number of initial nodes at value a\
adapt_stepsize: flag for Tangent Continuation (testbif.m)\
