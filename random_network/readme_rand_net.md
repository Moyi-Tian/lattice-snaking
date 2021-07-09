LDS_Driver_rand.m:\
Create a random symmetric adjacency matrix (emergence of each edge is 50%); \
Use Secant continuation; \
3 Figures: 1. Plot bifurcation, 2. Display 9 configurations across time, and 3. Show the graph in circular layout\
\
Driver_Rand_p.m:\
All the same except for controlling the probability of each edge to be present is p\
\
parameters:\
N: number of nodes\
p0(1): coupling strength d\
p0(2): value of mu\
steps: continuation steps\
stepsize: continuation stepsize (default: 0.0005)\
u0: initial conditions (N*1 vector)
