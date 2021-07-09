LDS_Driver_all_to_all.m:\
All-to-all coupling;\
Use Secant continuation;\
3 Figures: 1. Plot bifurcation, 2. Display 9 configurations across time, and 3. Plot bifurcation with initial point circled\
Print information: absolute step change in mu, norm and node value\
\
parameters: N
N = 8;	% number of nodes

p0(1) = 0.001;	% coupling strength d
p0(2) = 0.7;	% value of mu

steps = 10000;		% #continuation steps
stepsize = -0.01;	% continuation stepsize (default: 0.0005)

% initial conditions
a = sqrt(1 + sqrt(1 - p0(2)));
b = sqrt(1 - sqrt(1 - p0(2)));
u0 = zeros(N,1);
u0(1:floor(N/2)) = a*ones(floor(N/2),1);
% u0(1) = a;
