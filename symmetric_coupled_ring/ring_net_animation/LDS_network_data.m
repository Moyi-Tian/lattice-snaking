tic

% parameters
N = 121;	% number of nodes

p0(1) = 0.0001;	% coupling strength d
p0(2) = 0.7;	% value of mu

steps = 100000;		% #continuation steps
stepsize = -0.001;	% continuation stepsize (default: 0.0005)

m = 60; %nearest m-neighbors to right and left
k = 1; %number of initial nodes valued a

% initial conditions
a = sqrt(1 + sqrt(1 - p0(2)));
b = sqrt(1 - sqrt(1 - p0(2)));
u0 = zeros(N,1);
for i = 1:k
    u0(i) = a;
end

%nearest m-neighbors coupling
A = zeros(N);
if m > (N-1)/2
    error('Number of neighbors exceeds number of nodes');
end
if m ~= 0
    for j=1:m
        A(1,1+j) = 1; A(1,N-j+1) = 1;
    end
end
A = circulant(A(1,:),1);
degrees = sum(A);

C = A - diag(degrees);

% continuation
u0 = u0(:);
sh = @(u,p) LDS_RHS(u,p,C); % function handle for right-hand side
[bd,sol] = SecantContinuationLDS(sh,u0,p0,2,stepsize,steps);

par.N = N;
par.p0 = p0;
par.steps = steps;
par.stepsize = stepsize;
par.m = m;
par.k = k;
par.u0 = u0;
par.A = A;
par.C = C;

fname = sprintf('N=%u, m=%u, k=%u, d=%.3g, mu=%.3g, s=%u, ssize=%.5g.mat',N,m,k,p0(1),p0(2),steps,stepsize);
save(fname,'par','bd','sol');

toc