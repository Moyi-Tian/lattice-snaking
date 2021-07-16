muvec;
inipattern;

% parameters
N = 11;	% number of nodes

p0(1) = 0.001;	% coupling strength d
p0(2) = muvec(5,2);	% value of mu

steps = 9;		% #continuation steps
stepsize = -0.0005;	% continuation stepsize (default: 0.0005)

m = 5; %nearest m-neighbors to right and left
k = 5; %number of initial nodes valued a

% initial conditions
a = sqrt(1 + sqrt(1 - p0(2)));
b = sqrt(1 - sqrt(1 - p0(2)));
u0 = inipattern(5,:);
%u0 = zeros(N,1);
%for i = 1:k
%    u0(i) = a;
%end

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
[bd,sol] = SecantContinuationLDS(sh,u0,p0,1,stepsize,steps);

figure(1)
plot(bd(1:steps,2), bd(1:steps,3).^2, '.-');
title('Bifurcation diagram');
xlabel('d'); ylabel('norm u');