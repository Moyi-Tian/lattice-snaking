% parameters
N = 121;	% number of nodes

p0(1) = 0.001;	% coupling strength d
p0(2) = 0.7;	% value of mu

steps = 100000;		% #continuation steps
stepsize = -0.001;	% continuation stepsize (default: 0.0005)

m = 60; %nearest m-neighbors to right and left

% initial conditions
a = sqrt(1 + sqrt(1 - p0(2)));
b = sqrt(1 - sqrt(1 - p0(2)));
u0 = zeros(N,1);
u0(1) = a;

%nearest m-neighbors coupling
A = zeros(N);
if m > (N-1)/2
    error('Number of neighbors exceeds number of nodes');
end
if m ~= 0
    for k=1:m
        A(1,1+k) = 1; A(1,N-k+1) = 1;
    end
end
A = circulant(A(1,:),1);
degrees = sum(A);

C = A - diag(degrees);

% continuation
u0 = u0(:);
sh = @(u,p) LDS_RHS(u,p,C); % function handle for right-hand side
[bd,sol] = SecantContinuationLDS(sh,u0,p0,2,stepsize,steps);

% postprocessing
figure(1)
plot(bd(4:steps,2), bd(4:steps,3).^2, '.-');
title('Bifurcation diagram');
xlabel('mu'); ylabel('norm u');

figure(2)
for k=1:9
	subplot(3,3,k);
    g = graph(A);
    p = plot(g,'Layout','circle');
    colormap(flipud(copper));
    w = sol(k*floor(steps/9),:);
    p.NodeCData = (w-min(w))/(max(w)-min(w));
    p.NodeLabel = {};
    p.LineStyle = 'none';
    colorbar;
end