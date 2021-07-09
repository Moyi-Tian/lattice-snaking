% parameters
N = 30;	% number of nodes

p0(1) = 0.01;	% coupling strength d
p0(2) = 0.7;	% value of mu

steps = 100;		% #continuation steps
stepsize = -0.001;	% continuation stepsize (default: 0.0005)

% initial conditions
a = sqrt(1 + sqrt(1 - p0(2)));
b = sqrt(1 - sqrt(1 - p0(2)));
u0 = zeros(N,1);
u0(1) = a;

%Random laplacian matrix C
degrees = zeros(1,N);
while min(degrees) == 0
    A = round(rand(N));
    A = triu(A) + triu(A,1)';
    A = A - diag(diag(A));
    degrees = sum(A);   
end

C = A - diag(degrees);

% continuation
u0 = u0(:);
sh = @(u,p) LDS_RHS_rand(u,p,C); % function handle for right-hand side
[bd,sol] = SecantContinuationLDS(sh,u0,p0,2,stepsize,steps);

% postprocessing
figure(1)
plot(bd(4:steps,2), bd(4:steps,3).^2, '.-');
title('Bifurcation diagram');
xlabel('mu'); ylabel('norm u');

figure(2)
for k=1:9
	subplot(3,3,k);
    p = plot(graph(A),'Layout','circle');
    colormap(flipud(copper));
    w = sol(k*floor(steps/9),:);
    p.NodeCData = (w-min(w))/(max(w)-min(w));
    p.NodeLabel = {};
    p.LineStyle = 'none';
    colorbar;
end

figure(3)
plot(graph(A),'Layout','circle');