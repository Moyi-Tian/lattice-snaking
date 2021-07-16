tic

clear all

% parameters
N = 20;	% number of nodes

p0(1) = 0.05;	% coupling strength d
p0(2) = 0.7;	% value of mu
%p0(2) = muval;

steps = 10000;		% #continuation steps
stepsize = -0.01;	% continuation stepsize (default: 0.0005)
adapt_stepsize = 1;	% =1 adapts stepsize; =0 keeps ds unchanged

m = 9; %nearest m-neighbors to right and left
k = 1; %number of initial nodes valued a

% initial conditions
a = sqrt(1 + sqrt(1 - p0(2)));
b = sqrt(1 - sqrt(1 - p0(2)));
%u0 = initial;
u0 = zeros(N,1);
for w = 1:k
    u0(w) = a;
end

% %nearest m-neighbors coupling
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
% [bd,sol] = SecantContinuationLDS(sh,u0,p0,2,stepsize,steps);
[bd,sol,spec] = TangentContinuationLDS(sh,u0,p0,2,stepsize,steps,adapt_stepsize);

% Actual steps# run
steps = size(bd,1)-1;


% Difference between consecutive steps
% mu and norm
dif1 = abs(diff(bd));
dif1(:,1) = [];
difmax1 = max(dif1);
difmin1 = min(dif1);

% Node value
dif2 = abs(diff(sol));
dif2(:,1) = [];
difmax2 = max(max(dif2));
difmin2 = min(min(dif2));

% Print
fprintf('Absolute step change in mu is between %.3g and %.3g.\n',difmin1(1),difmax1(1));
fprintf('Absolute step change in norm is between %.3g and %.3g.\n',difmin1(2),difmax1(2));
fprintf('Absolute step change in node value is between %.3g and %.3g.\n',difmin2,difmax2);


% Plot Bifurcation
figure(2)
plot(bd(2:steps+1,2), bd(2:steps+1,3).^2, '.-');
hold on
plot(bd(1,2), bd(1,3).^2, 'ro','LineWidth',3,'MarkerSize',15);
title('Bifurcation diagram');
xlabel('mu'); ylabel('norm u');

subplot(1,2,1)
plot(bd(1:steps,2), bd(1:steps,3).^2, 'b.');
hold on
plot(bd(1,2), bd(1,3).^2, 'ro','LineWidth',3,'MarkerSize',15);
plot(bd(end,2), bd(end,3).^2, 'mo','LineWidth',3,'MarkerSize',15);
axis('square');
title('Bifurcation diagram');
xlabel('mu'); ylabel('norm u');

subplot(1,2,2)
hold on;
plot(spec(1:steps,1), 'r.');
plot(spec(1:steps,2), 'b.');
plot(spec(1:steps,3), 'm.');
grid on;
title('Spectra');
xlabel('n'); ylabel('eigenvalues');
axis('square');

% Save Data
par.N = N;
par.p0 = p0;
par.steps = steps;
par.stepsize = stepsize;
par.adapt_stepsize = adapt_stepsize;
par.m = m;
par.k = k;
par.u0 = u0;
par.A = A;
par.C = C;

fname = sprintf('N=%u, m=%u, k=%u, d=%.3g, mu=%.3g, s=%u, ssize=%.5g, TangentCon.mat',N,m,k,p0(1),p0(2),steps,stepsize);
save(fname,'par','bd','sol','spec');

toc