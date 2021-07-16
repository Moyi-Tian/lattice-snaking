tic

% clear all
% close all

% parameters
N = 20;	% number of nodes

p0(1) = 0.005;	% coupling strength d

steps = 1000000;		% #continuation steps
stepsize = -1e-04; %1e-06;	% continuation stepsize (default: 0.0005)
adapt_stepsize = 0;	% =1 adapts stepsize; =0 keeps ds unchanged

m = 1; %nearest m-neighbors to right and left
k = 20; %number of initial nodes valued a or b

% nearest m-neighbors coupling
A = zeros(N);
if m > ceil((N-1)/2)
    error('Number of neighbors exceeds number of nodes');
end
% All-to-all coupling
if m == ceil((N-1)/2)
    A = ones(N);
    A = A - diag(diag(A));
    degrees = sum(A);
    C = A - diag(degrees);
end
if m ~= 0
    for j=1:m
        A(1,1+j) = 1; A(1,N-j+1) = 1;
    end
    A = circulant(A(1,:),1);
    degrees = sum(A);
    C = A - diag(degrees);
end

% Analytic calc lower-left-corner homongeneous state evals and mu vals at bif
[mu1,mu2,eval] = hom_bif_mu(N,m,p0(1));

% Numeric calc of eval
Dlin = eig(p0(1)*C);
[V,D] = eig(p0(1)*C);

p0(2) = mu1(2);	% value of mu

% initial conditions
a = sqrt(1 + sqrt(1 - p0(2)));
b = sqrt(1 - sqrt(1 - p0(2)));
eps = 1e-03; % tolerance range for initial solution
% u0 = ustart;
u0 = zeros(N,1);
for w = 1:k
    u0(w) = b;
end
u0 = u0 + 1e-04*V(:,19);

% Check if starting point is a solution
[RHS_0,~] = LDS_RHS(u0,p0,C);
if RHS_0 < eps
    disp('Starting point is a steady state solution.');
else
    error('Starting point is NOT a steady state solution.');
end


% continuation
u0 = u0(:);
sh = @(u,p) LDS_RHS(u,p,C); % function handle for right-hand side
% [bd,sol] = SecantContinuationLDS(sh,u0,p0,2,stepsize,steps);
% [bd,sol,spec] = TangentContinuationLDS(sh,u0,p0,2,stepsize,steps,adapt_stepsize);
[bd,sol,spec] = TangentContinuationLDS_bjorn(sh,u0,p0,2,stepsize,steps,adapt_stepsize,A);


%Update actual steps
steps = size(bd,1) - 1;

% Compute values for figure
v1 = sol(:,10);
uplus = sqrt(1 + sqrt(1 - bd(:,2)));
uminus = sqrt(1 - sqrt(1 - bd(:,2)));
diff_ad_nodes = max(abs(diff(sol,1,2)),[],2);

% Plot difference in values with u-plus and u-minus
figure(3)
plot(diff_ad_nodes, 'b.');
hold on
plot(v1-uplus, 'r.');
hold on
plot(v1-uminus, 'm.');
title('Differences in Nodes'' Values');
xlabel('steps'); ylabel('value');
xlim([0 steps+1]);
legend('max of absolute differences in values between adjacent nodes in each step','v10-u+(mu)','v10-u-(mu)');
% legend('max of absolute differences in values between adjacent nodes in each step','v1-u-(mu)');

% Plot Bifurcation
figure('Position', [50 200 1200 400]);
subplot(1,2,1)
plot(bd(1:steps+1,2), bd(1:steps+1,3).^2, '.-');
hold on
plot(bd(1,2), bd(1,3).^2, 'ro','LineWidth',3,'MarkerSize',15);
hold on
plot(bd(steps+1,2), bd(steps+1,3).^2, 'mo','LineWidth',3,'MarkerSize',15);
titlename = sprintf('Bifurcation diagram with (N,m) = (%u,%u), d = %.3g, end at mu = %.5g',N,m,p0(1),bd(steps+1,2));
title(titlename);
xlabel('mu'); ylabel('norm u');

subplot(1,2,2)
plot(spec(1:steps,1), 'r.');
hold on;
plot(spec(1:steps,2), 'b.');
hold on;
plot(spec(1:steps,3), 'm.');
hold on;
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
% Lower-left-corner homogeneous branching evals data
par.mu1 = mu1;
par.mu2 = mu2;
par.eval = eval;
par.Dlin = Dlin;
par.V = V;
par.D = D;

fname = sprintf('N=%u, m=%u, k=%u, d=%.3g, mu=%.3g, s=%u, ssize=%.5g, TangentCon.mat',N,m,k,p0(1),p0(2),steps,stepsize);
save(fname,'par','bd','sol','spec');

toc
