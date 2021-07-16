tic

% clear all
% close all

% parameters
N = 20;	% number of nodes

p0(1) = 0.005;	% coupling strength d
% p0(2) = 2e-04; % 0.1;	% value of mu
% p0(2) = muval;

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


% % Actual steps# run
% steps = size(bd,1)-1;
% 
% 
% % Difference between consecutive steps
% % mu and norm
% dif1 = abs(diff(bd));
% dif1(:,1) = [];
% difmax1 = max(dif1);
% difmin1 = min(dif1);
% 
% % Node value
% dif2 = abs(diff(sol));
% dif2(:,1) = [];
% difmax2 = max(max(dif2));
% difmin2 = min(min(dif2));
% 
% % Print
% fprintf('Absolute step change in mu is between %.3g and %.3g.\n',difmin1(1),difmax1(1));
% fprintf('Absolute step change in norm is between %.3g and %.3g.\n',difmin1(2),difmax1(2));
% fprintf('Absolute step change in node value is between %.3g and %.3g.\n',difmin2,difmax2);
% 
% % Find Lower-Left and Upper-Right Corner
% diffbd = diff(bd);
% mat1 = diffbd(1:size(diffbd,1)-1,:);
% mat2 = diffbd(2:size(diffbd,1),:);
% compare = mat1(:,2).*mat2(:,2);
% index = find(compare < 0);
% turn_index = index + 1;
% turnbd = bd(turn_index,:);
% turnsol = sol(turn_index,:);
% [LLnorm,LLindex] = min(turnbd(:,3));
% index_maxturnbd = find(turnbd(:,2)>0.9);
% turnbd_max = turnbd(index_maxturnbd,:);
% [URnorm,URindex] = max(turnbd_max(:,3));
% 
% fprintf('Lower left corner happens at step %u with (mu,norm) = (%.3g,%.3g) \n',turnbd(LLindex,1),turnbd(LLindex,2),turnbd(LLindex,3));
% fprintf('The neighboring states of LL corner have node values (v1,v2) = (%.3g,%.3g) and (v1,v2) = (%.3g,%.3g) \n',sol(turnbd(LLindex,1),1),sol(turnbd(LLindex,1),2),sol(turnbd(LLindex,1)+2,1),sol(turnbd(LLindex,1)+2,2));
% % fprintf('Upper right happens at step %u with (mu,norm) = (%.3g,%.3g) \n',turnbd(URindex,1),turnbd(URindex,2),turnbd(URindex,3));
% % fprintf('The neighboring states of UR corner have node values (v1,v2) = (%.3g,%.3g) and (v1,v2) = (%.3g,%.3g) \n',sol(turnbd(URindex,1),1),sol(turnbd(URindex,1),2),sol(turnbd(URindex,1)+2,1),sol(turnbd(URindex,1)+2,2));
% fprintf('Upper right happens at step %u with (mu,norm) = (%.3g,%.3g) \n',turnbd_max(URindex,1),turnbd_max(URindex,2),turnbd_max(URindex,3));
% fprintf('The neighboring states of UR corner have node values (v1,v2) = (%.3g,%.3g) and (v1,v2) = (%.3g,%.3g) \n',sol(turnbd_max(URindex,1),1),sol(turnbd_max(URindex,1),2),sol(turnbd_max(URindex,1)+2,1),sol(turnbd_max(URindex,1)+2,2));

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
% str1 = sprintf('ending point near LL corner at (mu,norm) = (%.3g,%.3g)',bd(steps+1,2),bd(steps+1,3).^2);
% pt1 = [bd(steps+1,2),bd(steps+1,3).^2];
% text(pt1,str1);
% str2 = sprintf('starting point at (mu,norm) = (%.3g,%.3g)',bd(2,2),bd(2,3).^2);
% pt2 = [bd(2,2),bd(2,3).^2];
% text(pt2,str2);
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


% subplot(1,2,1)
% plot(bd(1:steps,2), bd(1:steps,3).^2, 'b.');
% hold on
% plot(bd(1,2), bd(1,3).^2, 'ro','LineWidth',3,'MarkerSize',15);
% plot(bd(end,2), bd(end,3).^2, 'mo','LineWidth',3,'MarkerSize',15);
% axis('square');
% title('Bifurcation diagram');
% xlabel('mu'); ylabel('norm u');
% 
% subplot(1,2,2)
% hold on;
% plot(spec(1:steps,1), 'r.');
% plot(spec(1:steps,2), 'b.');
% plot(spec(1:steps,3), 'm.');
% grid on;
% title('Spectra');
% xlabel('n'); ylabel('eigenvalues');
% axis('square');

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