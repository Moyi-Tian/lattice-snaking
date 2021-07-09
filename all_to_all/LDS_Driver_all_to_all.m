close all

% parameters
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

%All_to_all coupling
A = ones(N);
A = A - diag(diag(A));
degrees = sum(A);

C = A - diag(degrees);

% continuation
u0 = u0(:);
sh = @(u,p) LDS_RHS(u,p,C); % function handle for right-hand side
[bd,sol] = SecantContinuationLDS(sh,u0,p0,2,stepsize,steps);

% % postprocessing
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
figure(3)
plot(bd(2:steps+1,2), bd(2:steps+1,3).^2, '.-');
hold on
plot(bd(1,2), bd(1,3).^2, 'ro','LineWidth',3,'MarkerSize',15);
title('Bifurcation diagram');
xlabel('mu'); ylabel('norm u');