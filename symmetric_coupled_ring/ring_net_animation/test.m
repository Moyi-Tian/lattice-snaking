tic

clear all

%%Prepare Data
% parameters
N = 11;	% number of nodes

p0(1) = 0.001;	% coupling strength d
p0(2) = 0.7;	% value of mu
%p0(2) = muval;

steps = 10000;		% #continuation steps
stepsize = -0.01;	% continuation stepsize (default: 0.0005)

m = 4; %nearest m-neighbors to right and left
k = 1; %number of initial nodes valued a

% initial conditions
a = sqrt(1 + sqrt(1 - p0(2)));
b = sqrt(1 - sqrt(1 - p0(2)));
%u0 = initial;
u0 = zeros(N,1);
for w = 1:k
    u0(w) = a;
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

%% Draw

n = 300; %number of frames
%fName = sprintf('frame #=%u, N=%u, m=%u, k=%u, d=%.3g, mu=%.3g',n,N,m,k,p0(1),p0(2));
%figure('Name',fName,'NumberTitle','off','Position', [300 300 500 450]);
figure('Position', [350 350 650 450]);

q = ceil(max(sol,[],'all')*10)/10; %lifted max in sol
for i = 1:n
    clf
    j = (i-1)*floor(steps/n)+1; %kth step
    
    % bifuracation diagram
    plot(bd(1:steps,2), bd(1:steps,3).^2, '.-');
    hold on
    plot(bd(j,2), bd(j,3).^2, 'ro','LineWidth',3,'MarkerSize',15);
    title('Bifurcation diagram');
    xlabel('mu'); ylabel('norm u');
    
    %store frame
    movieVector(i) = getframe(gcf);

end
%% Save movie

r = 5; %Frame Rate
VideoName = sprintf('frame #=%u,FrameRate=%u, N=%u, m=%u, k=%u, d=%.3g, mu=%.3g, step=%u, stepsize=%.3g',n,r,N,m,k,p0(1),p0(2),steps,stepsize);
myWriter = VideoWriter(VideoName,'MPEG-4');

myWriter.FrameRate = r;

open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);

toc