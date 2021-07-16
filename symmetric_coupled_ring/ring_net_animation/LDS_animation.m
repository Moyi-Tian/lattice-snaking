tic

clear
clc
close all
%% Load data

%choose parameters to load
N = 11;
m = 4;	% number of nodes
k = 1; % number of initial nodes valued a
p0(1) = 0.0001;	% coupling strength d
p0(2) = 0.7;	% value of mu
steps = 10000000;		% #continuation steps
stepsize = -0.00001;	% continuation stepsize (default: 0.0005)
 %nearest m-neighbors to right and left
filename = sprintf('N=%u, m=%u, k=%u, d=%.3g, mu=%.3g, s=%u, ssize=%.5g.mat',N,m,k,p0(1),p0(2),steps,stepsize);

load(filename);
%% Draw

n = 500; %number of frames
fName = sprintf('frame #=%u, N=%u, m=%u, k=%u, d=%.3g, mu=%.3g',n,N,m,k,p0(1),p0(2));
figure('Name',fName,'NumberTitle','off','Position', [100 100 1500 450]);

q = ceil(max(sol,[],'all')*10)/10; %lifted max in sol
for i = 1:n
    clf
    j = (i-1)*floor(par.steps/n)+1; %kth step
    
    % bifuracation diagram
    subplot(1,2,1)
    plot(bd(4:par.steps,2), bd(4:par.steps,3).^2, '.-');
    hold on
    plot(bd(j,2), bd(j,3).^2, 'ro','LineWidth',3,'MarkerSize',15);
    title('Bifurcation diagram');
    xlabel('mu'); ylabel('norm u');
    
    % pattern diagram
    subplot(1,2,2)
    w = sol(j,:);    
    g = graph(par.A);
    p = plot(g,'Layout','circle');
    colormap(flipud(copper));
    p.NodeCData = w/q;
    p.NodeLabel = {};
    p.LineStyle = 'none';
    colorbar;
    caxis([0 q]); %colorbar range
    str1 = sprintf('%u',j);
    str2 = sprintf('%u',par.steps);
    title(['step # =  ',str1,' out of ',str2,'  steps'])
    
    %store frame
    movieVector(i) = getframe(gcf);

end
%% Save movie

r = 10; %Frame Rate
VideoName = sprintf('frame #=%u,FrameRate=%u, N=%u, m=%u, k=%u, d=%.3g, mu=%.3g, s=%u',n,r,N,m,k,p0(1),p0(2),steps);
myWriter = VideoWriter(VideoName,'MPEG-4');

myWriter.FrameRate = r;

open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);

toc