tic

%% Draw

n = 600; %number of frames
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
    plot(bd(j,2), bd(j,3).^2, 'ro','LineWidth',3,'MarkerSize',10);
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