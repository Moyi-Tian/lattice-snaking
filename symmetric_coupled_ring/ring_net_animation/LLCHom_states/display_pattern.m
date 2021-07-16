figure(9)
% for k=1:9
% 	subplot(3,3,k);
%     pl = plot(graph(A),'Layout','circle');
%     colormap(flipud(copper));
%     w = sol(k*floor(steps/9),:);
%     pl.NodeCData = (w-min(w))/(max(w)-min(w));
%     pl.NodeLabel = {};
%     pl.LineStyle = 'none';
%     colorbar;
% end

% pl = plot(graph(A),'Layout','circle');
% colormap(flipud(copper));
% w = sol(12558,:); %+19th&-19th,+stepsize,s = 9594;+19th&-19th,-stepsize,s =12558 
% pl.NodeCData = w/max(max(sol));
% pl.NodeLabel = {};
% pl.LineStyle = 'none';
% colorbar;
% title('hom branch N=20, m=1, +19th evec,ssize=-1e-04, mu = 0.5');

pl = plot(graph(A),'Layout','circle');
colormap(flipud(copper));
%w = sol(9594,:); %+19th&-19th,+stepsize,s = 9594; 
pl.NodeCData = V(:,16);
% pl.NodeCData = w/max(max(sol));
pl.NodeLabel = {};
pl.LineStyle = 'none';
colorbar;
title('16th evec');

