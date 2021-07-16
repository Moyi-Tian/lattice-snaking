figure(9)

pl = plot(graph(A),'Layout','circle');
colormap(flipud(copper));
pl.NodeCData = V(:,16);
pl.NodeLabel = {};
pl.LineStyle = 'none';
colorbar;
title('16th evec');
