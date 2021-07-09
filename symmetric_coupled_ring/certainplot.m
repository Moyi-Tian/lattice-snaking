N = 121;
m = 60;

A = zeros(N);
if m > (N-1)/2
    error('Number of neighbors exceeds number of nodes');
end
if m ~= 0
    for k=1:m
        A(1,1+k) = 1; A(1,N-k+1) = 1;
    end
end
A = circulant(A(1,:),1);

col = find(bd(:,2)<0.0081);
newsol = sol(col,:);
q = size(newsol,1);
n = ceil(q/36);

for i = 1:n
    figure(i)
    for k=1:36
        subplot(6,6,k);
        g = graph(A);
        p = plot(g,'Layout','circle');
        colormap(flipud(copper));
        w = newsol(min([q (i-1)*36+k]),:);
        p.NodeCData = (w-min(w))/(max(w)-min(w));
        p.NodeLabel = {};
        p.LineStyle = 'none';
        str = sprintf('%u',col(min([q (i-1)*36+k])));
        title(['step # =',str])
    end
end