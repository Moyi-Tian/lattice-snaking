function [mu1,mu2,eval] = hom_bif_mu(N,m,d)

if m > ceil((N-1)/2)
    error('Number of neighbors exceeds number of nodes');
end

if m == ceil((N-1)/2) && mod(N,2) == 0
    error('The network is an all-to-all coupling. This function does not work for all-to-all with even number of nodes.');
end

c = -m*d*ones(N,1);

for j = 0:N-1
    for i = 1:m
        c(j+1) = c(j+1) + d*cos(2*pi*i*j/N);
    end
end

eval = 2*c;
mu1 = 1 - (c+1)/2 - sqrt(2*c+1)/2;
mu2 = 1 - (c+1)/2 + sqrt(2*c+1)/2;