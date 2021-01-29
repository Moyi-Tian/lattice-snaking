tic

m = 2; % number of connected neighbors to each side

for N = 26:2:36 % number of nodes
    syms U [N N]; % matrix of eigenvectors; each row is an eigenvector to be specified
    syms L [N 1]; % vector of eigenvalues
    syms d;
    
    % Eigenvectors
    for p = 1:N % iterate through the N eigenvectors
        for q = 1:N
            U(p,q) = exp(j*2*sym(pi)*p*q/N); % entries of pth eigenvector
        end
    end
    
    % Eigenvalues
    for l = 1:N % iterate through the N eigenvalues
        L(l) = 2*d*(-2+cos(2*sym(pi)*l/N)+cos(4*sym(pi)*l/N));
    end

    
    [C,ind_c,ind_l] = unique(L,'stable'); 
    % C: unique values in eigenvalues vector L shown in order of appearance
    % ind_c: show where values in C occuring in L (default: first occurrence)
    % ind_l: show where values in L occuring in C

    l_counts = accumarray(ind_l,1); % 1 as the second input to count repeated subscripts in ind_l
    value_counts = [C, l_counts]; % [each unique eigenvalue,frequency]
    
    
    % Store indices corresponding to each unique eigenvalue as a string
    len = length(value_counts);
    display_ind = strings([len,1]); 
    for t = 1:len
        display_ind(t) = strjoin(string(find(L==C(t))));
    end

    
    % Make a Table
    indices_of_occurrence = display_ind;
    eigenvalue = string(value_counts(:,1));
    frequency = string(value_counts(:,2));

    fprintf('(N,m) = (%u,%u)',N,m)
    table(indices_of_occurrence,eigenvalue,frequency)
end

toc