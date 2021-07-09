function [F,J] = LDS_RHS_rand(u,p,C)
  
% rename parameters
d = p(1);
mu = p(2);
N = size(u,1);

% right-hand side
F = d*C*u -mu*u + 2*u.^3 - u.^5;

% Jacobian
if nargout > 1
	J = d*C + spdiags(-mu + 6*u.^2 - 5*u.^4, 0, N, N);
end

end