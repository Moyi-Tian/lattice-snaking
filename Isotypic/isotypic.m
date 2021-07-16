function ndim = isotypic(n)

ndim = [];

if rem(n,2)==0
	% n is even
    v = zeros(1,n);
    v(2) = 1;
    A = circulant(v); % rotation
    R = zeros(n); % reflection
    for l = 1:n/2
        R(l,n+1-l) = 1;
        R(n+1-l,l) = 1;
    end

	% one-dimensional representations
	ni = 1; % degree

	Pi = zeros(n);
	for j=0:(n-1)
		Pi = Pi + A^j + A^j*R;
	end
	Pi = ni/(2*n)*Pi;
	ndim = [ndim, trace(Pi)];
	
	Pi = zeros(n);
	for j=0:(n-1)
		Pi = Pi + A^j - A^j*R;
	end
	Pi = ni/(2*n)*Pi;
	ndim = [ndim, trace(Pi)];

	Pi = zeros(n);
	for j=0:(n-1)
		Pi = Pi + (-1)^j*(A^j + A^j*R);
	end
	Pi = ni/(2*n)*Pi;
	ndim = [ndim, trace(Pi)];
	
	Pi = zeros(n);
	for j=0:(n-1)
		Pi = Pi + (-1)^j*(A^j - A^j*R);
	end
	Pi = ni/(2*n)*Pi;
	ndim = [ndim, trace(Pi)];

	% two-dimensional representations
	ni = 2;
    for k=1:(n-2)/2 % iterate over representations
		Pi = zeros(n);
		for j=0:n-1 % # of rotation
			char = 2*cos(2*pi*j*k/n);
			Pi = Pi + A^j*char;
		end
		Pi = ni/(2*n)*Pi;
		ndim = [ndim, trace(Pi)];
	end
else
	% n is odd
	v = zeros(1,n);
	v(2) = 1;
	A = circulant(v); % rotation
	R = zeros(n); % reflection (w.r.t. axis passing first node)
	R(1,1) = 1;
	for l = 1:(n-1)/2
		R(l+1,n+1-l) = 1;
		R(n+1-l,l+1) = 1;
	end

	% one-dimensional representations
	ni = 1;

	Pi = zeros(n);
	for j=0:(n-1)
		Pi = Pi + A^j + A^j*R;
	end
	Pi = ni/(2*n)*Pi;
	ndim = [ndim, trace(Pi)];
	
	Pi = zeros(n);
	for j=0:(n-1)
		Pi = Pi + A^j - A^j*R;
	end
	Pi = ni/(2*n)*Pi;
	ndim = [ndim, trace(Pi)];

	% two-dimensional representations
	ni = 2;
	for k=1:(n-1)/2 % iterate over representations
		Pi = zeros(n);
		for j=0:(n-1) % # of ratations
			char = 2*cos(2*pi*j*k/n);
			Pi = Pi + A^j*char;
		end
		Pi = ni/(2*n)*Pi;
		ndim = [ndim, trace(Pi)];
	end
end

ndim = round(ndim);

if sum(ndim)~=n
	disp('dimensions do not match')
end


