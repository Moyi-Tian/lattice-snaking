function check_isotypic(N)

for n=2:N
	ne = isotypic(2*n); % even number from 4
	ne = [ne(1:4)*2, ne(5:end)];
	if sum(ne(ne>2))>0
		disp(2*n)
	end
	
	no = isotypic(2*n-1); % odd number from 3
	no = [no(1:2)*2, no(3:end)];
	if sum(no(no>2))>0
		disp(2*n-1)
	end
end
