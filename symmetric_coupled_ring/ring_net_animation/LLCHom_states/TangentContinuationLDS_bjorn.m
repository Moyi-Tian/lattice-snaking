function [bd,sol,spec] = TangentContinuationLDS_bjorn(fhandle,u0,pars,icp,ds,nmx,flag,adjmtx)

  % Initialise variables
  chkst = 5000; % step number between plots
  stopt = 0.8;
  counter = 0; % count how many times solution passes vertical line of stopt
  presetnum = 1; % stop when counter > presetnum
  
  ndim = size(u0,1);
  bd = zeros(nmx+1,3);
  sol = zeros(nmx+1,ndim);
  spec = zeros(nmx+1,ndim);
  v1 = zeros(ndim+1,1);
  ds1 = ds;

  % Options to the nonlinear solver
  opts = optimset('Display','off','Jacobian','on');

  % Converge initial guess
  v1(ndim+1) = pars(icp);
  [v1(1:ndim),~,exitflag,~] = fsolve( @(u) fhandle(u,pars), u0, opts );  
  bd(1,:)  = [0 v1(ndim+1) norm(v1(1:ndim))];
  sol(1,:) = v1(1:ndim)';
  
  % if fsolve does not find a solution, print
  if exitflag<=0
		fprintf('%2d %14.4e\n', exitflag, bd(1,2));
  end

  %  fprintf('%9s %14s %16s\n','STEP','PAR','2-NORM');
  %  fprintf('%9d %14.4e %16.4e\n',bd(1,:));

  % Predictor
  pars(icp) = v1(ndim+1);
  [F,DFDU] = fhandle(v1(1:ndim),pars);
  spec(1,:) = sort(real(eig(full(DFDU))),'descend');
  sec = null(full([DFDU, -v1(1:ndim)]));
  sec = sign(sec(ndim+1))*sign(ds)*sec;
  if flag==1
	  ds1 = ds*min(1,500/condest(DFDU));
  end
  v = v1 + abs(ds1)*sec;

  % Start tangent continuation
  ct2 = 1;
  for n = 1:nmx

	% Corrector
    [v,~,exitflag,~] = fsolve( @(v) SecantCorrector(v), v, opts );

    % Book-keeping
    bd(n+1,:)  = [n v(ndim+1) norm(v(1:ndim))];
    sol(n+1,:) = v(1:ndim)';
	v1 = v;
    
    % if fsolve does not find a solution, print
    if exitflag<=0
		fprintf('%2d %14.4e\n', exitflag, bd(n+1,2));
	end
    
    % Plot after every chkst steps
    if n+1 >= chkst
        ct1 = floor(n/chkst);
        if ct2 < ct1
            figure(1)
            plot(bd(1+chkst*(ct1-1):chkst*ct1,2), bd(1+chkst*(ct1-1):chkst*ct1,3).^2, 'b.');
            hold on
            
%             % Pattern Diagram
%             q = ceil(max(sol,[],'all')*10)/10; %lifted max in sol
%             figure(2)
%             w = sol(n+1,:);    
%             g = graph(adjmtx);
%             p = plot(g,'Layout','circle');
%             colormap(flipud(copper));
%             p.NodeCData = w/q;
%             p.NodeLabel = {};
%             p.LineStyle = 'none';
%             colorbar;
%             caxis([0 q]); %colorbar range
%             str1 = sprintf('%u',n+1);
%             str2 = sprintf('%u',nmx);
%             title(['step # =  ',str1,' out of ',str2,'  steps'])
            
            fprintf('%9d %14.4e %16.4e\n',bd(n+1,:));
        end
        ct2 = ct1;
    end
	
	% Stopping criterion
	par_stop = stopt;
	if (bd(n,2)-par_stop)*(bd(n+1,2)-par_stop)<=0
		counter = counter+1;
		if counter > presetnum
            bd(n+2:nmx+1,:) = [];
            sol(n+2:nmx+1,:) = [];
			break;
		end
	end
	
%  	if ((bd(n,2)-0.1)*(bd(n+1,2)-0.1)<=0)
%  		data.sol = sol(n+1,:);
%  		data.par = bd(n+1,2);
%  		save('data.mat','data');
%  		return
%  	end

	% Predictor
    pars(icp) = v1(ndim+1);
	[F,DFDU] = fhandle(v1(1:ndim),pars);
	spec(n+1,:) = sort(real(eig(full(DFDU))),'descend');
    sec1 = null(full([DFDU, -v1(1:ndim)]));
	sec = sign(sec'*sec1)*sec1;
	if flag==1
		ds1 = ds*min(1,1.e3/condest(DFDU));
		ds1 = max(ds1,1.e-4);
	end
	v = v1 + abs(ds1)*sec;
  
	% Output
%	b = sqrt(1-sqrt(1-pars(2)));
%	fprintf( '%9d %14.4e %16.4e %14.4e %14.4e %14.4e\n', ...
%		bd(n+1,:), sol(n+1,1), norm(diff(sol(n+1,:))), abs(sol(n+1,1)-b) );

  end
  
  function [G,DG] = SecantCorrector(v)

    % Compute F
    pars(icp) = v(ndim+1);
	F = fhandle(v(1:ndim),pars);

    % Extended system
    G = [F; sec'*(v-v1) - abs(ds1)];
    
    if nargout > 1
      % Jacobian of F
      [F,DFDU] = fhandle(v(1:ndim),pars);

      % Derivative with respect to p
%     epsi = 1e-8;
% 	  parsn = pars;
% 	  parsn(icp) = v(ndim+1) + epsi;
%       DF = fhandle(v(1:ndim),parsn);
% 	  DFDP = (DF-F)/epsi;
	  
	  DFDP = -v(1:ndim);

      % Jacobian of the extended system
      DG = [DFDU, DFDP; sec'];
    end

  end
end
