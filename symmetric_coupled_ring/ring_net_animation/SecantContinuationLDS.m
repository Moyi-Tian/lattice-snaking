function [bd,sol] = SecantContinuationLDS(fhandle,u0,pars,icp,ds,nmx)

  eflag = [0,0];
  chkst = 1000; % step number between plots
  %jump = 1000; % steps jumped over in period check
  %eps1 = 4e-6; % epsilon range for mu
  %eps2 = 5e-6; % epsilon range for norm or node value
  

  % Initialise variables
  ndim = size(u0,1); bd = zeros(nmx+1,3); sol = zeros(nmx+1,ndim);
  v0 = zeros(ndim+1,1); v1 = zeros(ndim+1,1); 

  % Options to the nonlinear solver
  %opts = optimset('Display','off','Jacobian','on');
  opts = optimset('Display','off','Jacobian','on','MaxIter',10000);

  % Prepare screen output
  fprintf('%9s %14s %16s\n','STEP','PAR','2-NORM');

 
  % Converge initial guess
  v0(ndim+1) = pars(icp);
  v0(1:ndim) = fsolve( @(u) fhandle(u,pars), u0, opts );
  bd(1,:)  = [0 v0(ndim+1) norm(v0(1:ndim))];
  sol(1,:) = v0(1:ndim)';
  fprintf('%9d %14.4e %16.4e\n',bd(1,:));

  % Poor-man continuation step
  pars(icp)  = pars(icp) + ds/(10*sqrt(ndim)); %divide by norm
  v1(ndim+1) = pars(icp);
  v1(1:ndim) = fsolve( @(u) fhandle(u,pars), v0(1:ndim), opts );
  bd(2,:)  = [1 v0(ndim+1) norm(v0(1:ndim))];
  sol(2,:) = v1(1:ndim)';
  fprintf('%9d %14.4e %16.4e\n',bd(2,:));

  % Start secant continuation
  ct2 = 1;
  for n = 2:nmx
      
    %ds = ds/2;

    % Prediction in the secant direction
    sec = (v1-v0)/norm(v1-v0); v = v1 + abs(ds) * sec;

    % Correction with Newton iteration
    [v,fval,exitflag,output] = fsolve( @(v) SecantCorrector(v), v, opts );
    
    % Adjust mu when mu becomes negative
    %count = 0;
    %while v(ndim+1)<0
    %   v(ndim+1) = count+0.001;
    %   sec = (v1-v0)/norm(v1-v0); v = v1 + abs(ds) * sec;
    %   [v,fval,exitflag,output] = fsolve( @(v) SecantCorrector(v), v, opts );
    %   count = count+0.001;
    %end
    
    if exitflag <= 0
        eflag = [eflag;n,exitflag];
    end

    % Book-keeping
    v0 = v1; v1 = v; 
    bd(n+1,:)  = [n v(ndim+1) norm(v(1:ndim))];
    sol(n+1,:) = v(1:ndim)';
    
    % Plot after every chkst steps
    if n+1 >= chkst
        ct1 = floor(n/chkst);
        if ct2 < ct1
            figure(1)
            plot(bd(1+chkst*(ct1-1):chkst*ct1,2), bd(1+chkst*(ct1-1):chkst*ct1,3).^2, 'b.');
            hold on
            
            fprintf('%9d %14.4e %16.4e\n',bd(n+1,:));
        end
        ct2 = ct1;
    end

    % Print 
    %fprintf('%9d %14.4e %16.4e\n',bd(n+1,:));
    
    % Stop after 1 period
    %if n+1 > jump
        % Compare mu and all solution points
        %if bd(n+1,2)>bd(1,2)-eps1 && bd(n+1,2)<bd(1,2)+eps1 && all(sol(n+1,:)>sol(1,:)-eps2) && all(sol(n+1,:)<sol(1,:)+eps2)
    
        % Compare mu and norm
        %if bd(n+1,2)>bd(1,2)-eps1 && bd(n+1,2)<bd(1,2)+eps1 && bd(n+1,3)>bd(1,3)-eps2 && bd(n+1,3)<bd(1,3)+eps2
            %disp('Reached 1 period. Returned to initial point.');
            %eflag(1,:) = [];
            %eflag
            %bd(n+2:nmx+1,:) = [];
            %sol(n+2:nmx+1,:) = [];
            %return
        %end
    %end
  end

  function [G,DG] = SecantCorrector(v)

    % Compute F
    pars(icp) = v(ndim+1); F = fhandle(v(1:ndim),pars);

    % Extended system
    G = [F; sec' * (v-v1) - abs(ds)];
    
    if nargout > 1

      % Jacobian of F
      [F,DFDU] = fhandle(v(1:ndim),pars);

      % Derivative with respect to p
      epsi = 1e-6; pars(icp) = v(ndim+1) + epsi;
      DF = fhandle(v(1:ndim),pars); DFDP = (DF - F)/epsi;

      % Jacobian of the extended system
      DG = [DFDU, DFDP; sec'];

    end

  end

eflag(1,:) = [];
eflag
end
