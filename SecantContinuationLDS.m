function [bd,sol] = SecantContinuationLDS(fhandle,u0,pars,icp,ds,nmx)

  % Initialise variables
  ndim = size(u0,1); bd = zeros(nmx+1,3); sol = zeros(nmx+1,ndim);
  v0 = zeros(ndim+1,1); v1 = zeros(ndim+1,1); 

  % Options to the nonlinear solver
  opts = optimset('Display','off','Jacobian','on');

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
  for n = 2:nmx

    % Prediction in the secant direction
    sec = (v1-v0)/norm(v1-v0); v = v1 + abs(ds) * sec;

    % Correction with Newton iteration
    v = fsolve( @(v) SecantCorrector(v), v, opts );

    % Book-keeping
    v0 = v1; v1 = v; 
    bd(n+1,:)  = [n v(ndim+1) norm(v(1:ndim))];
    sol(n+1,:) = v(1:ndim)';

    % Print 
    fprintf('%9d %14.4e %16.4e\n',bd(n+1,:));

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
end
