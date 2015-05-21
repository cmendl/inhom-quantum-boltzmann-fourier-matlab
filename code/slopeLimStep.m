function Un1 = slopeLimStep(h,k,A,Un)
%SLOPELIMSTEP - Numerically solve u_t + A u_x = 0 with periodic boundary conditions using a slope limiter method
%
%    Un1 = SLOPELIMSTEP(h,k,A,Un)
%
%    Input:    h        finite volume cell size
%              k        time step
%              A        the 'A' matrix
%              Un       state variable (m x N matrix, with N the number of finite volumes)
%
%    Output:   Un1      state variable at next time step
%
%    Reference: Randall LeVeque. Numerical Methods for Conservation Laws (1992)

% Eq. (13.16)
[R,lambda] = eig(A);
Aabs = R*abs(lambda)/R;		% R*diag(abs(lambda))*R^{-1}

% extract diagonal part
lambda = diag(lambda);
% nu_p = k*lambda_p/h, Eq. (16.36)
nu = k*lambda/h;

% for slope-limiter term
S    = R*diag(0.5*    lambda .*(sign(nu)-nu))/R;
Sabs = R*diag(0.5*abs(lambda).*(sign(nu)-nu))/R;

% Eq. (16.52)
minmod = @(a,b) 0.5 * (sign(a) + sign(b)) .* min(abs(a),abs(b));

% U^n_{j+1} using periodic boundary conditions
Unxt = circshift(Un,[0,-1]);
% U^n_{j-1} using periodic boundary conditions
Uprv = circshift(Un,[0,1]);

% Eq. (16.34)
% alpha_j = R^{-1}*(U^n_{j+1} - U^n_j)
alpha = R\(Unxt - Un);

% h * beta_j, Eq. (16.56)
h_beta = minmod(alpha,circshift(alpha,[0,1]));

% sum of h*sigma_j in Eq. (16.58) with respect to p
sn = R*h_beta;
snxt = circshift(sn,[0,-1]);
sprv = circshift(sn,[0,1]);

% first term in square brackets of Eq. (16.57) is basically Godunov step in Eq. (13.17)
Un1 = Un + k/(2*h)*( ...
	- A*(Unxt - Uprv) + Aabs*(Unxt - 2*Un + Uprv) ...
	- S*(snxt - sprv) + Sabs*(snxt - 2*sn + sprv));
