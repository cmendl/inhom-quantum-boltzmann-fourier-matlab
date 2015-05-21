function Un1 = godunovStep(h,k,A,Un)
%GODUNOVSTEP - Numerically solve u_t + A u_x = 0 with periodic boundary conditions using Godunov's method
%
%    Un1 = GODUNOVSTEP(h,k,A,Un)
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

% U^n_{j+1} using periodic boundary conditions
Unxt = circshift(Un,[0,-1]);
% U^n_{j-1} using periodic boundary conditions
Uprv = circshift(Un,[0,1]);

r = k/(2*h);

% Eq. (13.17)
Un1 = Un ...
	- r*A   *(Unxt - Uprv) ...
	+ r*Aabs*(Unxt - 2*Un + Uprv);
