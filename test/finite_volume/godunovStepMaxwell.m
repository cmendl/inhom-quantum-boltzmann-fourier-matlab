function Un1 = godunovStepMaxwell(h,k,A,alpha,UmaxwL,UmaxwR,Un)
%GODUNOVSTEPMAXWELL - Numerically solve u_t + A u_x = 0 with Maxwell boundary conditions using Godunov's method
%
%    Un1 = GODUNOVSTEPMAXWELL(h,k,A,alpha,UmaxwL,UmaxwR,Un)
%
%    Input:    h        finite volume cell size
%              k        time step
%              A        the 'A' matrix
%              alpha    accommodation coefficient for Maxwell reflection operator
%              UmaxwL   left  boundary state for Maxwell boundary condition
%              UmaxwR   right boundary state for Maxwell boundary condition
%              Un       state variable (m x N matrix, with N the number of finite volumes)
%
%    Output:   Un1      state variable at next time step
%
%    Reference: Randall LeVeque. Numerical Methods for Conservation Laws (1992)

% A = R*diag(lambda)*R^{-1}
[R,lambda] = eig(A);
lambda = diag(lambda);

% sort eigenvalues
[lambda,IX] = sort(lambda,'descend');
R = R(:,IX);

% consistency checks
assert(all(isreal(lambda)));
% eigenvalues must appear in pairs with opposite sign
assert(all(abs(lambda + flip(lambda)) < 8*eps*length(lambda)));
% first half is assumed to be positive
assert(all(lambda(1:end/2)>0));
% second half is assumed to be negative
assert(all(lambda(end/2+1:end)<0));

% Eq. (10.61)
Apos = R*diag(max(lambda,0))/R;
Aneg = R*diag(min(lambda,0))/R;

Aabs = R*diag(abs(lambda))/R;

% normalize prescribed boundary states
UmaxwL = UmaxwL / sum( Apos*UmaxwL);
UmaxwR = UmaxwR / sum(-Aneg*UmaxwR);

% select positive and negative parts
selectPos = R*diag(lambda >= 0)/R;
selectNeg = R*diag(lambda <  0)/R;

% flip positive and negative parts
flipPos = R*flip(diag(lambda >= 0))/R;
flipNeg = R*flip(diag(lambda <  0))/R;

% Maxwell boundary conditions for incoming part
% U^n_{-1}
Um1 = selectNeg*Un(:,1)   + (1-alpha)*flipNeg*Un(:,1)   + alpha*sum(-Aneg*Un(:,1))  *selectPos*UmaxwL;
% U^n_{N+1}
UN1 = selectPos*Un(:,end) + (1-alpha)*flipPos*Un(:,end) + alpha*sum( Apos*Un(:,end))*selectNeg*UmaxwR;

% U^n_{j+1}
Unxt = [Un(:,2:end),UN1];
% U^n_{j-1}
Uprv = [Um1,Un(:,1:end-1)];

r = k/(2*h);

% Eq. (13.17)
Un1 = Un ...
	- r*A   *(Unxt - Uprv) ...
	+ r*Aabs*(Unxt - 2*Un + Uprv);
