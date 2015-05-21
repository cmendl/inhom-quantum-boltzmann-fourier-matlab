function Un1 = slopeLimStepMaxwell(h,k,A,alpha,fluxL,fluxR,UmaxwL,UmaxwR,Un)
%SLOPELIMSTEPMAXWELL - Numerically solve u_t + A u_x = 0 with Maxwell boundary conditions using a slope limiter method
%
%    Un1 = SLOPELIMSTEPMAXWELL(h,k,A,alpha,fluxL,fluxR,UmaxwL,UmaxwR,Un)
%
%    Input:    h        finite volume cell size
%              k        time step
%              A        the 'A' matrix
%              alpha    accommodation coefficient for Maxwell reflection operator
%              fluxL    predetermined total outgoing flux at left  boundary
%              fluxR    predetermined total outgoing flux at right boundary
%              UmaxwL   normalized left  boundary state for Maxwell boundary condition
%              UmaxwR   normalized right boundary state for Maxwell boundary condition
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

% select positive and negative parts
selectPos = R*diag(lambda >= 0)/R;
selectNeg = R*diag(lambda <  0)/R;

% flip positive and negative parts
flipPos = R*flip(diag(lambda >= 0))/R;
flipNeg = R*flip(diag(lambda <  0))/R;

% Maxwell boundary conditions for incoming part
% U^n_{-1}
Um1 = selectNeg*Un(:,1)   + (1-alpha)*flipNeg*Un(:,1)   + alpha*fluxL*selectPos*UmaxwL;
% U^n_{N+1}
UN1 = selectPos*Un(:,end) + (1-alpha)*flipPos*Un(:,end) + alpha*fluxR*selectNeg*UmaxwR;

% pad copies of Maxwell boundary states on the left and right
Un1 = slopeLimStep(h,k,A,[Um1,Um1,Un,UN1,UN1]);
Un1 = Un1(:,3:end-2);	% remove auxiliary boundary values
