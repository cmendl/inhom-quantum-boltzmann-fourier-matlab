function Un1 = godunovStepDirichlet(h,k,A,Un)
%GODUNOVSTEPDIRICHLET - Numerically solve u_t + A u_x = 0 with Dirichlet boundary conditions using Godunov's method
%
%    Un1 = GODUNOVSTEPDIRICHLET(h,k,A,Un)
%
%    Input:    h        finite volume cell size
%              k        time step
%              A        the 'A' matrix
%              Un       state variable (m x N matrix, with N the number of finite volumes)
%
%    Output:   Un1      state variable at next time step
%
%    Reference: Randall LeVeque. Numerical Methods for Conservation Laws (1992)

% pad copies of boundary values on the left and right
% and run standard Godunov step
Un1 = godunovStep(h,k,A,[Un(:,1),Un,Un(:,end)]);

% remove padded copies
Un1 = Un1(:,2:end-1);

% A = R*diag(lambda)*R^{-1}
[R,lambda] = eig(A);
lambda = diag(lambda);

selectPos = R*diag(lambda >= 0)/R;		% select positive part
selectNeg = R*diag(lambda <  0)/R;		% select negative part

% substitute incoming part by given Dirichlet boundary states
% and keep the outgoing part
Un1(:,1)   = selectPos*Un(:,1)   + selectNeg*Un1(:,1);
Un1(:,end) = selectNeg*Un(:,end) + selectPos*Un1(:,end);
