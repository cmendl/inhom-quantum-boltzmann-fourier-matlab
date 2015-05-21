function W = createFD(N,L,beta,mu,vel,U)
%CREATEFD - Create Fermi-Dirac equilibrium state
%
%    W = CREATEFD(N,L,beta,mu,vel,U)
%
%    Create Fermi-Dirac equilibrium state (N x N x 4 tensor),
%    with the last dimension corresponding to coefficients
%    of the Pauli sigma basis

% velocity grid
vgrid = [0:N/2-1,-N/2:-1]*(2*L)/N;
[vx,vy] = meshgrid(vgrid,vgrid);

W = zeros(N^2,4);

% row vector
vel = vel(:).';

% evaluate in physical velocity space
for j=1:numel(vx)
	Wloc = U*FDmat(beta,mu,[vx(j),vy(j)]-vel)*U';
	% represent in Pauli sigma basis
	W(j,:) = matrixToPauli(Wloc);
end

W = reshape(W,[N,N,4]);


function f = FD(beta,mu,v)
%
% Fermi-Dirac distribution

f = 1/(exp(beta*(0.5*norm(v)^2 - mu)) + 1);


function f = FDmat(beta,mu,v)
%
% Fermi-Dirac distribution matrix

f = diag([FD(beta,mu(1),v),FD(beta,mu(2),v)]);
