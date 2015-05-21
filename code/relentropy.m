function KL = relentropy(rho,tau)
%RELENTROPY - Quantum relative entropy
%
%    KL = RELENTROPY(rho,tau) calculates the quantum relative entropy between rho and tau

[U,mu] = schur(rho);
[V,nu] = schur(tau);
assert(isequal(mu,diag(diag(mu))));
assert(isequal(nu,diag(diag(nu))));

% eigenvalues on diagonal
mu = diag(mu);
nu = diag(nu);

% clamp to valid range
mu(mu <= 0) = eps;
nu(nu <= 0) = eps;
mu(mu >= 1) = 1-eps;
nu(nu >= 1) = 1-eps;

% represent in eigenbasis of rho
W = U'*V;

rho = diag(mu);

log_rho =   diag(log(mu));
log_tau = W*diag(log(nu))*W';

log_1rho =   diag(log(1-mu));
log_1tau = W*diag(log(1-nu))*W';

% remove imaginary part resulting from numerical roundoff errors
KL = real(trace(rho*(log_rho - log_tau)) + trace((eye(2)-rho)*(log_1rho - log_1tau)));

