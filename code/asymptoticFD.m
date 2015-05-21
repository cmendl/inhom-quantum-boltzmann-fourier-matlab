function [WFD,beta,mu,vel,U] = asymptoticFD(L,W)
%ASYMPTOTICFD - Construct asymptotic Fermi-Dirac equilibrium state 
%
%    [WFD,beta,mu,vel,U] = ASYMPTOTICFD(L,W)
%
%    Construct the Fermi-Dirac equilibrium state corresponding to the
%    conserved moments (density, velocity and energy) of the input Wigner state W
%
%    Input:    L		periodic domain [-L,L]^2
%              W		Wigner state (N x N x 4 tensor)
%
%    Output:   WFD		Fermi-Dirac Wigner state with same moments as W
%              beta     inverse temperature
%              mu       chemical potential for spin-up and spin-down
%              vel      average velocity
%              U        eigenbasis

N = size(W,1);

% conserved moments
[rho,vel,en] = wignerMoments(L,W);

% diagonalize density matrix
[U,lambda] = eig(rho);
lambda = diag(lambda);

% reconstruct beta and mu
[beta,mu] = calculateFDparams(lambda,en);

% global Fermi-Dirac state corresponding to conserved momenta
WFD = createFD(N,L,beta,mu,vel,U);
