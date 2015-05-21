function FDparams_test()
%FDPARAMS_TEST - Test file for reconstructing moments and thermodynamic parameters from a Fermi-Dirac state
%
%    FDPARAMS_TEST()

N = 32;
R = 7.5;
% periodic volume [-L,L]^2
L = ceil(0.5*(3+sqrt(2))/sqrt(2) * R);

% parameters for Fermi-Dirac state
beta = 1.0;
mu = [-1.5,1.5];
v0 = [2,3];

% theoretical (analytic) values of density and energy
rho_exact = [densityFD(exp(beta*mu(1)),beta),densityFD(exp(beta*mu(2)),beta)];
 en_exact  =  energyFD(exp(beta*mu(1)),beta)+ energyFD(exp(beta*mu(2)),beta);
% normalization
en_exact = en_exact / sum(rho_exact);

% create Fermi-Dirac state
W = createFD(N,L,beta,mu,v0,eye(2));

% calculate moments from W
[rho,vel,en] = wignerMoments(L,W);

% spin density
rho = eig(rho);
fprintf('spin density matrix eigenvalues: [%g,%g]\n',rho(1),rho(2));
fprintf('                          exact: [%g,%g]\n',rho_exact(1),rho_exact(2));
fprintf('                     difference: [%g,%g]\n',rho(1)-rho_exact(1),rho(2)-rho_exact(2));
fprintf('\n');

% momentum
fprintf('average velocity: [%g,%g]\n',vel(1),vel(2));
fprintf('           exact: [%g,%g]\n',v0(1),v0(2));
fprintf('      difference: [%g,%g]\n',vel(1)-v0(1),vel(2)-v0(2));
fprintf('\n');

% energy
fprintf('average energy: %g\n',en);
fprintf('         exact: %g\n',en_exact);
fprintf('    difference: %g\n',en-en_exact);
fprintf('\n');

% reconstruct beta and mu
[beta1,mu1] = calculateFDparams(rho,en);

fprintf('reconstructed beta and mu: %g, (%g,%g)\n',beta1,mu1(1),mu1(2));
fprintf('       actual beta and mu: %g, (%g,%g)\n',beta,mu(1),mu(2));
fprintf('               difference: %g, (%g,%g)\n',beta1-beta,mu1(1)-mu(1),mu1(2)-mu(2));


% density of the Fermi-Dirac distribution function in 2D
% for dispersion omega(v) = v^2/2 and fugacity z = exp(beta*mu)
function [rho,drho] = densityFD(z,beta)

rho = 2*pi/beta * log(1+z);

% derivative with respect to z and beta
drho = [2*pi/(beta*(1+z)),-rho/beta];


% energy of the Fermi-Dirac distribution function in 2D
% for dispersion omega(v) = v^2/2 and fugacity z = exp(beta*mu)
function [en,den] = energyFD(z,beta)

en = -2*pi/beta^2 * polylog2(-z);

% derivative with respect to z and beta
den = [densityFD(z,beta)/(beta*z),-2*en/beta];


% use symbolic math toolbox to evaluate the polylogarithm function Li_2(x)
function y = polylog2(x)

y = double(feval(symengine,'dilog',1-x));
