function [beta,mu] = calculateFDparams(rho,en)
%CALCULATEFDPARAMS - Calculate Fermi-Dirac parameters
%
%    [beta,mu] = CALCULATEFDPARAMS(rho,en)
%
%    Calculate the beta and mu parameters of the Fermi-Dirac distribution
%    matching a given spin density 'rho' and normalized energy 'en'

% column vector
rho = rho(:);

% starting values
z0 = [1,1];
beta0 = 1;

% use un-normalized bare energy
en = en*sum(rho);

options = optimset('Jacobian','on','TolFun',1e-12,'Display','none');
[x,md,~,exitflag] = lsqnonlin(@(x) mdiff(x(1),x(2),x(3),rho,en),[z0,beta0],[],[],options);
md = sqrt(md);
if md>1e-4
	% try a different starting value
	beta0 = 0.5;
	[x,md,~,exitflag] = lsqnonlin(@(x) mdiff(x(1),x(2),x(3),rho,en),[z0,beta0],[],[],options);
	md = sqrt(md);
end
if exitflag<0
	warning('Least squares fitting failed, exitflag: %i\n',exitflag);
end
beta = x(3);
mu   = log(x(1:2))/beta;

% fprintf('distance: %g\n',md);
if md > 1e-10
	warning('Large residual %g\n',md);
end


%%
% difference of the moments, i.e., distance to 'rho' and 'en'
function [y,J] = mdiff(zUp,zDw,beta,rho_t,en_t)

[rhoUp,drhoUp] = densityFD(zUp,beta);	% density spin-up
[rhoDw,drhoDw] = densityFD(zDw,beta);	% density spin-down
[enUp, denUp]  =  energyFD(zUp,beta);	% energy spin-up
[enDw, denDw]  =  energyFD(zDw,beta);	% energy spin-down

% difference vector
y = [rhoUp;rhoDw;enUp+enDw]-[rho_t;en_t];

% Jacobi matrix
if nargout > 1
	J = [drhoUp(1),       0 ,drhoUp(2);...
		        0 ,drhoDw(1),drhoDw(2);...
		  denUp(1), denDw(1),denUp(2)+denDw(2)];
end


%%
% density of the Fermi-Dirac distribution function in 2D
% for dispersion omega(v) = v^2/2 and fugacity z = exp(beta*mu)
function [rho,drho] = densityFD(z,beta)

rho = 2*pi/beta * log(1+z);

% derivative with respect to z and beta
drho = [2*pi/(beta*(1+z)),-rho/beta];


%%
% energy of the Fermi-Dirac distribution function in 2D
% for dispersion omega(v) = v^2/2 and fugacity z = exp(beta*mu)
function [en,den] = energyFD(z,beta)

en = -2*pi/beta^2 * polylog2(-z);

% derivative with respect to z and beta
den = [densityFD(z,beta)/(beta*z),-2*en/beta];


%%
% use symbolic math toolbox to evaluate the polylogarithm function Li_2(x)
function y = polylog2(x)

y = double(feval(symengine,'dilog',1-x));
