function quadw = fourierI3(N,L,J,M,R)
%FOURIERI3 - Quadrature formula for the I3 integral
%
%    quadw = FOURIERI3(N,L,J,M,R)
%
%    In the output structure, psiR1, psiR2 and weight encode the
%    delta-function using Carleman's representation in Fourier space
%
%    Reference:
%      Jingwei Hu and Lexing Ying
%      A fast spectral algorithm for the quantum Boltzmann collision operator
%      Commun. Math. Sci. 10, pp 989-999 (2012)
%
%    Note that we flip J and M from Hu-Ying paper

kgrid = [0:N/2-1,-N/2:-1];
[kx,ky] = meshgrid(kgrid,kgrid);
% extended grid
kgrid2 = [0:N-1,-N:-1];
[kx2,ky2] = meshgrid(kgrid2,kgrid2);

% uniform discretization of the interval [0,pi]
theta1 = pi*(0:J-1)/J;
% shift by pi/2
theta2 = theta1 + pi/2;

% Gauss-Legendre quadrature
[rho,w] = lgwt(M,-R,R);

psiR1 = cell(1,J*M);
psiR2 = cell(1,J*M);
weight = zeros(1,J*M);

for j=1:J
	for m=1:M
		% integration weights for radial and angular integral
		weight((j-1)*M + m) = pi/J * w(m);

		% Eq. (2.30)
		psiR1{(j-1)*M + m} = exp(1i*pi/L*rho(m) * (kx *cos(theta1(j))+ky *sin(theta1(j))));
		psiR2{(j-1)*M + m} = 2*R * sinc(pi*R/L  * (kx2*cos(theta2(j))+ky2*sin(theta2(j))));
	end
end

% store in structure
quadw.psiR1  = psiR1;
quadw.psiR2  = psiR2;
quadw.weight = weight;
