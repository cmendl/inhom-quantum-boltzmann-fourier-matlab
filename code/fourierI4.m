function quadw = fourierI4(N,L,J,R)
%FOURIERI4 - Quadrature formula for the I4 integral
%
%    quadw = FOURIERI4(N,L,J,R)
%
%    In the output structure, psiR1, psiR2 and weight encode the
%    delta-function using Carleman's representation in Fourier space

kgrid = [0:N/2-1,-N/2:-1];
[kx,ky] = meshgrid(kgrid,kgrid);
% extended grid
kgrid2 = [0:N-1,-N:-1];
[kx2,ky2] = meshgrid(kgrid2,kgrid2);

% uniform discretization of the interval [0,pi]
theta1 = pi*(0:J-1)/J;
% shift by pi/2
theta2 = theta1 + pi/2; 

psiR1 = cell(1,J);
psiR2 = cell(1,J);
weight = zeros(1,J);

for j=1:J
	% integration weights
	weight(j) = pi/J;

	psiR1{j}  = 2*R * sinc(pi*R/L * (kx2*cos(theta1(j))+ky2*sin(theta1(j))));
	psiR2{j}  = 2*R * sinc(pi*R/L * (kx *cos(theta2(j))+ky *sin(theta2(j))));
end

% store in structure
quadw.psiR1  = psiR1;
quadw.psiR2  = psiR2;
quadw.weight = weight;
