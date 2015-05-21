function quadw = fourierI1(N,L,J,R)
%FOURIERI1 - Quadrature formula for the I1 integral
%
%    quadw = FOURIERI1(N,L,J,R)
%
%    In the output structure, psiR1, psiR2 and weight encode the
%    1/(u \cdot u') principal value using Carleman's representation

kgrid1 = [0:N/2-1,-N/2:-1];
[kx1,ky1] = meshgrid(kgrid1,kgrid1);
% extended grid
kgrid2 = [0:N-1,-N:-1];
[kx2,ky2] = meshgrid(kgrid2,kgrid2);

% uniform discretization of the interval [0,pi]
theta1 = pi*(0:J-1)/J;
% shift by pi/(2*J), corresponds to midpoint rule
theta2 = theta1 + pi/(2*J);

psiR1 = cell(1,J^2);
psiR2 = cell(1,J^2);
weight = zeros(1,J^2);

for j2=1:J
	for j1=1:J
		j = j1+J*(j2-1);
		% principal value singularity, minus sign due to i^2
		weight(j) = -(2*pi/J)^2/cos(theta1(j1)-theta2(j2));
		% note that entries appear repeatedly
		psiR1{j}  = R * f(pi*R/(2*L) * (kx2*cos(theta1(j1))+ky2*sin(theta1(j1))));
		psiR2{j}  = R * f(pi*R/(2*L) * (kx1*cos(theta2(j2))+ky1*sin(theta2(j2))));
	end
end

% store in structure
quadw.psiR1  = psiR1;
quadw.psiR2  = psiR2;
quadw.weight = weight;


function y = f(x)
%
% sin(x)*sinc(x), same as imaginary part of (exp(2i*x)-1)/(2i*x)

y = sin(x) .* sinc(x);
