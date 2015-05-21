function [Cd,Cc] = C_test(N,L,J,M,R)
%C_TEST - Test file for the dissipative and conservative collision operators
%
%    [Cd,Cc] = C_TEST()
%    [Cd,Cc] = C_TEST(N,L,J,M,R)

if (nargin < 4)
	% default parameters
	N = 32;
	R = 7.5;
	% periodic volume [-L,L]^2
	L = ceil(0.5*(3+sqrt(2)) * R);
	J = 32;
	M = 32;
end
fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('M: %d\n',M);
fprintf('R: %g\n',R);

% create Fermi-Dirac state and represent in Fourier space
beta = 1.2;
mu = [-2.5,-1.5];
v0 = [0,0];
U  = [cos(pi/5),-1i*sin(pi/5);sin(pi/5),1i*cos(pi/5)];
Wfd = createFD(N,L,beta,mu,v0,U);
for j=1:4
	% divide by N^2 due to normalization convention
	Wfd(:,:,j) = fft2(Wfd(:,:,j))/N^2;
end

vgrid = [0:N/2-1,-N/2:-1]*(2*L)/N;
[vx,vy] = meshgrid(vgrid,vgrid);
% visualize trace component
figure();
mesh(vx,vy,Wfd(:,:,1));
xlabel('\xi_1');
ylabel('\xi_2');
title('W_{tr}(\xi)');
% visualize sigma_x component
figure();
mesh(vx,vy,Wfd(:,:,2));
xlabel('\xi_1');
ylabel('\xi_2');
title('W_{x}(\xi)');

% represent as cell array
Wfd = {Wfd(:,:,1),Wfd(:,:,2),Wfd(:,:,3),Wfd(:,:,4)};

% dissipative collision operator in Fourier representation
fprintf('Calculating dissipative collision operator...\n');
quadwI2 = fourierI2(N,L,J,R);
quadwI3 = fourierI3(N,L,J,M,R);
quadwI4 = fourierI4(N,L,J,R);
Cd = CdInt(Wfd,quadwI2,quadwI3,quadwI4);
fprintf('Cd{j}(1,1): %g, %g, %g, %g (should be zero due to spin conservation)\n',Cd{1}(1,1),Cd{2}(1,1),Cd{3}(1,1),Cd{4}(1,1));

% conservative collision operator
fprintf('Calculating conservative collision operator...\n');
quadwI1 = fourierI1(N,L,J,R);
Cc = CcInt(Wfd,quadwI1);
fprintf('Cc{j}(1,1): %g, %g, %g, %g (should be zero due to spin conservation)\n',Cc{1}(1,1),Cc{2}(1,1),Cc{3}(1,1),Cc{4}(1,1));
