function [err,H2,H3] = fourierI23_test(N,L,J,M,R)
%FOURIERI23_TEST - Compare quadrature formulas for I2 and I3 integrals with each other
%
%    [err,H2,H3] = FOURIERI23_TEST(N,L,J,M,R)

% default parameters
if (nargin < 5)
	N = 8;
	L = 17;
	J = 32;
	M = 32;
	R = 7.5;
end
fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('M: %d\n',M);
fprintf('R: %g\n',R);


%%
% I2 integral

quadwI2 = fourierI2(N,L,J,R);

H2 = zeros(N^2);
for j=1:length(quadwI2.weight)
	H2 = H2 + quadwI2.weight(j) * quadwI2.psiR1{j}(:) * quadwI2.psiR2{j}(:).';
end


%%
% I3 integral

quadwI3 = fourierI3(N,L,J,M,R);

IX = [1:N/2,3*N/2+1:2*N];

H3 = zeros(N^2);
for j=1:length(quadwI3.weight)
	psiR2_mod = quadwI3.psiR2{j}(IX,IX);	% extract [-N/2:N/2-1]^2 part
	H3 = H3 + quadwI3.weight(j) * quadwI3.psiR1{j}(:) * psiR2_mod(:).';
end


%%
% compare

err = norm(H2 - H3);
