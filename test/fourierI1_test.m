function [err,G,G_ref] = fourierI1_test(N,L,J,R)
%FOURIERI1_TEST - Compare quadrature formula for I1 integral with reference implementation
%
%    [err,G,G_ref] = FOURIERI1_TEST()
%    [err,G,G_ref] = FOURIERI1_TEST(N,L,J,R)

% default parameters
if (nargin < 4)
	N = 4;
	L = 17;
	J = 32;
	R = 7.5;
end
fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('R: %g\n',R);

% reference implementation
G_ref = fourierI1_ref(N,L,R);

quadw = fourierI1(N,L,J,R);

G = zeros(size(G_ref));
for j=1:length(quadw.weight)
	G = G + quadw.weight(j) * quadw.psiR1{j}(:) * quadw.psiR2{j}(:).';
end

err = norm(G - G_ref);
