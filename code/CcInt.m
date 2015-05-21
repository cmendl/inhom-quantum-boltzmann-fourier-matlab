function Cc = CcInt(W,quadwI1)
%CCINT - Calculate the integral of the conservative collision operator Cc
%
%    Cc = CCINT(W,quadwI1)
%
%    Input:    W        Wigner state in Fourier space (1 x 4 cell of N x N matrices)
%              quadwI1  quadrature formula for I1 integral
%
%    Output:   Cc       collision state

Cc = cell(1,4);

% no trace terms
Cc{1} = zeros(size(W{1}));

% (w_{3,tr} + w_{3,tr} - 1) i [ W_2, W_1 ] =
% - 2 (2*w_{3,tr} - 1) ( \vec w_2 cross \vec w_1 )\cdot\sigma

% 2*w_{3,tr} - 1
Wtr2m1 = 2*W{1};
% -ones(N) in Fourier space, divided by N^2 due to normalization
Wtr2m1(1,1) = Wtr2m1(1,1) - 1;

crossIX = [...
	3,4;...	% y, z
	4,2;... % z, x
	2,3 ];	% x, y
for j=2:4
	Cc{j} = -2*(...
		  I4integral(Wtr2m1,W{crossIX(j-1,1)},W{crossIX(j-1,2)},quadwI1.psiR1,quadwI1.psiR2,quadwI1.weight) ...
		- I4integral(Wtr2m1,W{crossIX(j-1,2)},W{crossIX(j-1,1)},quadwI1.psiR1,quadwI1.psiR2,quadwI1.weight));
end
