function Cd = CdInt(W,quadwI2,quadwI3,quadwI4)
%CDINT - Calculate the integral of the dissipative collision operator Cd
%
%    Cd = CDINT(W,quadwI2,quadwI3,quadwI4)
%
%    Input:    W        Wigner state in Fourier space (1 x 4 cell of N x N matrices)
%              quadwI2  quadrature formula for I2 integral
%              quadwI3  quadrature formula for I3 integral
%              quadwI4  quadrature formula for I4 integral
%
%    Output:   Cd       collision state

Cd = cell(1,4);

eta = [1,-1,-1,-1];


%%
% 2 <w_3, w_4>_{\eta} (id - W_1)
Wt = tildeWigner(W);
for j=1:4
	Cd{j} = 2*(...
		I2integral(W{1},eta(1)*W{1},Wt{j},quadwI2.psiR1,quadwI2.psiR2,quadwI2.weight) +...
		I2integral(W{2},eta(2)*W{2},Wt{j},quadwI2.psiR1,quadwI2.psiR2,quadwI2.weight) +...
		I2integral(W{3},eta(3)*W{3},Wt{j},quadwI2.psiR1,quadwI2.psiR2,quadwI2.weight) +...
		I2integral(W{4},eta(4)*W{4},Wt{j},quadwI2.psiR1,quadwI2.psiR2,quadwI2.weight));
end

%%
% - 2 <w_3, w_4>_{\eta} (\eta w_2)\cdot\sigma
for j=1:4
	Cd{j} = Cd{j} - 2*(...
		I3integral(W{1},eta(1)*W{1},eta(j)*W{j},quadwI3.psiR1,quadwI3.psiR2,quadwI3.weight) +...
		I3integral(W{2},eta(2)*W{2},eta(j)*W{j},quadwI3.psiR1,quadwI3.psiR2,quadwI3.weight) +...
		I3integral(W{3},eta(3)*W{3},eta(j)*W{j},quadwI3.psiR1,quadwI3.psiR2,quadwI3.weight) +...
		I3integral(W{4},eta(4)*W{4},eta(j)*W{j},quadwI3.psiR1,quadwI3.psiR2,quadwI3.weight));
end

%%
% (w_{3,tr} + w_{3,tr} - 1) anticomm( (\eta w_2)\cdot\sigma, W_1 ) =
% 2 (2*w_{3,tr} - 1) ( (w_{1,tr}(eta w_2) + w_{2,tr} w_1)\cdot\sigma - <w_1, w_2>_{id} id)

Wtr2m1 = 2*W{1};
% -ones(N) in Fourier space, divided by N^2 due to normalization
Wtr2m1(1,1) = Wtr2m1(1,1) - 1;

% sub-expression 2 (2*w_{3,tr} - 1) (w_{1,tr}(eta w_2))\cdot\sigma
for j=1:4
	Cd{j} = Cd{j} + 2*I4integral(Wtr2m1,eta(j)*W{j},W{1},quadwI4.psiR1,quadwI4.psiR2,quadwI4.weight);
end

% sub-expression 2 (2*w_{3,tr} - 1) (w_{2,tr} w_1)\cdot\sigma
for j=1:4
	Cd{j} = Cd{j} + 2*I4integral(Wtr2m1,W{1},W{j},quadwI4.psiR1,quadwI4.psiR2,quadwI4.weight);
end

% sub-expression - 2 (2*w_{3,tr} - 1) <w_1, w_2>_{id} id
for j=1:4
	Cd{1} = Cd{1} - 2*I4integral(Wtr2m1,W{j},W{j},quadwI4.psiR1,quadwI4.psiR2,quadwI4.weight);
end




%%
function Wt = tildeWigner(W)
% implements \tilde{W} = id - W

Wt = cell(1,4);
% FFT of ones(N)
Wt{1} = zeros(size(W{1}));
Wt{1}(1,1) = 1;		% 1 instead of N^2 due to normalization
Wt{1} = Wt{1} - W{1};
for j=2:4
	Wt{j} = -W{j};
end
