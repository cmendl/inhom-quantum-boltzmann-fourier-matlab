function dW = magneticRot(W,B)
%MAGNETICROT - Calculate time derivative of Wigner state W effected by external magnetic field B
%
%    dW = MAGNETICROT(W,B)
%
%    W is a 1 x 4 cell of N x N matrices, representing the Fourier-transformed
%    Wigner matrices in the Pauli-spin basis (including identity as first component)
%
%    B is the external magnetic field stored as 3D vector

dW = cell(1,4);

% no trace terms
dW{1} = zeros(size(W{1}));

% -i [ B, W ] = 2 ( B cross \vec w_1 )\cdot\sigma

% y, z
dW{2} = 2*(B(2)*W{3+1} - B(3)*W{2+1});
% z, x
dW{3} = 2*(B(3)*W{1+1} - B(1)*W{3+1});
% x, y
dW{4} = 2*(B(1)*W{2+1} - B(2)*W{1+1});
