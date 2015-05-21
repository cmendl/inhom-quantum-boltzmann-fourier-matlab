function W = pauliToMatrix(w)
%PAULITOMATRIX - Convert representation in Pauli sigma basis to matrix form
%
%    W = PAULITOMATRIX(w)

sigma_x = [0, 1 ;1 , 0];
sigma_y = [0,-1i;1i, 0];
sigma_z = [1, 0 ;0 ,-1];

W = w(1)*eye(2) + w(2)*sigma_x + w(3)*sigma_y + w(4)*sigma_z;
