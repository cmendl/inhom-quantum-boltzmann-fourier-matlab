function w = matrixToPauli(W)
%MATRIXTOPAULI - Represent matrix W in Pauli sigma basis
%
%    w = MATRIXTOPAULI(W)

sigma_x = [0, 1 ;1 , 0];
sigma_y = [0,-1i;1i, 0];
sigma_z = [1, 0 ;0 ,-1];

w = 0.5*real([trace(W),trace(W*sigma_x),trace(W*sigma_y),trace(W*sigma_z)]);
