function [rho,vel,en] = wignerMoments(L,W)
%WIGNERMOMENTS - Calculate the moments of a given Wigner state W
%
%    [rho,vel,en] = WIGNERMOMENTS(L,W)
%
%    Input:    L        length of the velocity grid
%              W        Wigner state (N x N x 4 array)
%
%    Output:   rho      spin density matrix
%              vel      normalized velocity
%              en       normalized (intrinsic) energy

N = size(W,1);

% spin density matrix
w = zeros(1,4);
for j=1:4
	w(j) = sum(sum(W(:,:,j)))*(2*L/N)^2;
end
rho = pauliToMatrix(w);

% velocity grid
vgrid = [0:N/2-1,-N/2:-1]*(2*L)/N;
[vx,vy] = meshgrid(vgrid,vgrid);

% normalized momentum;
% multiply pointwise by [vx,vy] and average entries;
% factor 2 to compensate for factor 1/2 of trace coefficient in Pauli sigma basis
tx = vx.*W(:,:,1);
ty = vy.*W(:,:,1);
vel = 2*[sum(tx(:));sum(ty(:))]*(2*L/N)^2 / trace(rho);

% normalized intrinsic energy;
% multiply by v^2/2 and average entries;
% factor 2 to compensate for factor 1/2 of trace coefficient in Pauli sigma basis
ten = 0.5*((vx-vel(1)).^2 + (vy-vel(2)).^2).*W(:,:,1);
en = 2*sum(ten(:))*(2*L/N)^2 / trace(rho);
