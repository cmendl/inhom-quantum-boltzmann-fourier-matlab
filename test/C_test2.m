function [Cd,Cc] = C_test2(N,L,J,M,R)
%C_TEST2 - Test file for the dissipative and conservative collision operators
%
%    [Cd,Cc] = C_TEST2()
%    [Cd,Cc] = C_TEST2(N,L,J,M,R)

if (nargin < 4)
	% default parameters
	N = 32;
	R = 7.5;
	% periodic volume [-L,L]^2
	L = ceil(0.5*(3+sqrt(2)/sqrt(2)) * R);
	J = 32;
	M = N;
end
fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('M: %d\n',M);
fprintf('R: %g\n',R);


%%
% define Wigner states using Gaussians

vgrid = [0:N/2-1,-N/2:-1]*(2*L)/N;
[vx,vy] = meshgrid(vgrid,vgrid);

% define Wigner components in physical velocity space
Wtr = @(u,v) 1/(3*pi)*(1 + u)*exp(-(u - 1/2)^2/2 - 1/2*(v + u)^2);
Wx  = @(u,v) (3*sqrt(6/5))/(55*pi)*(1 - u - v/6)^2*exp(-(2*u + v)^2/12 - v^2/10);
Wy  = @(u,v) 1/(sqrt(2)*4*pi)*exp(-1/8*(u-1)^2 - 1/4*v^2);
Wz  = @(u,v) 9/(56*pi)*(1/3+u)^2 * exp(-1/8*(u-v)^2 - 1/8*(u+2*v)^2);
% evaluate on grid
W = cell(1,4);
for j=1:4
	W{j} = zeros(N);	% N x N matrix
end
for j=1:numel(W{1})
	W{1}(j) = Wtr(vx(j),vy(j));
	W{2}(j) = Wx (vx(j),vy(j));
	W{3}(j) = Wy (vx(j),vy(j));
	W{4}(j) = Wz (vx(j),vy(j));
end

% visualize trace component
vx_shift = circshift(vx,[N/2,N/2]);
vy_shift = circshift(vy,[N/2,N/2]);
figure();
mesh(vx_shift,vy_shift,circshift(W{1},[N/2,N/2]));
xlabel('v_1');
ylabel('v_2');
title('W_{tr}(v)');
% visualize sigma_x component
figure();
mesh(vx_shift,vy_shift,circshift(W{2},[N/2,N/2]));
xlabel('v_1');
ylabel('v_2');
title('W_{x}(v)');

% represent in Fourier space
for j=1:4
	W{j} = fft2(W{j})/N^2;
	% fprintf('W{%d}(1:5,1:5):\n',j);
	% disp(W{j}(1:5,1:5));
end


%%

% dissipative collision operator
quadwI2 = fourierI2(N,L,J,R);
quadwI3 = fourierI3(N,L,J,M,R);
quadwI4 = fourierI4(N,L,J,R);
Cd = CdInt(W,quadwI2,quadwI3,quadwI4);

% conservative collision operator
quadwI1 = fourierI1(N,L,J,R);
Cc = CcInt(W,quadwI1);

% check spin conservation
fprintf('abs(Cd{j}(1,1)): %g, %g, %g, %g (should be zero due to spin conservation)\n',...
	abs(Cd{1}(1,1)),abs(Cd{2}(1,1)),abs(Cd{3}(1,1)),abs(Cd{4}(1,1)));
fprintf('abs(Cc{j}(1,1)): %g, %g, %g, %g (should be zero due to spin conservation)\n',...
	abs(Cc{1}(1,1)),abs(Cc{2}(1,1)),abs(Cc{3}(1,1)),abs(Cc{4}(1,1)));
fprintf('\n');

%%
% transform back to physical velocity space
fprintf('physical velocity space:\n');
for j=1:4
	Cd{j} = N^2*ifft2(Cd{j});
	Cc{j} = N^2*ifft2(Cc{j});
	fprintf('norm(imag(Cd{%d})) / norm(real(Cd{%d})): %g (should be zero)\n',...
		j,j,norm(imag(Cd{j})) / (norm(real(Cd{j})) + eps));
	fprintf('norm(imag(Cc{%d})) / norm(real(Cc{%d})): %g (should be zero)\n',...
		j,j,norm(imag(Cc{j})) / (norm(real(Cc{j})) + eps));
	Cd{j} = real(Cd{j});
	Cc{j} = real(Cc{j});
end
fprintf('\n');

%%
% test momentum conservation

% multiply by [k1,k2] and sum entries
kgrid = [0:N/2-1,-N/2:-1];
[k1,k2] = meshgrid(kgrid,kgrid);
int_v1_trCd = mean(mean(k1 .* Cd{1}));
int_v2_trCd = mean(mean(k2 .* Cd{1}));
fprintf('\\int v tr[Cd(v)] dv: (%g, %g) (should be zero due to momentum conservation)\n',int_v1_trCd,int_v2_trCd);

%%
% test energy conservation

% multiply by v^2/2 and sum entries
kgrid = [0:N/2-1,-N/2:-1];
[k1,k2] = meshgrid(kgrid,kgrid);
int_vsq_trCd = mean(mean(0.5*(k1.^2 + k2.^2) .* Cd{1}));
int_vsq_trCc = mean(mean(0.5*(k1.^2 + k2.^2) .* Cc{1}));
fprintf('\\int v^2/2 tr[Cd(v)] dv: %g (should be zero due to energy conservation)\n',int_vsq_trCd);
fprintf('\\int v^2/2 tr[Cc(v)] dv: %g (should be zero due to energy conservation)\n',int_vsq_trCc);
