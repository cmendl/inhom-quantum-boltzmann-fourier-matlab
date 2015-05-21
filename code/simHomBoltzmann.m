function W = simHomBoltzmann()
%SIMHOMBOLTZMANN - Main file for the simulation of the spatially homogeneous quantum Boltzmann equation
%
%    W = SIMHOMBOLTZMANN()

%%
% simulation parameters

N = 32;
R = 7.5;
% periodic volume [-L,L]^2
L = ceil(0.5*(3+sqrt(2))/sqrt(2) * R);
J = 32;
M = N;
fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('M: %d\n',M);
fprintf('R: %g\n',R);

% time steps
dt = 0.001;
tmax = 0.1;
t = 0:dt:tmax;


%%
% define initial Wigner state

% parameters for Fermi-Dirac state
beta = 0.8;
mu = [1,1.5];
v0 = [0.4,-0.1];
U  = [cos(pi/5),-1i*sin(pi/5);sin(pi/5),1i*cos(pi/5)];
Winit = createFD(N,L,beta,mu,v0,U);
Winit = reshape(Winit,[N^2,4]);

% Pauli sigma matrices
sigma_x = [0, 1 ;1 , 0];
sigma_y = [0,-1i;1i, 0];
sigma_z = [1, 0 ;0 ,-1];

% perturb Fermi-Dirac state by k-dependent rotation
vgrid = [0:N/2-1,-N/2:-1]*(2*L)/N;
[vx,vy] = meshgrid(vgrid,vgrid);
for j=1:numel(vx)
	w = zeros(1,4);
	for m=1:4
		w(m) = Winit(j,m);
	end
	Wloc = pauliToMatrix(w);
	% rotate
	X = sin(vy(j)^2)*eye(2) + (vx(j)-1)*sigma_x + 2/(vx(j)^2 + abs(vy(j)) + 1)*sigma_y + cos(vx(j)+vy(j)/2)*sigma_z;
	U = expm(1i*X);
	Wloc = U'*Wloc*U;
	% store in Pauli sigma representation
	w = matrixToPauli(Wloc);
	for m=1:4
		Winit(j,m) = w(m);
	end
end

Winit = reshape(Winit,[N,N,4]);

% Wigner states for all time steps
W = cell(length(t),4);
% copy initial state
for j=1:4
	W{1,j} = Winit(:,:,j);
end

%%
% calculate quadrature weights
quadwI1 = fourierI1(N,L,J,R);
quadwI2 = fourierI2(N,L,J,R);
quadwI3 = fourierI3(N,L,J,M,R);
quadwI4 = fourierI4(N,L,J,R);

%%
% run actual simulation
Wtmp = cell(1,4);
Wast = cell(1,4);
for it=2:length(t)

	fprintf('t = %g\n',t(it));

	% represent in Fourier space
	for j=1:4
		% divide by N^2 due to normalization convention
		Wtmp{j} = fft2(W{it-1,j})/N^2;
	end

	% midpoint rule
	Cc = CcInt(Wtmp,quadwI1);
	Cd = CdInt(Wtmp,quadwI2,quadwI3,quadwI4);
	for j=1:4
		Wast{j} = Wtmp{j} + dt/2 * (Cc{j} + Cd{j});
	end

	Cc = CcInt(Wast,quadwI1);
	Cd = CdInt(Wast,quadwI2,quadwI3,quadwI4);
	for j=1:4
		Wtmp{j} = Wtmp{j} + dt * (Cc{j} + Cd{j});
	end

	% transform back to physical velocity space
	for j=1:4
		% multiply by N^2 due to normalization convention
		W{it,j} = N^2*real(ifft2(Wtmp{j}));
	end
end

%%
% save final results to disk
fname = ['../output/simdata_hom_N',num2str(N),'_',strrep(strrep(datestr(clock()),' ','_'),':','_'),'.mat'];
save(fname ,'N','L','R','J','M','W','dt','t');
fprintf('results saved to file "%s"\n',fname);
