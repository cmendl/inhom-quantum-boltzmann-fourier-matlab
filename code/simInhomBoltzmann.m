function W = simInhomBoltzmann()
%SIMINHOMBOLTZMANN - Main file for the simulation of the spatially inhomogeneous quantum Boltzmann equation
%
%    W = SIMINHOMBOLTZMANN()

%%
% simulation parameters

N = 32;
R = 7.5;
% periodic volume [-L,L]^2
L = ceil(0.5*(3+sqrt(2))/sqrt(2) * R);
J = 32;
M = N;

% spatial mesh width
h = 0.1;
xmax = 1;
% x = 0:h:(xmax-h);	% for periodic boundary conditions
x = 0:h:xmax;

% time steps
dt = h/L;
dt = floor(1000*dt)/1000;	% round
assert(dt > 0);
% tmax = 0.1;
tmax = 4*dt;
t = 0:dt:tmax;

% boundary condition
btype = 'Dirichlet';
% btype = 'Maxwell';
% alpha = 0.4;	% accommodation coefficient for Maxwell reflection operator

% Courant number must be <= 1 according to CFL condition
courant = L*dt/h;
fprintf('Courant number: %g (should be <= 1)\n',courant);

% external magnetic field, assumed to be constant
% Bext = [0,0,0.2];
Bext = [0,0,0];
% Bext = [0.25,-0.8,0.1];

fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('M: %d\n',M);
fprintf('R: %g\n',R);
fprintf('h: %g\n',h);
fprintf('dt: %g\n',dt);
fprintf('number of time steps: %d\n',length(t));
fprintf('number of finite volumes: %d\n',length(x));
fprintf('boundary condition: %s\n',btype);
if strcmpi(btype,'Maxwell')
	fprintf('alpha: %g\n',alpha);
end

%%
% common functions and states

% velocity grid
vgrid = [0:N/2-1,-N/2:-1]*(2*L)/N;
[vx,vy] = meshgrid(vgrid,vgrid);

% Pauli sigma matrices
sigma_x = [0, 1 ;1 , 0];
sigma_y = [0,-1i;1i, 0];
sigma_z = [1, 0 ;0 ,-1];

% Fermi-Dirac distribution
FD = @(beta,mu,v) 1/(exp(beta*(0.5*norm(v)^2 - mu)) + 1);

% Fermi-Dirac distribution matrix
FDmat = @(beta,mu,v) diag([FD(beta,mu(1),v),FD(beta,mu(2),v)]);

%%
% define left and right Wigner states for Maxwell reflection operator

if strcmpi(btype,'Maxwell')
	WMaxwL = cell(1,4);
	WMaxwR = cell(1,4);

	for j=1:4
		WMaxwL{j} = zeros(N);
		WMaxwR{j} = zeros(N);
	end

	for k=1:numel(vx)
		% parameters for left Fermi-Dirac state
		beta = 0.6;
		mu = [-0.4,1.8];
		Wfd = FDmat(beta,mu,[vx(k),vy(k)]);
		% represent in Pauli sigma basis
		w = matrixToPauli(Wfd);
		for j=1:4
			WMaxwL{j}(k) = w(j);
		end
		% parameters for right Fermi-Dirac state
		beta = 1.1;
		mu = [1.3,-0.9];
		Wfd = FDmat(beta,mu,[vx(k),vy(k)]);
		% represent in Pauli sigma basis
		w = matrixToPauli(Wfd);
		for j=1:4
			WMaxwR{j}(k) = w(j);
		end
	end
end

%%

% Wigner states for all time steps
W = cell(length(t),length(x),4);
% initialize data structure
for it=1:length(t)
	for ix=1:length(x)
		for j=1:4
			W{it,ix,j} = zeros(N);
		end
	end
end

% define initial Wigner state

% perturb Fermi-Dirac state by k-dependent rotation
for ix=1:length(x)
	for k=1:numel(vx)
		% parameters for local Fermi-Dirac state
		beta = 0.8 + 0.5*sin(2*pi*x(ix));
		mu = [1,1.5+2*sin(vx(k)+vy(k) + cos(4*pi*x(ix)))];
		v0 = [0.4-0.5*erf(vx(k)),-0.1];
		% rotate
		X = sin(vy(k)^2)*eye(2) + (vx(k)-1)*sigma_x + 2/(vx(k)^2 + abs(vy(k)) + 1)*sigma_y + cos(vx(k)+vy(k)/2)*sigma_z;
		U = expm(1i*X);
		Wfd = U*FDmat(beta,mu,[vx(k),vy(k)]-v0)*U';
		% represent in Pauli sigma basis
		w = matrixToPauli(Wfd);
		for j=1:4
			W{1,ix,j}(k) = w(j);
		end
	end
end

%%
% calculate quadrature weights
quadwI1 = fourierI1(N,L,J,R);
quadwI2 = fourierI2(N,L,J,R);
quadwI3 = fourierI3(N,L,J,M,R);
quadwI4 = fourierI4(N,L,J,R);

%%
% run actual simulation

if strcmpi(btype,'Maxwell')
	% indices of velocities to be used for Maxwell boundary condition,
	% omitting velocity zero and repeating entry N/2+1 (corresponding to velocity -L)
	indMaxw = [2:N/2+1,N/2+1:N];
	% grid for transport with Maxwell boundary condition
	vxMaxw = vx(:,indMaxw);
	vxMaxw(:,N/2) = -vxMaxw(:,N/2);		% flip -L -> +L
	% extract velocity components for Maxwell boundary condition
	for j=1:4
		WMaxwL{j} = WMaxwL{j}(:,indMaxw);
		WMaxwR{j} = WMaxwR{j}(:,indMaxw);
	end
	% normalization
	normL =  max(vxMaxw(:).',0) * WMaxwL{1}(:);		% rightward flux at left boundary
	normR = -min(vxMaxw(:).',0) * WMaxwR{1}(:);		% leftward flux at right boundary
	for j=1:4
		WMaxwL{j} = WMaxwL{j} / normL;
		WMaxwR{j} = WMaxwR{j} / normR;
	end
end

for it=2:length(t)

	fprintf('t = %g\n',t(it));

	% temporary variable
	Wtmp = cell(length(x),4);
	for ix=1:length(x)
		for j=1:4
			Wtmp{ix,j} = zeros(N);
		end
	end

	% half step of transport term (Trotter splitting)
	for j=1:4
		if strcmpi(btype,'Maxwell')
			Wx = zeros(N,N,length(x));
			% copy x values into array
			for ix=1:length(x)
				Wx(:,:,ix) = W{it-1,ix,j}(:,indMaxw);
			end
			% sort velocity components into column vectors
			Wx = reshape(Wx,[N^2,length(x)]);
			% outgoing fluxes based on trace component
			if j==1
				fluxL = -min(vxMaxw(:).',0) * Wx(:,1);		% leftward  flux at left  boundary
				fluxR =  max(vxMaxw(:).',0) * Wx(:,end);	% rightward flux at right boundary
			end
			% transport for half a time step
			Wx = slopeLimStepMaxwell(h,0.5*dt,diag(vxMaxw(:)),alpha,fluxL,fluxR,WMaxwL{j}(:),WMaxwR{j}(:),Wx);
			Wx = reshape(Wx,[N,N,length(x)]);
			% copy values to temporary variable
			for ix=1:length(x)
				Wtmp{ix,j}(:,indMaxw) = Wx(:,:,ix);
				% explicitly set -L component since N/2+1 appears twice in 'indMaxw'
				Wtmp{ix,j}(:,N/2+1) = Wx(:,N/2+1,ix);
				% components with velocity zero in x direction not changed by transport
				Wtmp{ix,j}(:,1) = W{it-1,ix,j}(:,1);
			end
		else
			for k=1:numel(vx)	% for each velocity component...
				% copy x values into array
				Wx = zeros(size(x));
				for ix=1:length(x)
					Wx(ix) = W{it-1,ix,j}(k);
				end
				% transport with velocity vx(k) for half a time step
				if strcmpi(btype,'periodic')
					Wx = slopeLimStep(h,0.5*dt,vx(k),Wx);
				elseif strcmpi(btype,'Dirichlet')
					Wx = slopeLimStepDirichlet(h,0.5*dt,vx(k),Wx);
				end
				% copy values to temporary variable
				for ix=1:length(x)
					Wtmp{ix,j}(k) = Wx(ix);
				end
			end
		end
	end

	effVol = 1:length(x);
	% for Dirichlet boundary conditions, omit boundary volumes from collision step
	if strcmpi(btype,'Dirichlet')
		effVol = effVol(2:end-1);
	end

	% represent in Fourier space
	for ix=effVol
		for j=1:4
			% divide by N^2 due to normalization convention
			Wtmp{ix,j} = fft2(Wtmp{ix,j})/N^2;
		end
	end

	% local collisions and magnetic rotation at each x
	for ix=effVol
		% midpoint rule
		Wloc = Wtmp(ix,:);
		Cc = CcInt(Wloc,quadwI1);
		Cd = CdInt(Wloc,quadwI2,quadwI3,quadwI4);
		CB = magneticRot(Wloc,Bext);
		% temporary variable for midpoint rule
		Wast = cell(1,4);
		for j=1:4
			Wast{j} = Wloc{j} + dt/2 * (Cc{j} + Cd{j} + CB{j});
		end

		Cc = CcInt(Wast,quadwI1);
		Cd = CdInt(Wast,quadwI2,quadwI3,quadwI4);
		CB = magneticRot(Wast,Bext);
		for j=1:4
			Wtmp{ix,j} = Wtmp{ix,j} + dt * (Cc{j} + Cd{j} + CB{j});
		end
	end

	% transform to physical velocity space
	for ix=effVol
		for j=1:4
			% multiply by N^2 due to normalization convention
			Wtmp{ix,j} = N^2*real(ifft2(Wtmp{ix,j}));
		end
	end

	% half step of transport term (Trotter splitting)
	for j=1:4
		if strcmpi(btype,'Maxwell')
			Wx = zeros(N,N,length(x));
			% copy x values into array
			for ix=1:length(x)
				Wx(:,:,ix) = Wtmp{ix,j}(:,indMaxw);
			end
			% sort velocity components into column vectors
			Wx = reshape(Wx,[N^2,length(x)]);
			% outgoing fluxes based on trace component
			if j==1
				fluxL = -min(vxMaxw(:).',0) * Wx(:,1);		% leftward  flux at left  boundary
				fluxR =  max(vxMaxw(:).',0) * Wx(:,end);	% rightward flux at right boundary
			end
			% transport for half a time step
			Wx = slopeLimStepMaxwell(h,0.5*dt,diag(vxMaxw(:)),alpha,fluxL,fluxR,WMaxwL{j}(:),WMaxwR{j}(:),Wx);
			Wx = reshape(Wx,[N,N,length(x)]);
			% copy values back
			for ix=1:length(x)
				W{it,ix,j}(:,indMaxw) = Wx(:,:,ix);
				% explicitly set -L component since N/2+1 appears twice in 'indMaxw'
				W{it,ix,j}(:,N/2+1) = Wx(:,N/2+1,ix);
				% components with velocity zero in x direction not changed by transport
				W{it,ix,j}(:,1) = Wtmp{ix,j}(:,1);
			end
		else
			for k=1:numel(vx)	% for each velocity component...
				% copy x values into array
				Wx = zeros(size(x));
				for ix=1:length(x)
					Wx(ix) = Wtmp{ix,j}(k);
				end
				% transport with velocity vx(k) for half a time step
				if strcmpi(btype,'periodic')
					Wx = slopeLimStep(h,0.5*dt,vx(k),Wx);
				elseif strcmpi(btype,'Dirichlet')
					Wx = slopeLimStepDirichlet(h,0.5*dt,vx(k),Wx);
				end
				% copy values back
				for ix=1:length(x)
					W{it,ix,j}(k) = Wx(ix);
				end
			end
		end
	end
end

%%
% save final results to disk
fname = ['../output/simdata_inhom_N',num2str(N),'_',btype,'_',strrep(strrep(datestr(clock()),' ','_'),':','_'),'.mat'];
if strcmpi(btype,'Maxwell')
	save(fname ,'N','L','R','J','M','W','h','x','dt','t','Bext','btype','alpha');
else
	save(fname ,'N','L','R','J','M','W','h','x','dt','t','Bext','btype');
end
fprintf('results saved to file "%s"\n',fname);
