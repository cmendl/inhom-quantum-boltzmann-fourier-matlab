function visualizeInhomBoltzmann(N,L,J,M,R,h,dt,btype,W,titlestr)
%VISUALIZEINHOMBOLTZMANN - Visualize simulation results for the spatially inhomogeneous Boltzmann equation
%
%    VISUALIZEINHOMBOLTZMANN(N,L,J,M,R,h,dt,btype,W,titlestr)
%
%    Input:    N        number of grid points in each direction
%              L        periodized interval is [-L,L] in each direction
%              J        number of angular discretization points
%              M        number of Legendre-Gauss quadrature nodes
%              R        Wigner functions should be supported on ball with radius R
%              h        spatial mesh width (finite volume cell size)
%              dt       time step
%              btype    boundary conditions: 'periodic', 'Dirichlet' or 'Maxwell'
%              W        Wigner state (simulation output); array with dimension N x N x 4 x numVol x numsteps,
%                       represented in physical velocity space
%              titlestr title string

% consistency checks
assert(length(size(W)) == 5);
assert(N == size(W,1));
assert(N == size(W,2));

% time steps
t = dt*(0:size(W,5)-1);

% spatial mesh
x = h*(0:size(W,4)-1);

titlestr = sprintf('N=%d, L=%g, J=%d, M=%d, R=%g, h=%g, dt=%g\n%s',N,L,J,M,R,h,dt,titlestr);

%%

% represent as matrix entries
fprintf('Representing as matrix...');
Wmat = zeros(N,N,2,2,length(x),length(t));
for it=1:length(t)
	for ix=1:length(x)
		% diagonal entries
		Wmat(:,:,1,1,ix,it) = W(:,:,1,ix,it) + W(:,:,4,ix,it);		% w_tr + w_z
		Wmat(:,:,2,2,ix,it) = W(:,:,1,ix,it) - W(:,:,4,ix,it);		% w_tr - w_z
		% off-diagonal entries
		Wmat(:,:,2,1,ix,it) = W(:,:,2,ix,it) + 1i*W(:,:,3,ix,it);	% w_x + i*w_y
		Wmat(:,:,1,2,ix,it) = W(:,:,2,ix,it) - 1i*W(:,:,3,ix,it);	% w_x - i*w_y
	end
end
fprintf(' Done.\n');

% rearrange dimensions for efficiency,
% resulting in a "2 x 2 x N x N x length(x) x length(t)" array
WmatP = permute(Wmat,[3,4,1,2,5,6]);

% eigenvalues
fprintf('Calculating eigenvalues...');
Weig = zeros(N,N,2,length(x),length(t));
for it=1:length(t)
	for ix=1:length(x)
		for k2=1:N
			for k1=1:N
				lambda = eig(WmatP(:,:,k1,k2,ix,it));
				Weig(k1,k2,1,ix,it) = lambda(1);
				Weig(k1,k2,2,ix,it) = lambda(2);
			end
		end
	end
end
fprintf(' Done.\n');

% Fourier representation
fprintf('Calculating Fourier representation...');
WF = zeros(size(W));
for it=1:length(t)
	for ix=1:length(x)
		for j=1:4
			WF(:,:,j,ix,it) = fft2(W(:,:,j,ix,it))/N^2;
		end
	end
end
fprintf(' Done.\n');

% global Wigner state: average over finite volumes
Wglob = zeros(N,N,4,length(t));
for it=1:length(t)
	for j=1:4
		for ix=1:length(x)
			Wglob(:,:,j,it) = Wglob(:,:,j,it) + W(:,:,j,ix,it);
		end
		Wglob(:,:,j,it) = Wglob(:,:,j,it)/length(x);
	end
end


%%
% conservation laws

% only for periodic boundary conditions
if strcmp(btype,'periodic')

	% differences to initial moments
	drho = zeros(1,length(t));
	dvel = zeros(1,length(t));
	deny = zeros(1,length(t));

	% initial moments
	[rho_init,vel_init,en_init] = wignerMoments(L,Wglob(:,:,:,1));
	for it=1:length(t)
		% calculate moments of spatial average
		[rho,vel,en] = wignerMoments(L,Wglob(:,:,:,it));

		% differences
		drho(it) = norm(rho - rho_init)/norm(rho_init);
		dvel(it) = norm(vel - vel_init)/norm(vel_init);
		deny(it) =  abs( en -  en_init)/ abs( en_init);
	end

	% density
	fprintf('initial average spin density:\n');
	disp(rho_init);
	lambda = eig(rho_init);
	fprintf('eigenvalues: (%g, %g)\n',lambda(1),lambda(2));
	figure();
	semilogy(t,drho);
	title(sprintf('spin density matrix conservation\n%s',titlestr));
	xlabel('t');
	ylabel('|\rho(t) - \rho(0)| / |\rho(0)|');

	% velocity
	fprintf('initial average velocity: [%g,%g]\n',vel_init(1),vel_init(2));
	figure();
	semilogy(t,dvel);
	title(sprintf('momentum conservation\n%s',titlestr));
	xlabel('t');
	ylabel('|v(t) - v(0)| / |v(0)|');

	% energy
	fprintf('initial average energy: %g\n',en_init);
	figure();
	semilogy(t,deny);
	title(sprintf('energy conservation\n%s',titlestr));
	xlabel('t');
	ylabel('|\epsilon(t) - \epsilon(0)| / |\epsilon(0)|');
end

%%
% total entropy

S = zeros(1,length(t));
for it=1:length(t)
	tmp = zeros(N);
	for ix=1:length(x)
		% sum of entropies given by eigenvalues
		tmp = tmp + entropy(Weig(:,:,1,ix,it)) + entropy(Weig(:,:,2,ix,it));
	end
	tmp = h*tmp;	% take spatial finite volume size into account
	S(it) = sum(tmp(:)) * (2*L/N)^2;
end
fprintf('initial and final entropy: %g, %g\n',S(1),S(end));
figure();
plot(t,S);
title(sprintf('entropy (spatial integral)\n%s',titlestr));
xlabel('t');
ylabel('S(t)');

%%
% convergence to Fermi-Dirac

% only for periodic boundary conditions
if strcmp(btype,'periodic')

	% use last global Wigner state as reference to construct
	% the asymptotic Fermi-Dirac equilibrium state
	WFDglob = asymptoticFD(L,Wglob(:,:,:,end));

	fprintf('Calculating quantum relative entropy...');
	KLglob = zeros(1,length(t));	% quantum relative entropy with respect to global asymptotic Fermi-Dirac equilibrium state
	KLloc  = zeros(1,length(t));	% quantum relative entropy with respect to local Fermi-Dirac equilibrium states
	KLeqfd = zeros(1,length(t));	% quantum relative entropy between local and global Fermi-Dirac equilibrium states
	for it=1:length(t)
		% fprintf('time step %d/%d:\n',it,length(t));
		for ix=1:length(x)

			% local Fermi-Dirac state corresponding to Wigner state at current spatial location
			WFDloc = asymptoticFD(L,W(:,:,:,ix,it));

			for k2=1:N
				for k1=1:N
					% quantum relative entropy with respect to global Fermi-Dirac state;
					% take spatial finite volume size into account and normalize velocity integral
					if WFDglob(k1,k2,1) > 1e-6		% ignore very small entries since they lead to artificially large values of the relative entropy
						KLglob(it) = KLglob(it) + relentropy(WmatP(:,:,k1,k2,ix,it),...
							pauliToMatrix(squeeze(WFDglob(k1,k2,:)))) * h*(2*L/N)^2;
					end
				
					if WFDloc(k1,k2,1) > 1e-6		% ignore very small entries since they lead to artificially large values of the relative entropy
						KLloc(it) = KLloc(it) + relentropy(WmatP(:,:,k1,k2,ix,it),...
							pauliToMatrix(squeeze(WFDloc(k1,k2,:)))) * h*(2*L/N)^2;
					end
					
					if WFDglob(k1,k2,1) > 1e-6		% ignore very small entries since they lead to artificially large values of the relative entropy
						KLeqfd(it) = KLeqfd(it) + relentropy(pauliToMatrix(squeeze(WFDloc(k1,k2,:))),...
							pauliToMatrix(squeeze(WFDglob(k1,k2,:)))) * h*(2*L/N)^2;
					end
				end
			end
		end
	end
	fprintf(' Done.\n');
	
	fprintf('initial and final global      quantum relative entropy (Kullback-Leibler divergence): %g, %g\n',KLglob(1),KLglob(end));
	fprintf('initial and final local       quantum relative entropy (Kullback-Leibler divergence): %g, %g\n',KLloc(1),KLloc(end));
	fprintf('initial and final equilibrium quantum relative entropy (Kullback-Leibler divergence): %g, %g\n',KLeqfd(1),KLeqfd(end));
	figure();
	semilogy(t,KLglob,'b',t,KLloc,'g',t,KLeqfd,'--r');
	title(sprintf('quantum relative entropy\n%s',titlestr));
	legend('S(W(t) || W_{FD,glob})','S(W(t) || W_{FD,loc}(t))','S(W_{FD,loc}(t) || W_{FD,glob})')
	xlabel('t');
end


%%
% stationary, spatially dependent spin density matrix

figure();
hold on;
clr = {'b','g','r','c'};
rho_glob = zeros(4,length(x));
for j=1:4
	% zero frequency entry in Fourier representation equals mean value
	rho_glob(j,:) = (2*L)^2 * squeeze(WF(1,1,j,:,end)).';
	if strcmp(btype,'periodic')
		plot([x,x(end)+h],[rho_glob(j,:),rho_glob(j,1)],clr{j});	% repeat entry at x = 0
	else
		plot(x,rho_glob(j,:),clr{j});
	end
	title(sprintf('final density\n%s',titlestr));
	xlabel('x');
	ylabel('\rho(t,x)');
end
legend('\rho_0','\rho_1','\rho_2','\rho_3');

%%
% stationary, spatially dependent (local) inverse temperature beta

beta_loc = zeros(1,length(x));
for ix=1:length(x)
	% local Fermi-Dirac state corresponding to Wigner state at current spatial location
	[~,beta_loc(ix)] = asymptoticFD(L,W(:,:,:,ix,end));
end

figure();
if strcmp(btype,'periodic')
	plot([x,x(end)+h],[beta_loc,beta_loc(1)]);	% repeat entry at x = 0
else
	plot(x,beta_loc);
end
title(sprintf('inverse temperature\n%s',titlestr));
xlabel('x');
ylabel('\beta(x)');

% also plot temperature T
T = 1./beta_loc;
figure();
if strcmp(btype,'periodic')
	plot([x,x(end)+h],[T,T(1)]);	% repeat entry at x = 0
else
	plot(x,T);
end
title(sprintf('temperature\n%s',titlestr));
xlabel('x');
ylabel('k_B T(x)');


%%
% stationary, spatially dependent (local) entropy

S_loc = zeros(1,length(x));
for ix=1:length(x)
	% sum of entropies given by eigenvalues
	tmp = entropy(Weig(:,:,1,ix,end)) + entropy(Weig(:,:,2,ix,end)) * (2*L/N)^2;
	S_loc(ix) = sum(tmp(:)) * (2*L/N)^2;
end

figure();
if strcmp(btype,'periodic')
	plot([x,x(end)+h],[S_loc,S_loc(1)]);	% repeat entry at x = 0
else
	plot(x,S_loc);
end
title(sprintf('local entropy\n%s',titlestr));
xlabel('x');
ylabel('S(x)');

%%
% spatially dependent spin density matrix

% rhoX = zeros(2,2,length(x),length(t));
% for it=1:length(t)
% 	for ix=1:length(x)
% 		w = (2*L)^2 * squeeze(WF(1,1,:,ix,it));
% 		rhoX(:,:,ix,it) = pauliToMatrix(w);
% 	end
% end

% animation of spin density matrix
hf = figure();
rho_mean_entries = (2*L)^2 * squeeze(WF(1,1,:,:,:));
ymin = min(rho_mean_entries(:)) - 1;
ymax = max(rho_mean_entries(:)) + 1;
for it=1:length(t)
	clf(hf);
	hold on;
	for j=1:4
		% zero frequency entry in Fourier representation equals mean value
		rho_glob = (2*L)^2 * squeeze(WF(1,1,j,:,it)).';
		if strcmp(btype,'periodic')
			plot([x,x(end)+h],[rho_glob,rho_glob(1)],clr{j});	% repeat entry at x = 0
		else
			plot(x,rho_glob,clr{j});
		end
	end
	ylim([ymin,ymax]);
	title(sprintf('density \\rho, t=%g\n%s',(it-1)*dt,titlestr));
	xlabel('x');
	ylabel('\rho(t,x)');
	pause(0.05);
end
