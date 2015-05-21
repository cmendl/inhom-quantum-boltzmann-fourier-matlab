function visualizeHomBoltzmann(N,L,J,M,R,dt,W)
%VISUALIZEHOMBOLTZMANN - Visualize simulation results for the spatially homogeneous Boltzmann equation
%
%    VISUALIZEHOMBOLTZMANN(N,L,J,M,R,dt,W)
%
%    Input:    N        number of grid points in each direction
%              L        periodized interval is [-L,L] in each direction
%              J        number of angular discretization points
%              M        number of Legendre-Gauss quadrature nodes
%              R        Wigner functions should be supported on ball with radius R
%              dt       time step
%              W        Wigner state (simulation output); array with dimension N x N x 4 x numsteps,
%                       represented in physical velocity space

% consistency checks
assert(length(size(W)) == 4);
assert(N == size(W,1));
assert(N == size(W,2));

% time steps
t = dt*(0:size(W,4)-1);

%%

% represent as matrix entries
Wmat = zeros(N,N,2,2,length(t));
for it=1:length(t)
	% diagonal entries
	Wmat(:,:,1,1,it) = W(:,:,1,it) + W(:,:,4,it);		% w_tr + w_z
	Wmat(:,:,2,2,it) = W(:,:,1,it) - W(:,:,4,it);		% w_tr - w_z
	% off-diagonal entries
	Wmat(:,:,2,1,it) = W(:,:,2,it) + 1i*W(:,:,3,it);	% w_x + i*w_y
	Wmat(:,:,1,2,it) = W(:,:,2,it) - 1i*W(:,:,3,it);	% w_x - i*w_y
end

% eigenvalues
Weig = zeros(N,N,2,length(t));
for it=1:length(t)
	for k2=1:N
		for k1=1:N
			lambda = eig(squeeze(Wmat(k1,k2,:,:,it)));
			Weig(k1,k2,1,it) = lambda(1);
			Weig(k1,k2,2,it) = lambda(2);
		end
	end
end


%%
% conservation laws

% differences to initial moments
drho = zeros(1,length(t));
dvel = zeros(1,length(t));
deny = zeros(1,length(t));

% initial moments
[rho_init,vel_init,en_init] = wignerMoments(L,W(:,:,:,1));
for it=1:length(t)
	% moments at current time
	[rho,vel,en] = wignerMoments(L,W(:,:,:,it));

	% differences
	drho(it) = norm(rho - rho_init)/norm(rho_init);
	dvel(it) = norm(vel - vel_init)/norm(vel_init);
	deny(it) =  abs( en -  en_init)/ abs( en_init);
end

% spin density conservation
fprintf('initial average spin density:\n');
disp(rho_init);
lambda = eig(rho_init);
fprintf('eigenvalues: (%g, %g)\n',lambda(1),lambda(2));
figure();
semilogy(t,drho);
title(sprintf('spin density matrix conservation\nN=%d, L=%g, J=%d, M=%d, R=%g',N,L,J,M,R));
xlabel('t');
ylabel('|\rho(t) - \rho(0)| / |\rho(0)|');

% momentum conservation
fprintf('initial average velocity: [%g,%g]\n',vel_init(1),vel_init(2));
figure();
semilogy(t,dvel);
title(sprintf('momentum conservation\nN=%d, L=%g, J=%d, M=%d, R=%g',N,L,J,M,R));
xlabel('t');
ylabel('|v(t) - v(0)| / |v(0)|');

% energy conservation
fprintf('initial average energy: %g\n',en_init);
figure();
semilogy(t,deny);
title(sprintf('energy conservation\nN=%d, L=%g, J=%d, M=%d, R=%g',N,L,J,M,R));
xlabel('t');
ylabel('|\epsilon(t) - \epsilon(0)| / |\epsilon(0)|');


%%
% entropy

S = zeros(1,length(t));
for it=1:length(t)
	% sum of entropies given by eigenvalues
	tmp = entropy(Weig(:,:,1,it)) + entropy(Weig(:,:,2,it));
	S(it) = sum(tmp(:)) * (2*L/N)^2;
end
fprintf('initial and final entropy: %g, %g\n',S(1),S(end));
figure();
plot(t,S);
title(sprintf('entropy\nN=%d, L=%g, J=%d, M=%d, R=%g',N,L,J,M,R));
xlabel('t');
ylabel('S(t)');

%%
% convergence to Fermi-Dirac (quantum relative entropy)

% use last Wigner state as reference to construct
% the asymptotic Fermi-Dirac equilibrium state
WFD = asymptoticFD(L,W(:,:,:,end));

% L^1-norm distance to asymptotic Fermi-Dirac equilibrium state
distFD = zeros(1,length(t));
for it=1:length(t)
	for j=1:4
		distFD(it) = distFD(it) + norm(W(:,:,j,it)-WFD(:,:,j),1) * (2*L/N)^2;
	end
end
fprintf('initial and final distance to asymptotic Fermi-Dirac equilibrium state: %g, %g\n',distFD(1),distFD(end));

figure();
semilogy(t,distFD);
% title(sprintf('distance to asymptotic Fermi-Dirac equilibrium state\n%s',titlestr));
title('distance to Fermi-Dirac equilibrium state');
xlabel('t');
ylabel('||W(t) - W_{FD}||');


% quantum relative entropy with respect to asymptotic Fermi-Dirac equilibrium state
KL = zeros(1,length(t));
for it=1:length(t)
	for k2=1:N
		for k1=1:N
			if WFD(k1,k2,1) > 1e-6		% ignore very small entries since they lead to artificially large values of the relative entropy
				KL(it) = KL(it) + relentropy(pauliToMatrix(squeeze(W(k1,k2,:,it))),...
					pauliToMatrix(squeeze(WFD(k1,k2,:)))) * (2*L/N)^2;
			end
		end
	end
end
fprintf('initial and final quantum relative entropy (Kullback-Leibler divergence): %g, %g\n',KL(1),KL(end));

figure();
semilogy(t,KL);
title('quantum relative entropy');
xlabel('t');
ylabel('S(W(t) || W_{FD})');

%%
% eigenvalue range

WeigMax = max(reshape(Weig,[N*N*2,length(t)]),[],1);
WeigMin = min(reshape(Weig,[N*N*2,length(t)]),[],1);

figure();
fill([t,t(end:-1:1)],[WeigMax,WeigMin(end:-1:1)],0.9*[1,1,1],'EdgeColor','none');
hold on;
plot(t,WeigMax,'k-');
plot(t,WeigMin,'r--');
ylim([-0.1,1]);
title(sprintf('eigenvalue range\nN=%d, L=%g, J=%d, M=%d, R=%g',N,L,J,M,R));
xlabel('t');
ylabel('eig(W(t,v))');


%%
% visualize final state

% velocity grid
vgrid = [0:N/2-1,-N/2:-1]*(2*L)/N;
[vx,vy] = meshgrid(vgrid,vgrid);

vx_shift = circshift(vx,[N/2,N/2]);
vy_shift = circshift(vy,[N/2,N/2]);

figure();
mesh(vx_shift,vy_shift,circshift(W(:,:,1,end),[N/2,N/2]));
xlabel('v_1');
ylabel('v_2');
title('W_{tr}^{final}(v)');
for j=1:3
	figure();
	mesh(vx_shift,vy_shift,circshift(W(:,:,j+1,end),[N/2,N/2]));
	xlabel('v_1');
	ylabel('v_2');
	title(sprintf('W_{%c}^{final}(v)','x'+j-1));
end


%%
% animation

figure();
for it=1:length(t)
	% mesh(vx_shift,vy_shift,circshift(W(:,:,1,it),[N/2,N/2]));
	% xlabel('v_1');
	% ylabel('v_2');
	% title(sprintf('W_{tr}(v), t=%g',(it-1)*dt));
	mesh(vx_shift,vy_shift,circshift(W(:,:,2,it),[N/2,N/2]));
	xlabel('v_1');
	ylabel('v_2');
	title(sprintf('W_{x}(v), t=%g',(it-1)*dt));
	pause(0.05);
end
