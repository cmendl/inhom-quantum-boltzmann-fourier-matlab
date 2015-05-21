function loadHomMatlab(datafile)
%LOADHOMMATLAB - Load and visualize Matlab simulation results for the spatially homogeneous Boltzmann equation
%
%    LOADHOMMATLAB(datafile)

% load simulation data from disk
load(datafile);

fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('M: %d\n',M);
fprintf('R: %g\n',R);
fprintf('dt: %g\n',dt);	% time step
fprintf('number of time steps: %d\n',length(t));

% represent all data as a single array
W1 = zeros(N,N,4,length(t));
for it=1:length(t)
	for j=1:4
		for k2=1:N
			for k1=1:N
				W1(k1,k2,j,it) = W{it,j}(k1,k2);
			end
		end
	end
end

visualizeHomBoltzmann(N,L,J,M,R,dt,W1);
