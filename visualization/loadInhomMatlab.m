function loadInhomMatlab(datafile)
%LOADINHOMMATLAB - Load and visualize Matlab simulation results for the spatially inhomogeneous Boltzmann equation
%
%    LOADINHOMMATLAB(datafile)

% load simulation data from disk
load(datafile);

fprintf('N: %d\n',N);
fprintf('L: %g\n',L);
fprintf('J: %d\n',J);
fprintf('M: %d\n',M);
fprintf('R: %g\n',R);
fprintf('h: %g\n',h);	% spatial mesh width
fprintf('dt: %g\n',dt);	% time step
fprintf('number of time steps: %d\n',length(t));
fprintf('number of finite volumes: %d\n',length(x));
fprintf('boundary type: %s\n',btype);

% represent all data as a single array
W1 = zeros(N,N,4,length(x),length(t));
for it=1:length(t)
	for ix=1:length(x)
		for j=1:4
			for k2=1:N
				for k1=1:N
					W1(k1,k2,j,ix,it) = W{it,ix,j}(k1,k2);
				end
			end
		end
	end
end

visualizeInhomBoltzmann(N,L,J,M,R,h,dt,btype,W1,sprintf('%s\n%s boundary conditions',strrep(datafile,'_','-'),btype));
