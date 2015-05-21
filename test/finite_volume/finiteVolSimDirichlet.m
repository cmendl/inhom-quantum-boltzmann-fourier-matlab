function finiteVolSimDirichlet()
%FINITEVOLSIMDIRICHLET - Finite volume simulation example for u_t + A u_x = 0 with Dirichlet boundary conditions
%
%	FINITEVOLSIMDIRICHLET()

% eigenvalues of A
lambda = [-1,0.4,1.2];

% pseudo-random eigenvectors
R = [ 0.4133   -0.2381   -0.3634;
      0.2962   -0.1646    0.2212;
     -0.4013    0.1797   -0.3932];

% R*diag(lambda)*R^{-1}
A = R*diag(lambda)/R;

% spatial mesh width
h = 0.01;
% time step
dt = 0.005;

% Courant number, Eq. (10.56), must be <= 1
courant = max(abs(eig(A)))*dt/h;
fprintf('Courant number: %g (should be <= 1)\n',courant);

% spatial grid
xmax = 1;
x = 0:h:xmax;
% time grid
tmax = 0.4;
t = 0:dt:tmax;

% initial finite volume state
Uinit = zeros(length(A),length(x));
Uinit(1,:) = x+0.4;
Uinit(2,:) = cos(pi*x);
Uinit(3,ceil(end/2):end) = -1;
Uinit(3,ceil(end/2):end) = 1;	% jump

% run simulation
Ugod = cell(1,length(t));
Uslm = cell(1,length(t));
Ugod{1} = Uinit;
Uslm{1} = Uinit;
for n=2:length(t)
	Ugod{n} =  godunovStepDirichlet(h,dt,A,Ugod{n-1});
	Uslm{n} = slopeLimStepDirichlet(h,dt,A,Uslm{n-1});
end

% visualize simulation
figure();
for n=1:length(t)
	plot(x,Ugod{n}(1,:),'-.b',x,Ugod{n}(2,:),'-.g',x,Ugod{n}(3,:),'-.r',...
		 x,Uslm{n}(1,:),'b',  x,Uslm{n}(2,:),'g',  x,Uslm{n}(3,:),'r');
	xlabel('x');
	ylabel('U');
	title(sprintf('dashed: Godunov, continuous: slope limiter\nt = %g',t(n)));
	xlim([0,xmax]);
	ylim([-10,10]);
	pause(0.1);
end
